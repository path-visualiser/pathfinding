#ifndef CPD_SEARCH_H
#define CPD_SEARCH_H

// cpd_search.h
//
// A* implementation using an upper bound (CPD really) to enable:
//  - bounded sub-optimal search;
//  - heuristic weighting;
//  - anytime search;
//  - k-move search.
//
// @author: amaheo
// @created: 26/02/20
//

#include "cpool.h"
#include "pqueue.h"
#include "problem_instance.h"
#include "search.h"
#include "solution.h"
#include "timer.h"
#include "vec_io.h"
#include "log.h"
#include "dummy_listener.h"

#include <functional>
#include <iostream>
#include <memory>
#include <vector>

namespace warthog
{

// H is a heuristic function
// E is an expansion policy
// Q is the open list
// L is a "listener" which is used for callbacks
template< class H,
          class E,
          class Q = warthog::pqueue_min,
          class L = warthog::dummy_listener >
class cpd_search : public warthog::search
{
  public:
    cpd_search(H* heuristic, E* expander, Q* queue, L* listener = 0) :
        heuristic_(heuristic), expander_(expander), open_(queue),
        listener_(listener)
    {
        cost_cutoff_ = DBL_MAX;
        exp_cutoff_ = UINT32_MAX;
        time_cutoff_ = DBL_MAX;
        max_k_moves_ = UINT32_MAX;
        pi_.instance_id_ = UINT32_MAX;
        quality_cutoff_ = 1.0;
        // TODO Check whether this is the number of nodes
        k_moves_ = std::vector<uint32_t>(expander_->get_node_pool_size(), 0);
    }

    virtual ~cpd_search() { }

    virtual void
    get_pathcost(
        warthog::problem_instance& instance, warthog::solution& sol)
    {
        sol.reset();
        pi_ = instance;

        warthog::search_node* target = search(sol);
        if(target)
        {
            sol.sum_of_edge_costs_ = target->get_g();
        }
    }

    virtual void
    get_path(warthog::problem_instance& instance, warthog::solution& sol)
    {
        sol.reset();
        pi_ = instance;

        warthog::search_node* target = search(sol);
        if(target)
        {
            sol.sum_of_edge_costs_ = target->get_g();

            // follow backpointers to extract the path
            assert(expander_->is_target(target, &pi_));
            warthog::search_node* current = target;
            while(true)
            {
                sol.path_.push_back(current->get_id());
                if(current->get_parent() == warthog::SN_ID_MAX) break;
                current = expander_->generate(current->get_parent());
            }
            assert(sol.path_.back() == pi_.start_id_);

            DO_ON_DEBUG_IF(pi_.verbose_)
            {
                for(auto& node_id : sol.path_)
                {
                    int32_t x, y;
                    expander_->get_xy(node_id, x, y);
                    std::cerr
                            << "final path: (" << x << ", " << y << ")...";
                    warthog::search_node* n =
                            expander_->generate(node_id);
                    assert(n->get_search_number() == pi_.instance_id_);
                    n->print(std::cerr);
                    std::cerr << std::endl;
                }
            }
        }
    }

    // return a list of the nodes expanded during the last search
    // @param coll: an empty list
    void
    closed_list(std::vector<warthog::search_node*>& coll)
    {
        for(size_t i = 0; i < expander_->get_node_pool_size(); i++)
        {
            warthog::search_node* current = expander_->generate(i);
            if(current->get_search_number() == pi_.instance_id_)
            {
                coll.push_back(current);
            }
        }
    }

    // return a pointer to the warthog::search_node object associated
    // with node @param id. If this node was not generate during the
    // last search instance, 0 is returned instead
    warthog::search_node*
    get_generated_node(warthog::sn_id_t id)
    {
        warthog::search_node* ret = expander_->generate(id);
        return ret->get_search_number() == pi_.instance_id_ ? ret : 0;
    }

    // apply @param fn to every node on the closed list
    void
    apply_to_closed(std::function<void(warthog::search_node*)>& fn)
    {
        for(size_t i = 0; i < expander_->get_nodes_pool_size(); i++)
        {
            warthog::search_node* current = expander_->generate(i);
            if(current->get_search_number() == pi_.instance_id_)
            { fn(current); }
        }
    }

    // set a cost-cutoff to run a bounded-cost A* search.
    // the search terminates when the target is found or the f-cost
    // limit is reached.
    inline void
    set_cost_cutoff(warthog::cost_t cutoff) { cost_cutoff_ = cutoff; }

    inline warthog::cost_t
    get_cost_cutoff() { return cost_cutoff_; }

    // set a cutoff on the maximum number of node expansions.
    // the search terminates when the target is found or when
    // the limit is reached
    inline void
    set_max_expansions_cutoff(uint32_t cutoff) { exp_cutoff_ = cutoff; }

    inline uint32_t
    get_max_expansions_cutoff() { return exp_cutoff_; }

    // Set a time limit cutoff
    inline void
    set_max_time_cutoff(uint32_t cutoff) { time_cutoff_ = cutoff; }

    inline void
    set_max_us_cutoff(uint32_t cutoff) { set_max_time_cutoff(cutoff * 1e3); }

    inline void
    set_max_ms_cutoff(uint32_t cutoff) { set_max_time_cutoff(cutoff * 1e6); }

    inline void
    set_max_s_cutoff(uint32_t cutoff) { set_max_time_cutoff(cutoff * 1e9); }

    inline uint32_t
    get_max_time_cutoff() { return time_cutoff_; }

    // Set a k-radius cut-off -- stop expanding nodes further than k moves away
    // from the start.
    inline void
    set_max_k_moves(uint32_t k_moves) { max_k_moves_ = k_moves; }

    inline uint32_t
    get_max_k_moves() { return max_k_moves_; }

    // Set a quality cut-off, if the LB is within xx% of the UB we can stop
    inline void
    set_quality_cutoff(double cutoff) { quality_cutoff_ = cutoff; }

    inline double
    get_quality_cutoff() { return quality_cutoff_; }

    void
    set_listener(L* listener)
    { listener_ = listener; }

    virtual inline size_t
    mem()
    {
        size_t bytes =
                // memory for the priority quete
                open_->mem() +
                // gridmap size and other stuff needed to expand nodes
                expander_->mem() +
                // heuristic uses some memory too
                heuristic_->mem() +
                // misc
                sizeof(*this);
        return bytes;
    }


  private:
    H* heuristic_;
    E* expander_;
    Q* open_;
    L* listener_;
    warthog::problem_instance pi_;

    // early termination limits
    warthog::cost_t cost_cutoff_;   // Fixed upper bound
    uint32_t exp_cutoff_;           // Number of iterations
    double time_cutoff_;            // Time limit in nanoseconds
    uint32_t max_k_moves_;          // Max "distance" from target
    std::vector<uint32_t> k_moves_; // "Distance" from target
    double quality_cutoff_;

    // no copy ctor
    cpd_search(const cpd_search& other) { }
    cpd_search&
    operator=(const cpd_search& other) { return *this; }

    warthog::cost_t
    get_cost_(warthog::search_node *current, warthog::sn_id_t nid)
    {
        warthog::cost_t cost_to_n = 0;
        warthog::search_node *n;

        expander_->expand(current, &pi_);
        for(expander_->first(n, cost_to_n);
            n != nullptr;
            expander_->next(n, cost_to_n))
        {
            if (n->get_id() == nid)
            {
                return cost_to_n;
            }
        }

        error(pi_.verbose_, "Could not find", nid, "in neighbours of",
              current->get_id());
        return warthog::COST_MAX;
    }

    /**
     * Determine whether we should be pruning a node or adding it to the open
     * list.
     */
    bool
    should_prune_(warthog::search_node *incumbent,
                  warthog::search_node *n,
                  std::string stage)
    {
        bool prune = false;

        // if not [incumbent = nil or f(n) < f(incumbent)]
        if (incumbent != nullptr)
        {
            warthog::cost_t bound;

            // If we have an UB, we need to use it as source-of-truth
            if (incumbent->get_ub() < warthog::COST_MAX)
            {
                bound = incumbent->get_ub();

                if (n->get_f() * quality_cutoff_ > bound)
                {
                    debug(pi_.verbose_, stage, "by quality:", *n);
                    prune = true;
                }
            }
            else
            {
                bound = incumbent->get_f();
            }

            if (n->get_f() >= bound)
            {
                debug(pi_.verbose_, stage, "by f-val:", *n);
                prune = true;
            }
            // TODO Do we need to make a case where we have an incumbent's UB
            // but the current node does not?
            // else if (n->get_ub() < warthog::COST_MAX && n->get_ub() >= bound)
            // {
            //     debug(pi_.verbose_, stage, "by UB:", *n);
            //     prune = true;
            // }
        }

        if (k_moves_.at(n->get_id()) >= max_k_moves_)
        {
            debug(pi_.verbose_, stage, "by maximum k.");
            prune = true;
        }

        return prune;
    }

    bool
    early_stop_(warthog::search_node* current,
                warthog::solution* sol,
                warthog::timer* mytimer)
    {
        bool stop = false;

        mytimer->stop();
        // early termination: in case we want bounded-cost
        // search or if we want to impose some memory limit
        if(current->get_f() > cost_cutoff_)
        {
            trace(pi_.verbose_, "Cost cutoff", current->get_f(), ">",
                  cost_cutoff_);
            stop = true;
        }
        if(sol->nodes_expanded_ >= exp_cutoff_)
        {
            trace(pi_.verbose_, "Expanded cutoff", sol->nodes_expanded_, ">",
                  exp_cutoff_);
            stop = true;
        }
        // Exceeded time limit
        if (mytimer->elapsed_time_nano() > time_cutoff_)
        {
            trace(pi_.verbose_, "Time cutoff", mytimer->elapsed_time_nano(),
                  ">", time_cutoff_);
            stop = true;
        }
        // Extra early-stopping criteria when we have an upper bound; in
        // CPD terms, we have an "unperturbed path."
        if (current->get_f() == current->get_ub())
        {
            trace(pi_.verbose_, "Early stop");
            stop = true;
        }

        // A bit of a travestite use here.
        return stop;
    }

    /**
     * Initialise a new 'search_node' for the ongoing search given the parent
     * node (@param current).
     */
    void
    generate_node_(warthog::search_node* current,
                   warthog::search_node* n,
                   warthog::cost_t gval)
    {
        warthog::cost_t hval;
        warthog::cost_t ub;

        heuristic_->h(n->get_id(), pi_.target_id_, hval, ub);

        if (ub < warthog::COST_MAX)
        {
            ub += gval;
        }

        // Should we check for overflow here?
        //
        // KLUDGE We set the $g$ value of newly generated nodes to
        // infinity so we can pick them up when deciding whether to
        // prune.
        n->init(current->get_search_number(), current->get_id(),
                warthog::COST_MAX, gval + hval, ub);

        debug(pi_.verbose_, "Generating:", n->get_id());

        listener_->relax_node(n);
    }

    /**
     * Expand the current node: generate its neighbours, update the incumbent
     * and prune accordingly.
     */
    void
    expand_node_(warthog::search_node* current,
                 warthog::search_node* &incumbent,
                 warthog::solution* sol,
                 warthog::timer* mytimer)
    {
        warthog::cost_t cost_to_n = 0;
        uint32_t edge_id = 0;
        warthog::search_node* n;

        current->set_expanded(true); // NB: set before generating
        assert(current->get_expanded());
        sol->nodes_expanded_++;

        expander_->expand(current, &pi_);
        listener_->expand_node(current);

        // Incorrect timings reported otherwise
        DO_ON_DEBUG
        {
            mytimer->stop();

            info(pi_.verbose_, "[", mytimer->elapsed_time_micro(),"]",
                 sol->nodes_expanded_, "- Expanding:", current->get_id());
        }

        // The first loop over the neighbours is used to update the
        // incumbent.
        for(expander_->first(n, cost_to_n);
            n != nullptr;
            expander_->next(n, cost_to_n))
        {
            warthog::cost_t gval = current->get_g() + cost_to_n;
            sol->nodes_touched_++;

            listener_->generate_node(current, n, gval, edge_id);

            // Generate new search nodes
            if(n->get_search_number() != current->get_search_number())
            {
                generate_node_(current, n, gval);
            }

            // such that g(n) + c(n, n_i) + h(n_i) < f(incumbent)
            if (should_prune_(incumbent, n, "Ignore"))
            {
                continue;
            }

            // Incumbent handling

            // if n_i is a goal node
            if(expander_->is_target(n, &pi_))
            {
                incumbent = n;
                relax_node_(incumbent, gval, current->get_id());
                incumbent->set_ub(n->get_g());
                trace(pi_.verbose_, "New path to target:", n->get_g());
                debug(pi_.verbose_, *n);
            }
            else if (n->get_ub() < warthog::COST_MAX)
            {
                // Found a new incumbent
                if (incumbent == nullptr)
                {
                    trace(pi_.verbose_, "Found UB:", n->get_ub());
                    debug(pi_.verbose_, *n);
                    incumbent = n;
                }
                // Better incumbent
                else if (n->get_ub() < incumbent->get_ub())
                {
                    trace(pi_.verbose_, "Update UB:", n->get_ub());
                    debug(pi_.verbose_, *n);
                    incumbent = n;
                }
                // incumbent = n;
            }

            // Add nodes to queue

            // if n_i \in OPEN u CLOSED and g(n_i) > g(n) + c(n, n_i)
            if (gval < n->get_g())
            {
                relax_node_(n, gval, current->get_id());

                // g(n_i) <- g(n) + c(n, n_i)
                if(open_->contains(n))
                {
                    open_->decrease_key(n);
                    sol->nodes_updated_++;
                    debug(pi_.verbose_, "Updating:", *n);
                }
                // if n_i \in CLOSED
                else
                {
                    open_->push(n);
                    sol->nodes_inserted_++;
                    debug(pi_.verbose_, "Insert:", *n);
                }
            }
            else
            {
                debug(pi_.verbose_, "Skip:", *n);
            }
        }
    }

    /**
     * Differentiate between nodes that need relaxation and newly generated
     * nodes (cf. KLUDGE).
     */
    void
    relax_node_(warthog::search_node* n,
                warthog::cost_t gval,
                warthog::sn_id_t pid)
    {
        if (pid == warthog::SN_ID_MAX)
        {
            // Start node has no parent
            k_moves_.at(n->get_id()) = 0;
        }
        else
        {
            k_moves_.at(n->get_id()) = k_moves_.at(pid) + 1;
        }

        debug(pi_.verbose_, "Node", n->get_id(), "set to k=",
              k_moves_.at(n->get_id()));

        if (n->get_g() < warthog::COST_MAX)
        {
            n->relax(gval, pid);
            listener_->relax_node(n);
        }
        else
        {
            n->set_g(gval);
        }
    }

    // TODO refactor node information inside Stats
    warthog::search_node*
    search(warthog::solution& sol)
    {
        warthog::timer mytimer;
        mytimer.start();
        open_->clear();

        warthog::search_node* start;
        warthog::search_node* incumbent = nullptr;

        // get the internal target id
        if(pi_.target_id_ != warthog::SN_ID_MAX)
        {
            warthog::search_node* target =
                    expander_->generate_target_node(&pi_);
            if(!target) { return nullptr; } // invalid target location
            pi_.target_id_ = target->get_id();
        }

        // initialise and push the start node
        if(pi_.start_id_ == warthog::SN_ID_MAX) { return nullptr; }
        start = expander_->generate_start_node(&pi_);
        if(!start) { return nullptr; } // invalid start location
        pi_.start_id_ = start->get_id();

        // This heuristic also returns its upper bound
        warthog::cost_t start_h;
        warthog::cost_t start_ub;
        heuristic_->h(pi_.start_id_, pi_.target_id_, start_h, start_ub);

        // `hscale` is contained in the heuristic
        start->init(pi_.instance_id_, warthog::SN_ID_MAX, 0, start_h, start_ub);

        listener_->generate_node(0, start, 0, UINT32_MAX);

        user(pi_.verbose_, pi_);
        info(pi_.verbose_, "cut-off =", cost_cutoff_, "- tlim =", time_cutoff_,
             "- k-move =", max_k_moves_);
        debug(pi_.verbose_, "Start node:", *start);

        open_->push(start);
        sol.nodes_inserted_++;

        if (start_ub < warthog::COST_MAX)
        {
            // Having an UB means having a *concrete* path.
            incumbent = start;
            trace(pi_.verbose_, "Set UB:", incumbent->get_ub());
        }

        // begin expanding
        while(open_->size())
        {
            warthog::search_node* current = open_->pop();

            if (early_stop_(current, &sol, &mytimer))
                 // Stop if the $f$ value of UB of the best candidate node is
                 // worse than the incumbent.
                 // should_prune_(incumbent, current, "Stop"))
            {
                break;
            }

            // generate successors
            expand_node_(current, incumbent, &sol, &mytimer);
        }

        mytimer.stop();
        sol.time_elapsed_nano_ = mytimer.elapsed_time_nano();
        sol.nodes_surplus_ = open_->size();

        DO_ON_DEBUG_IF(pi_.verbose_)
        {
            if(incumbent == nullptr)
            {
                warning(pi_.verbose_, "Search failed; no solution exists.");
            }
            else
            {
                user(pi_.verbose_, "Best incumbent", *incumbent);
            }
        }

        // Rebuild path from incumbent to solution
        while (incumbent != nullptr && !expander_->is_target(incumbent, &pi_))
        {
            warthog::sn_id_t pid = incumbent->get_id();
            warthog::sn_id_t nid = heuristic_->get_move(pid, pi_.target_id_);

            if (nid == warthog::SN_ID_MAX)
            {
                warning(pi_.verbose_, "Cannot rebuild path from", pid);
                incumbent = nullptr;
            }
            else
            {
                warthog::search_node* n = expander_->generate(nid);
                warthog::cost_t gval =
                        incumbent->get_g() + get_cost_(incumbent, nid);

                // Force node "generation" in case some of them were not
                // generated during the search.
                generate_node_(incumbent, n, gval);
                n->set_g(gval);  // Get rid of the silly $g=\inf$
                incumbent = n;
            }
        }

        return incumbent;
    }
};

}

#endif
