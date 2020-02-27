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

#include <functional>
#include <iostream>
#include <memory>
#include <vector>

namespace warthog
{

// H is a heuristic function
// E is an expansion policy
// S is a stats class
template< class H,
          class E,
          class Q = warthog::pqueue_min >
class cpd_search : public warthog::search
{
  public:
    cpd_search(H* heuristic, E* expander, Q* queue) :
            heuristic_(heuristic), expander_(expander)
    {
        open_ = queue;
        cost_cutoff_ = DBL_MAX;
        exp_cutoff_ = UINT32_MAX;
        time_lim_ = DBL_MAX;
        max_k_move_ = UINT32_MAX;
        on_relax_fn_ = 0;
        on_generate_fn_ = 0;
        on_expand_fn_ = 0;
        pi_.instance_id_ = UINT32_MAX;
    }

    virtual ~cpd_search() { }

    virtual void
    get_distance(
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
        for(size_t i = 0; i < expander_->get_nodes_pool_size(); i++)
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

    // apply @param fn every time a node is successfully relaxed
    void
    apply_on_relax(std::function<void(warthog::search_node*)>& fn)
    {
        on_relax_fn_ = &fn;
    }

    // apply @param fn every time a node is generated (equiv, reached)
    void
    apply_on_generate(
        std::function<void(
            warthog::search_node* succ,
            warthog::search_node* from,
            warthog::cost_t edge_cost,
            uint32_t edge_id)>& fn)
    {
        on_generate_fn_ = &fn;
    }

    // apply @param fn when a node is popped off the open list for
    // expansion
    void
    apply_on_expand(std::function<void(warthog::search_node*)>& fn)
    {
        on_expand_fn_ = &fn;
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
    warthog::problem_instance pi_;

    // early termination limits
    warthog::cost_t cost_cutoff_;  // Fixed upper bound
    warthog::cost_t upper_bound_;  // Dynamic upper bound (based on h-value)
    uint32_t exp_cutoff_;          // Number of iterations
    double time_lim_;              // Time limit in nanoseconds
    uint32_t max_k_move_;          // "Distance" from target

    // callback for when a node is relaxed
    std::function<void(warthog::search_node*)>* on_relax_fn_;

    // callback for when a node is reached / generated
    std::function<void(
        warthog::search_node*,
        warthog::search_node*,
        warthog::cost_t edge_cost,
        uint32_t edge_id)>* on_generate_fn_;

    // callback for when a node is expanded
    std::function<void(warthog::search_node*)>* on_expand_fn_;

    // no copy ctor
    cpd_search(const cpd_search& other) { }
    cpd_search&
    operator=(const cpd_search& other) { return *this; }

    warthog::search_node* search(warthog::solution& sol);
};

}

#endif
