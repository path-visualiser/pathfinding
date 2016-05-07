#ifndef WARTHOG_BIDIRECTIONAL_SEARCH_H
#define WARTHOG_BIDIRECTIONAL_SEARCH_H

// bidirectional_search.h
//
// A customisable variant of bidirectional best-first search.
// Users can pass in any heuristic and any (domain-specific) expansion policy.
//
// @author: dharabor
// @created: 2016-02-14
//

#include "constants.h"
#include "graph_expansion_policy.h"
#include "planar_graph.h"
#include "pqueue.h"
#include "search.h"
#include "search_node.h"
#include "timer.h"
#include "zero_heuristic.h"

#include "constants.h"
#include <cstdlib>
#include <stack>
#include <stdint.h>

namespace warthog
{

namespace graph
{
    class planar_graph;
}

class expansion_policy;
class pqueue;
class search_node;

typedef double (* heuristicFn)
(uint32_t nodeid, uint32_t targetid);

template<class H>
class bidirectional_search : public warthog::search
{
    public:
        bidirectional_search(
                warthog::graph::planar_graph* fwd_g,
                warthog::graph::planar_graph* bwd_g,
                H* heuristic) 
            : heuristic_(heuristic)
        {
            verbose_ = false;
            fexpander_ = new warthog::graph_expansion_policy(fwd_g);
            bexpander_ = new warthog::graph_expansion_policy(bwd_g);
            fopen_ = new pqueue(512, true);
            bopen_ = new pqueue(512, true);
            
            dijkstra_ = false;
            if(typeid(*heuristic_) == typeid(warthog::zero_heuristic))
            {
                dijkstra_ = true;
            }
        }

        ~bidirectional_search()
        {
            delete fopen_;
            delete bopen_;
            delete fexpander_;
            delete bexpander_;
        }

        double 
        get_length(uint32_t startid, uint32_t goalid)
        {
            this->search(startid, goalid);

#ifndef NDEBUG

#endif
            cleanup();
            return best_cost_ / (double)warthog::ONE;
        }
            
		inline bool
		get_verbose() { return verbose_; }

		inline void
		set_verbose(bool verbose) { verbose_ = verbose; } 

        size_t
        mem()
        {
            return sizeof(*this) + 
                fopen_->mem() +
                bopen_->mem() +
                fexpander_->mem();
                bexpander_->mem();
        }

    private:
        warthog::pqueue* fopen_;
        warthog::pqueue* bopen_;
        warthog::expansion_policy* fexpander_;
        warthog::expansion_policy* bexpander_;
        H* heuristic_;
        bool dijkstra_;

        // some temporary variables used during search
        warthog::search_node* v_;
        warthog::search_node* w_;
        warthog::cost_t best_cost_;
        warthog::problem_instance instance_;

        void 
        search(uint32_t startid, uint32_t goalid)
        {
            warthog::timer mytimer;
            mytimer.start();

            // init
            this->reset_metrics();
            best_cost_ = warthog::INF;
            v_ = w_ = 0;

            #ifndef NDEBUG
            if(verbose_)
            {
                std::cerr << "bidirectional_search: startid="
                    << startid<<" goalid=" <<goalid
                    << std::endl;
            }
            #endif

            instance_.set_goal(goalid);
            instance_.set_start(startid);
            instance_.set_searchid(++warthog::search::searchid_);

            warthog::search_node* start = fexpander_->generate(startid);
            start->reset(instance_.get_searchid());
            start->set_g(0);
            start->set_f(heuristic_->h(startid, goalid));
            fopen_->push(start);

            warthog::search_node* goal = bexpander_->generate(goalid);
            goal->reset(instance_.get_searchid());
            goal->set_g(0);
            goal->set_f(heuristic_->h(startid, goalid));
            bopen_->push(goal);

            // these variables can be used to balance the backward and
            // forward expansion. set bwd_modulo higher to expand nodes
            // less often in the backward direction
            int fwd_expansions = 0;
            int bwd_modulo = 1;

            // interleave search by expanding one node in the forward 
            // direction followed by a node in the backward direction
            while(fopen_->size() > 0 && bopen_->size() > 0)
            {
                // expand in the forward direction
                warthog::search_node* current = fopen_->pop();
                if(current->get_f() >= best_cost_)
                {
                    // early terminate if we cannot improve the best path
                    break;
                }
                expand(current, fopen_, fexpander_, bexpander_, goalid);
                //fwd_expansions++;

                if(dijkstra_ || !(fwd_expansions % bwd_modulo))
                {
                    // expand in the backward direction
                    current = bopen_->pop();
                    if(current->get_f() >= best_cost_)
                    {
                        // early terminate if we cannot improve the best path
                        break;
                    }
                    expand(current, bopen_, bexpander_, fexpander_, startid);
                }

                // special case for dijkstra; terminate when fringes meet
                if(dijkstra_ && best_cost_ < warthog::INF)
                {
                    break;
                }
            }

            assert(best_cost_ != warthog::INF ||
                    (v_->get_g() != warthog::INF && w_->get_g() != warthog::INF));

			mytimer.stop();
			search_time_ = mytimer.elapsed_time_micro();
        }

        void
        expand( warthog::search_node* current,
                warthog::pqueue* open,
                warthog::expansion_policy* expander,
                warthog::expansion_policy* reverse_expander, 
                uint32_t tmp_goalid)
        {
            // goal test
            if(current->get_id() == tmp_goalid) 
            {
                best_cost_ = current->get_g();
                v_ = current;
                w_ = 0;
                return;
            }

            // goal not found yet; expand as normal
            current->set_expanded(true);
            expander->expand(current, &instance_);
            nodes_expanded_++;

            #ifndef NDEBUG
            if(verbose_)
            {
                int32_t x, y;
                expander->get_xy(current, x, y);
                std::cerr << this->nodes_expanded_ << ". "
                    "expanding " << (instance_.get_goal() == tmp_goalid ? "(f)" : "(b)")
                    << " ("<<x<<", "<<y<<")...";
                current->print(std::cerr);
                std::cerr << std::endl;
            }
            #endif
            
            // prune the search if current is already expanded in the
            // reverse direction
            warthog::search_node* reverse_current = 
                reverse_expander->generate(current->get_id());
            if(reverse_current->get_searchid() == current->get_searchid() &&
               reverse_current->get_expanded())
            {
                // relax the best path if necessary
                if((current->get_g() + reverse_current->get_g()) < best_cost_)
                {
                    v_ = current;
                    w_ = reverse_current;
                    best_cost_ = current->get_g() + reverse_current->get_g();

                    #ifndef NDEBUG
                    if(verbose_)
                    {
                        int32_t x, y;
                        expander->get_xy(current, x, y);
                        std::cerr <<"new best solution!  cost=" << best_cost_<<std::endl;
                    }
                    #endif
                }
                #ifndef NDEBUG
                if(verbose_)
                {
                    int32_t x, y;
                    expander->get_xy(current, x, y);
                    std::cerr <<"closing ("<<x<<", "<<y<<")...";
                    current->print(std::cerr);
                    std::cerr << std::endl;
                }
                #endif
                return;
            }


            // generate all neighbours
            warthog::search_node* n = 0;
            warthog::cost_t cost_to_n = warthog::INF;
            for(expander->first(n, cost_to_n); n != 0; expander->next(n, cost_to_n))
            {
                nodes_touched_++;
                if(n->get_expanded())
                {
                    // skip neighbours already expanded
                    #ifndef NDEBUG
                    if(verbose_)
                    {
                        int32_t x, y;
                        expander->get_xy(n, x, y);
                        std::cerr << "  closed; (edgecost=" << cost_to_n << ") "
                            << "("<<x<<", "<<y<<")...";
                        n->print(std::cerr);
                        std::cerr << std::endl;

                    }
                    #endif
                    continue;
                }

                // relax (or generate) each neighbour
                warthog::cost_t gval = current->get_g() + cost_to_n;
                if(open->contains(n))
                {
                    // update a node from the fringe
                    if(gval < n->get_g())
                    {
                        n->relax(gval, current);
                        open->decrease_key(n);
                        #ifndef NDEBUG
                        if(verbose_)
                        {
                            int32_t x, y;
                            expander->get_xy(n, x, y);
                            std::cerr << "  open; updating "
                                << "(edgecost="<< cost_to_n<<") "
                                << "("<<x<<", "<<y<<")...";
                            n->print(std::cerr);
                            std::cerr << std::endl;
                        }
                        #endif
                    }
                    else
                    {
                        #ifndef NDEBUG
                        if(verbose_)
                        {
                            int32_t x, y;
                            expander->get_xy(n, x, y);
                            std::cerr << "  open; not updating "
                                << "(edgecost=" << cost_to_n<< ") "
                                << "("<<x<<", "<<y<<")...";
                            n->print(std::cerr);
                            std::cerr << std::endl;
                        }
                        #endif
                    }
                }
                else
                {
                    // add a new node to the fringe
                    n->set_g(gval);
                    n->set_f(gval + heuristic_->h(n->get_id(), tmp_goalid));
                    n->set_parent(current);
                    n->set_searchid(current->get_searchid());
                    open->push(n);
                    #ifndef NDEBUG
                    if(verbose_)
                    {
                        int32_t x, y;
                        expander->get_xy(n, x, y);
                        std::cerr << "  generating "
                            << "(edgecost=" << cost_to_n<<") " 
                            << "("<<x<<", "<<y<<")...";
                        n->print(std::cerr);
                        std::cerr << std::endl;
                    }
                    #endif
                    nodes_generated_++;
                }
            }

            #ifndef NDEBUG
            if(verbose_)
            {
                int32_t x, y;
                expander->get_xy(current, x, y);
                std::cerr <<"closing ("<<x<<", "<<y<<")...";
                current->print(std::cerr);
                std::cerr << std::endl;
            }
            #endif
        }

        void
        cleanup()
        {
            fopen_->clear();
            bopen_->clear();
            fexpander_->clear();
            bexpander_->clear();
        }

};

}

#endif

