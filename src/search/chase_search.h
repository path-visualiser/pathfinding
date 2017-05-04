#ifndef WARTHOG_CHASE_SEARCH_H
#define WARTHOG_CHASE_SEARCH_H

// search/chase_search.h
//
// An algorithm resembling CHASE. For theoretical details see:
//
// [Bauer, Delling, Sanders, Schieferdecker, Schultes and Wagner, 
// Combining Hierarchical and Goal-directed Speed-up Techniques 
// for Dijkstra's Algorithm, Journal of Experimental Algorithms,
// vol 15, 2010]
//
// NB: NOT FULLY IMPLEMENTED!! Only stalls nodes higher than 
// some cutoff and resumes the search from these nodes if no
// optimal path has been found.
//
// @author: dharabor
// @created: 2016-09-10
//

#include "constants.h"
#include "ch_expansion_policy.h"
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
class chase_search : public warthog::search
{
    public:
        chase_search(
                warthog::ch_expansion_policy* fexp,
                warthog::ch_expansion_policy* bexp,
                H* heuristic) 
            : fexpander_(fexp), bexpander_(bexp), heuristic_(heuristic)
        {
            fopen_ = new pqueue(512, true);
            bopen_ = new pqueue(512, true);
            
            dijkstra_ = false;
            if(typeid(*heuristic_) == typeid(warthog::zero_heuristic))
            {
                dijkstra_ = true;
            }
            
            max_phase1_rank_ = fexpander_->get_num_nodes()*0.95;
        }

        ~chase_search()
        {
            delete fopen_;
            delete bopen_;
        }

        void
        get_path(warthog::problem_instance& pi, warthog::solution& sol)
        {
            pi_ = pi;
            this->search(sol);

            sol.sum_of_edge_costs_ = best_cost_;
            if(best_cost_ != warthog::INF) { reconstruct_path(sol.path_); }
            cleanup();

#ifndef NDEBUG
            if(pi_.verbose_)
            {
                warthog::search_node* apex = v_;
                for(uint32_t i = 0; i < sol.path_.size(); i++)
                {
                    warthog::search_node* n = 
                        fexpander_->get_ptr(
                            sol.path_.at(i), pi_.instance_id_);
                    int32_t x, y;
                    fexpander_->get_xy(n->get_id(), x, y);
                    if(&*n == &*apex)
                    {
                        std::cerr << "(apex)";
                    }
                    std::cerr 
                        << "final path: (" << x << ", " << y << ")...";
                    n->print(std::cerr);
                    std::cerr << std::endl;
                }
            }
#endif
        }
            
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
        warthog::ch_expansion_policy* fexpander_;
        warthog::ch_expansion_policy* bexpander_;
        H* heuristic_;
        bool dijkstra_;
        bool forward_next_;

        // CHASE-specific stuff
        uint32_t phase_;
        uint32_t max_phase1_rank_;
        std::vector<warthog::search_node*> fwd_norelax;
        std::vector<warthog::search_node*> bwd_norelax;

        // v is the section of the path in the forward
        // direction and w is the section of the path
        // in the backward direction. need parent pointers
        // of both to extract the actual path
        warthog::search_node* v_;
        warthog::search_node* w_;
        double best_cost_;
        warthog::problem_instance pi_;

        void
        reconstruct_path(std::vector<uint32_t>& path)
        {
            if(v_ && (&*v_ == &*bexpander_->generate(v_->get_id())))
            {
                warthog::search_node* tmp = v_;
                v_ = w_;
                w_ = tmp;
            }

            warthog::search_node* current = w_;
            while(current)
            {  
               path.push_back(current->get_id());
               current = current->get_parent();
            }

            current = v_->get_parent();
            while(current)
            {
               path.push_back(current->get_id());
               current = current->get_parent();
            }
        }

        // modify this function to balance the search
        // by default the search expands one node in
        // each direction then switches to the other direction
        bool
        forward_next()
        {
            forward_next_ = !forward_next_;
            return forward_next_;
        }

        void 
        search(warthog::solution& sol)
        {
            warthog::timer mytimer;
            mytimer.start();

            // init
            v_ = w_ = 0;
            fwd_norelax.clear();
            bwd_norelax.clear();
            forward_next_ = true;

            #ifndef NDEBUG
            if(pi_.verbose_)
            {
                std::cerr << "chase_search. ";
                pi_.print(std::cerr);
                std::cerr << std::endl;
            }
            #endif

            warthog::search_node *start, *goal;
            start = fexpander_->generate_start_node(&pi_);
            goal = bexpander_->generate_target_node(&pi_);

            start->init(pi_.instance_id_, 0, 0, 
                    heuristic_->h(pi_.start_id_, pi_.target_id_));
            fopen_->push(start);

            goal->init(pi_.instance_id_, 0, 0, 
                    heuristic_->h(pi_.start_id_, pi_.target_id_));
            bopen_->push(goal);

            // interleave search; we terminate when the best lower bound
            // is larger than the cost of the best solution so far (or when
            // we exhaust both open lists)
            //
            // NB: in a standard formulation we can terminate when either
            // open list is exhausted. however in a more general context
            // bi-directional search might choose to ignore an edge so it
            // can be explored in the opposite direction. thus one open list
            // can be exhausted and it still makes sense to continue the
            // search. such an algorithm might be more efficient
            // (e.g. CH with arc flags) 
            warthog::search_node *ftop = 0, *btop = 0;
            phase_ = 1;
            bool cannot_improve = false;
            while(true)
            {
                if(cannot_improve || !(fopen_->size() || bopen_->size()))
                {
                    if(phase_ == 1)
                    {
                        for(uint32_t i = 0; i < fwd_norelax.size(); i++)
                        {
                            fopen_->push(fwd_norelax.at(i));
                        }
                        for(uint32_t i = 0; i < bwd_norelax.size(); i++)
                        {
                            bopen_->push(bwd_norelax.at(i));
                        }
                        phase_++;
                        cannot_improve = false;
#ifndef NDEBUG
                    if(pi_.verbose_)
                    {
                        std::cerr << "=== PHASE2 ===" << std::endl;
                    }
#endif

                    }
                    else { break; }  // no solution
                }

                if(fopen_->size() > 0) { ftop = fopen_->peek(); } 
                if(bopen_->size() > 0) { btop = bopen_->peek(); } 

                // the way we calculate the lower-bound on solution cost 
                // differs when we have a heuristic available vs not.
                uint32_t best_bound_ = dijkstra_ ? 
                    // no heuristic
                    (ftop->get_g() + btop->get_g()) : 
                    // with a heuristic
                    (std::min( 
                        ftop->get_f(), btop->get_f()));

                // terminate if we cannot improve the best solution found so far
                if(best_bound_ > best_cost_)
                {
#ifndef NDEBUG
                    if(pi_.verbose_)
                    {
                        std::cerr << "provably-best solution found; cost=" << 
                            best_cost_ << std::endl;
                    }
#endif
                    break;
                }

                // ok, we still have hope. let's keep expanding. 
                if(forward_next())
                {
                    if(ftop->get_f() < best_cost_ && fopen_->size()) 
                    {
                        warthog::search_node* current = fopen_->pop();
                        expand(current, fopen_, fexpander_, bexpander_, 
                                pi_.target_id_, fwd_norelax, sol);
                    }
                    else if(btop->get_f() < best_cost_ && bopen_->size())
                    {
                        // we can't improve the best solution in the forward
                        // direction so we ignore the interleave policy and 
                        // expand a node in the backward direction instead
                        warthog::search_node* current = bopen_->pop();
                        expand(current, bopen_, bexpander_, fexpander_, 
                                pi_.start_id_, bwd_norelax, sol);
                    }
                    else 
                    {
                        cannot_improve = true;
                        // search is finished in both directions
                    }
                }
                else 
                {
                    if(btop->get_f() < best_cost_ && bopen_->size())
                    {
                        warthog::search_node* current = bopen_->pop();
                        expand(current, bopen_, bexpander_, fexpander_, 
                                pi_.start_id_, bwd_norelax, sol);
                    }
                    else if(ftop->get_f() < best_cost_ && fopen_->size())
                    {
                        // we can't improve the best solution in the backward
                        // direction so we ignore the interleave policy and 
                        // expand a node in the forward direction instead
                        warthog::search_node* current = fopen_->pop();
                        expand(current, fopen_, fexpander_, bexpander_, 
                                pi_.target_id_, fwd_norelax, sol);
                    }
                    else
                    {
                        // search is finished in both directions
                        cannot_improve = true;
                    }
                }
            }

            assert(best_cost_ != warthog::INF || (v_ == 0 && w_ == 0));

			mytimer.stop();
			sol.time_elapsed_micro_= mytimer.elapsed_time_micro();
        }

        void
        expand( warthog::search_node* current,
                warthog::pqueue* open,
                warthog::ch_expansion_policy* expander,
                warthog::ch_expansion_policy* reverse_expander, 
                uint32_t tmp_goalid,
                std::vector<warthog::search_node*>& norelax, 
                warthog::solution& sol)
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
            expander->expand(current, &pi_);
            sol.nodes_expanded_++;

            #ifndef NDEBUG
            if(pi_.verbose_)
            {
                int32_t x, y;
                expander->get_xy(current->get_id(), x, y);
                std::cerr 
                    << sol.nodes_expanded_ 
                    << ". expanding " 
                    << (pi_.target_id_ == tmp_goalid ? "(f)" : "(b)")
                    << " ("<<x<<", "<<y<<")...";
                current->print(std::cerr);
                std::cerr << std::endl;
            }
            #endif
            
            // generate all neighbours
            warthog::search_node* n = 0;
            double cost_to_n = warthog::INF;
            for(expander->first(n, cost_to_n); n != 0; expander->next(n, cost_to_n))
            {
                sol.nodes_touched_++;
                if(n->get_expanded())
                {
                    // skip neighbours already expanded
                    #ifndef NDEBUG
                    if(pi_.verbose_)
                    {
                        int32_t x, y;
                        expander->get_xy(n->get_id(), x, y);
                        std::cerr << "  closed; (edgecost=" << cost_to_n << ") "
                            << "("<<x<<", "<<y<<")...";
                        n->print(std::cerr);
                        std::cerr << std::endl;

                    }
                    #endif
                    continue;
                }

                // relax (or generate) each neighbour
                double gval = current->get_g() + cost_to_n;
                if(open->contains(n))
                {
                    // update a node from the fringe
                    if(gval < n->get_g())
                    {
                        sol.nodes_updated_++;
                        n->relax(gval, current);
                        open->decrease_key(n);
                        #ifndef NDEBUG
                        if(pi_.verbose_)
                        {
                            int32_t x, y;
                            expander->get_xy(n->get_id(), x, y);
                            std::cerr << " updating "
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
                        if(pi_.verbose_)
                        {
                            int32_t x, y;
                            expander->get_xy(n->get_id(), x, y);
                            std::cerr << " not updating "
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
                    if(phase_ == 1 && 
                            n->get_search_id() == current->get_search_id())
                    {
                        // relax the g-value of the nodes not being
                        // expanded in phase1
                        // (these are not added to open yet)
                        if(gval < n->get_g())
                        {
                            n->relax(gval, current);
                            #ifndef NDEBUG
                            if(pi_.verbose_)
                            {
                                int32_t x, y;
                                expander->get_xy(n->get_id(), x, y);
                                std::cerr << " phase1 norelax update "
                                    << "(edgecost="<< cost_to_n<<") "
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
                        sol.nodes_inserted_++;
                        n->init(current->get_search_id(), 
                                current, 
                                gval,
                                gval + heuristic_->h(n->get_id(), tmp_goalid));

                        if( phase_ == 2 || 
                            expander->get_rank(n->get_id()) < max_phase1_rank_)
                        {
                            open->push(n);
                        }
                        else
                        {
                            norelax.push_back(n);
                            #ifndef NDEBUG
                            if(pi_.verbose_)
                            {
                                std::cerr << " phase1 norelax;";
                            }
                            #endif
                        }
                    }

                    #ifndef NDEBUG
                    if(pi_.verbose_)
                    {
                        int32_t x, y;
                        expander->get_xy(n->get_id(), x, y);
                        std::cerr << " generating "
                            << "(edgecost=" << cost_to_n<<") " 
                            << "("<<x<<", "<<y<<")...";
                        n->print(std::cerr);
                        std::cerr << std::endl;
                    }
                    #endif
                }

                // update the best solution if possible
                warthog::search_node* reverse_n = 
                    reverse_expander->generate(n->get_id());
                if(reverse_n->get_search_id() == n->get_search_id())
//                        && reverse_n->get_expanded())
                {
                    if((current->get_g() + cost_to_n + reverse_n->get_g()) < best_cost_)
                    {
                        v_ = current;
                        w_ = reverse_n;
                        best_cost_ = current->get_g() + cost_to_n + reverse_n->get_g();

                        #ifndef NDEBUG
                        if(pi_.verbose_)
                        {
                            int32_t x, y;
                            expander->get_xy(current->get_id(), x, y);
                            std::cerr <<"new best solution!  cost=" << best_cost_<<std::endl;
                        }
                        #endif
                    }
                }
            }

            #ifndef NDEBUG
            if(pi_.verbose_)
            {
                int32_t x, y;
                expander->get_xy(current->get_id(), x, y);
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

