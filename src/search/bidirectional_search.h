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
                warthog::expansion_policy* fexp,
                warthog::expansion_policy* bexp,
                H* heuristic) 
            : fexpander_(fexp), bexpander_(bexp), heuristic_(heuristic)
        {
            verbose_ = false;
            fopen_ = new pqueue(512, true);
            bopen_ = new pqueue(512, true);
            
            dijkstra_ = false;
            forward_next_ = true;
            if(typeid(*heuristic_) == typeid(warthog::zero_heuristic))
            {
                dijkstra_ = true;
            }
        }

        ~bidirectional_search()
        {
            delete fopen_;
            delete bopen_;
        }

        double 
        get_length(warthog::problem_instance pi)
        {
            this->search(pi);

#ifndef NDEBUG
            if(get_verbose())
            {
                std::vector<warthog::search_node*> path;
                warthog::search_node* apex = v_;
                reconstruct_path(path);

                for(uint32_t i = 0; i < path.size(); i++)
                {
                    warthog::search_node* n = path.at(i);
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
            cleanup();
            return best_cost_;
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
        bool forward_next_;

        // v is the section of the path in the forward
        // direction and w is the section of the path
        // in the backward direction. need parent pointers
        // of both to extract the actual path
        warthog::search_node* v_;
        warthog::search_node* w_;
        double best_cost_;
        warthog::problem_instance pi_;

        void
        reconstruct_path(std::vector<warthog::search_node*>& path)
        {
            if(v_ && (&*v_ == &*bexpander_->generate(v_->get_id())))
            {
                warthog::search_node* tmp = v_;
                v_ = w_;
                w_ = tmp;
            }

            while(v_)
            {
               path.push_back(v_);
               v_ = v_->get_parent();
            }
            std::reverse(path.begin(), path.end());

            while(w_)
            {
               path.push_back(w_);
               w_ = w_->get_parent();
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
        search(warthog::problem_instance pi)
        {
            warthog::timer mytimer;
            mytimer.start();

            // init
            this->reset_metrics();
            best_cost_ = warthog::INF;
            v_ = w_ = 0;
            forward_next_ = false;

            pi_ = pi;
            pi_.set_search_id(++warthog::search::searchid_);

            #ifndef NDEBUG
            if(verbose_)
            {
                std::cerr << "bidirectional_search. ";
                pi_.print(std::cerr);
                std::cerr << std::endl;
            }
            #endif


            // generate the start and target nodes. also, we update the 
            // problem instance with internal ids (makes debugging easier)
            warthog::search_node* start = fexpander_->generate_start_node(&pi_);
            warthog::search_node* target = bexpander_->generate_target_node(&pi_);
            pi_.set_start_id(start->get_id());
            pi_.set_target_id(target->get_id());

            // compute an initial priority for both the start and target
            start->init(pi_.get_search_id(), 0, 0, 
                        heuristic_->h(start->get_id(), target->get_id()));
            target->init(pi_.get_search_id(), 0, 0, 
                         heuristic_->h(start->get_id(), target->get_id()));
            fopen_->push(start);
            bopen_->push(target);

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
            while(fopen_->size() || bopen_->size())
            {
                if(fopen_->size() > 0) { ftop = fopen_->peek(); } 
                if(bopen_->size() > 0) { btop = bopen_->peek(); } 

                // The code below can terminate bi-dijkstra search 
                // more quickly by computing the bound as the sum of g-values
                // at the top of the open list in both directions. This works
                // for UNDIRECTED GRAPHS ONLY
                //uint32_t best_bound = dijkstra_ ? 
                //     //no heuristic
                //    (ftop->get_g() + btop->get_g()) : 
                //     //with a heuristic
                //    (std::min( 
                //        ftop->get_f(), btop->get_f()));

                // the code below computes the bound more simply and works
                // on directed planar graphs with asymmetric costs
                uint32_t best_bound = std::min( 
                    ftop->get_f(), btop->get_f());

                // terminate if we cannot improve the best solution so far.
                // NB: bidirectional dijkstra stops when the two search 
                // frontiers expand the same node; bidirectional A* stops 
                // when the best bound is larger than the best solution
                if(best_bound > best_cost_ || 
                        (dijkstra_ && v_ && w_ && 
                         v_->get_expanded() && w_->get_expanded()))
                {
#ifndef NDEBUG
                    if(verbose_)
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
                                pi_.get_target_id());
                    }
                    else if(btop->get_f() < best_cost_ && bopen_->size())
                    {
                        // we can't improve the best solution in the forward
                        // direction so we ignore the interleave policy and 
                        // expand a node in the backward direction instead
                        warthog::search_node* current = bopen_->pop();
                        expand(current, bopen_, bexpander_, fexpander_, 
                                pi_.get_start_id());
                    }
                    else 
                    {
                        // search is finished in both directions
                        break;
                    }
                }
                else 
                {
                    if(btop->get_f() < best_cost_ && bopen_->size())
                    {
                        warthog::search_node* current = bopen_->pop();
                        expand(current, bopen_, bexpander_, fexpander_, 
                                pi_.get_start_id());
                    }
                    else if(ftop->get_f() < best_cost_ && fopen_->size())
                    {
                        // we can't improve the best solution in the backward
                        // direction so we ignore the interleave policy and 
                        // expand a node in the forward direction instead
                        warthog::search_node* current = fopen_->pop();
                        expand(current, fopen_, fexpander_, bexpander_, 
                                pi_.get_target_id());
                    }
                    else
                    {
                        // search is finished in both directions
                        break;
                    }
                }
            }

            assert(best_cost_ != warthog::INF || (v_ == 0 && w_ == 0));

			mytimer.stop();
			search_time_ = mytimer.elapsed_time_micro();
        }

        void
        expand( warthog::search_node* current,
                warthog::pqueue* open,
                warthog::expansion_policy* expander,
                warthog::expansion_policy* reverse_expander, 
                uint32_t tmp_targetid)
        {
            // target test
            if(current->get_id() == tmp_targetid) 
            {
                best_cost_ = current->get_g();
                v_ = current;
                w_ = 0;
                return;
            }

            // target not found yet; expand as normal
            current->set_expanded(true);
            expander->expand(current, &pi_);
            nodes_expanded_++;

            #ifndef NDEBUG
            if(verbose_)
            {
                int32_t x, y;
                expander->get_xy(current->get_id(), x, y);
                std::cerr << this->nodes_expanded_ << ". "
                    "expanding " << (pi_.get_target_id() == tmp_targetid ? "(f)" : "(b)")
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
                nodes_touched_++;
                if(n->get_expanded())
                {
                    // skip neighbours already expanded
                    #ifndef NDEBUG
                    if(verbose_)
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
                        n->relax(gval, current);
                        open->decrease_key(n);
                        #ifndef NDEBUG
                        if(verbose_)
                        {
                            int32_t x, y;
                            expander->get_xy(n->get_id(), x, y);
                            std::cerr << "  open; updating "
                                << "(edgecost="<< cost_to_n<<") "
                                << "("<<x<<", "<<y<<")...";
                            n->print(std::cerr);
                            std::cerr << std::endl;
                        }
                        #endif
                    }
                    #ifndef NDEBUG
                    else
                    {
                        if(verbose_)
                        {
                            int32_t x, y;
                            expander->get_xy(n->get_id(), x, y);
                            std::cerr << "  open; not updating "
                                << "(edgecost=" << cost_to_n<< ") "
                                << "("<<x<<", "<<y<<")...";
                            n->print(std::cerr);
                            std::cerr << std::endl;
                        }
                    }
                    #endif
                }
                else
                {
                    // add a new node to the fringe
                    n->init(current->get_search_id(), 
                            current, 
                            gval,
                            gval + heuristic_->h(n->get_id(), tmp_targetid));
                    open->push(n);
                    #ifndef NDEBUG
                    if(verbose_)
                    {
                        int32_t x, y;
                        expander->get_xy(n->get_id(), x, y);
                        std::cerr << "  generating "
                            << "(edgecost=" << cost_to_n<<") " 
                            << "("<<x<<", "<<y<<")...";
                        n->print(std::cerr);
                        std::cerr << std::endl;
                    }
                    #endif
                    nodes_generated_++;
                }

                // update the best solution if possible
                warthog::search_node* reverse_n = 
                    reverse_expander->generate(n->get_id());
                if(reverse_n->get_search_id() == n->get_search_id())
//                        && reverse_n->get_expanded())
                {
                    if((current->get_g() + cost_to_n + reverse_n->get_g()) < best_cost_)
                    {
                        //v_ = current;
                        v_ = n;
                        w_ = reverse_n;
                        best_cost_ = current->get_g() + cost_to_n + reverse_n->get_g();

                        #ifndef NDEBUG
                        if(verbose_)
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
            if(verbose_)
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

