#ifndef FLEXIBLE_ASTAR_H
#define FLEXIBLE_ASTAR_H

// flexible_astar.h
//
// A* implementation that allows arbitrary combinations of 
// (weighted) heuristic functions and node expansion policies.
// This implementation uses a binary heap for the open_ list
// and a bit array for the closed_ list.
//
// TODO: is it better to store a separate closed list and ungenerate nodes
// or use more memory and not ungenerate until the end of search??
// 32bytes vs... whatever unordered_map overhead is a two integer key/value pair
// 
// @author: dharabor
// @created: 21/08/2012
//

#include "cpool.h"
#include "dummy_filter.h"
#include "pqueue.h"
#include "problem_instance.h"
#include "search.h"
#include "search_node.h"
#include "timer.h"

#include <iostream>
#include <memory>
#include <stack>

namespace warthog
{

// H is a heuristic function
// E is an expansion policy
// F is a node filtering (== pruning) policy
template <class H, class E, class F = warthog::dummy_filter>
class flexible_astar : public warthog::search
{
	public:
		flexible_astar(H* heuristic, E* expander, 
                F* filter = new warthog::dummy_filter())
			: heuristic_(heuristic), expander_(expander), nf_(filter)
		{
			open_ = new warthog::pqueue(1024, true);
			verbose_ = false;
            cost_cutoff_ = warthog::INF;
            exp_cutoff_ = warthog::INF;
		}

		virtual ~flexible_astar()
		{
			cleanup();
			delete open_;
		}

		inline std::stack<uint32_t>
		get_path(uint32_t startid, uint32_t goalid)
		{
			std::stack<uint32_t> path;
			warthog::search_node* goal = search(startid, goalid);
			if(goal)
			{
				// follow backpointers to extract the path
				assert(goal->get_id() == goalid);
				for(warthog::search_node* cur = goal;
						cur != 0;
					    cur = cur->get_parent())
                {
					path.push(cur->get_id());
				}
				assert(path.top() == startid);
			}
			//cleanup();
			return path;
		}
        
        // return a list of the nodes expanded during the last search
        // @param coll: an empty list
        void
        closed_list(std::vector<warthog::search_node*>& coll)
        {
            for(uint32_t i = 0; i < expander_->get_num_nodes(); i++)
            {
                warthog::search_node* current = expander_->get_ptr(i, searchid_);
                if(current) { coll.push_back(current); }
            }
        }

        // TODO: need a better and more general way to interact with the closed list
        double
        max_g_on_closed()
        {
            double max = 0;
            for(uint32_t i = 0; i < expander_->get_num_nodes(); i++)
            {
                warthog::search_node* current = expander_->get_ptr(i, searchid_);
                if(current && current->get_g() > max)
                {
                    max = current->get_g();
                }
            }
            return max;
        }

        // no cleanup after search
		double
		get_length(uint32_t startid, uint32_t goalid)
		{
			warthog::search_node* goal = search(startid, goalid);
			double len = warthog::INF;
			if(goal)
			{
				assert(goal->get_id() == goalid);
				len = goal->get_g();
			}

#ifndef NDEBUG

			if(verbose_)
			{
				std::stack<warthog::search_node*> path;
				warthog::search_node* current = goal;
				while(current != 0)	
				{
					path.push(current);
					current = current->get_parent();
				}

				while(!path.empty())
				{
					warthog::search_node* n = path.top();
                    int32_t x, y;
                    expander_->get_xy(n, x, y);
					std::cerr 
                        << "final path: (" << x << ", " << y << ")...";
					n->print(std::cerr);
					std::cerr << std::endl;
					path.pop();
				}
			}
#endif
            //cleanup();
			return len;
		}

        // set a cost-cutoff to run a bounded-cost A* search.
        // the search terminates when the goal is found or the f-cost 
        // limit is reached.
        inline void
        set_cost_cutoff(double cutoff) { cost_cutoff_ = cutoff; }

        inline double
        get_cost_cutoff() { return cost_cutoff_; }

        // set a cutoff on the maximum number of node expansions.
        // the search terminates when the goal is found or when
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
        F* nf_;
		warthog::pqueue* open_;
        double cost_cutoff_; // early termination
        uint32_t exp_cutoff_;
        uint32_t last_searchid_;

		// no copy
		flexible_astar(const flexible_astar& other) { } 
		flexible_astar& 
		operator=(const flexible_astar& other) { return *this; }

		warthog::search_node*
		search(uint32_t startid, uint32_t goalid)
		{
            cleanup();
			nodes_expanded_ = nodes_generated_ = 0;
            nodes_touched_ = heap_ops_ = 0;
			search_time_ = 0;

			warthog::timer mytimer;
			mytimer.start();

			warthog::problem_instance instance;
			instance.set_goal(goalid);
			instance.set_start(startid);
			instance.set_searchid(++warthog::search::searchid_);

			#ifndef NDEBUG
			if(verbose_)
			{
				std::cerr << "search: startid="<<startid<<" goalid=" <<goalid 
                    << " (searchid: " << instance.get_searchid() 
                    << ")" << std::endl;
			}
			#endif

			warthog::search_node* goal = 0;
			warthog::search_node* start = expander_->generate(startid);
			start->init(instance.get_searchid(), 0, 0, 
                    heuristic_->h(startid, goalid));
			open_->push(start);

			while(open_->size())
			{
				nodes_touched_++;
				if(open_->peek()->get_id() == goalid)
				{
					#ifndef NDEBUG
					if(verbose_)
					{
						int32_t x, y;
						warthog::search_node* current = open_->peek();
                        expander_->get_xy(current, x, y);
						std::cerr << "goal found ("<<x<<", "<<y<<")...";
						current->print(std::cerr);
						std::cerr << std::endl;
					}
					#endif
					goal = open_->peek();
					break;
				}
                // bounded-cost search runs until the goal is found or until
                // a maximum cost limit is reached
                if(open_->peek()->get_f() > cost_cutoff_) { break; } 

                if(nodes_expanded_ >= exp_cutoff_) { break; }

				warthog::search_node* current = open_->pop();
                heap_ops_++;
				nodes_expanded_++;

				#ifndef NDEBUG
				if(verbose_)
				{
					int32_t x, y;
                    expander_->get_xy(current, x, y);
					std::cerr << this->nodes_expanded_ << ". expanding ("<<x<<", "<<y<<")...";
					current->print(std::cerr);
					std::cerr << std::endl;
				}
				#endif
				current->set_expanded(true); // NB: set this before calling expander_ 
				assert(current->get_expanded());
				expander_->expand(current, &instance);

				warthog::search_node* n = 0;
				double cost_to_n = warthog::INF;
				for(expander_->first(n, cost_to_n); 
						n != 0;
					   	expander_->next(n, cost_to_n))
				{
					nodes_touched_++;
					if(n->get_expanded())
					{
						// skip neighbours already expanded
                        #ifndef NDEBUG
                        if(verbose_)
                        {
                            int32_t x, y;
                            expander_->get_xy(n, x, y);
                            std::cerr << "  closed; (edgecost=" << cost_to_n << ") ("<<x<<", "<<y<<")...";
                            n->print(std::cerr);
                            std::cerr << std::endl;

                            //double gval = current->get_g() + cost_to_n;
                            //double hval = heuristic_->h(n->get_id(), goalid) * hscale_;
                            //double fval = gval + hval;
                            //std::cerr << "  alt g: " << gval << " h: " << hval << " f: " << fval << std::endl;
                            //if(gval < n->get_g())
                            //{
                            //    assert(gval >= n->get_g());
                            //}

                        }
                        #endif
						continue;
					}

					if(open_->contains(n))
					{
						// update a node from the fringe
						double gval = current->get_g() + cost_to_n;
						if(gval < n->get_g())
						{
							n->relax(gval, current);
							open_->decrease_key(n);
                            heap_ops_++;
							#ifndef NDEBUG
							if(verbose_)
							{
								int32_t x, y;
                                expander_->get_xy(n, x, y);
								std::cerr << "  open; updating (edgecost="<< cost_to_n<<") ("<<x<<", "<<y<<")...";
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
                                expander_->get_xy(n, x, y);
								std::cerr << "  open; not updating (edgecost=" << cost_to_n<< ") ("<<x<<", "<<y<<")...";
								n->print(std::cerr);
								std::cerr << std::endl;
							}
							#endif
						}
					}
					else
					{
						// add a new node to the fringe
						double gval = current->get_g() + cost_to_n;
                        n->init(instance.get_searchid(), current, 
                            gval, gval + heuristic_->h(n->get_id(), goalid));
                        
                        // but only if the node is not provably redundant
                        if(nf_->filter(n))
                        {
#ifndef NDEBUG
                            if(verbose_)
                            {
                                int32_t x, y;
                                expander_->get_xy(n, x, y);
                                std::cerr << "  filtered-out (edgecost=" << cost_to_n<<") ("<<x<<", "<<y<<")...";
                                n->print(std::cerr);
                                std::cerr << std::endl;
                            }
#endif
                            continue;
                        }

                        open_->push(n);
                        nodes_generated_++;
                        heap_ops_++;

                        #ifndef NDEBUG
                        if(verbose_)
                        {
                            int32_t x, y;
                            expander_->get_xy(n, x, y);
                            std::cerr << "  generating (edgecost=" << cost_to_n<<") ("<<x<<", "<<y<<")...";
                            n->print(std::cerr);
                            std::cerr << std::endl;
                        }
                        #endif
					}
				}
//				#ifndef NDEBUG
//				if(verbose_)
//				{
//					int32_t x, y;
//                    expander_->get_xy(current, x, y);
//					std::cerr <<"closing ("<<x<<", "<<y<<")...";
//					current->print(std::cerr);
//					std::cerr << std::endl;
//			    }
//			    #endif
			}

#ifndef NDEBUG
            if(verbose_)
            {
                if(goal == 0) 
                {
                    std::cerr << "search failed; no solution exists " << std::endl;
                }
            }
#endif


			mytimer.stop();
			search_time_ = mytimer.elapsed_time_micro();
			return goal;
		}

		void
		cleanup()
		{
			open_->clear();
			expander_->clear();
		}
};

}

#endif

