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
template <class H, class E>
class flexible_astar : public warthog::search
{
	public:
		flexible_astar(H* heuristic, E* expander)
			: heuristic_(heuristic), expander_(expander)
		{
			open_ = new warthog::pqueue(1024, true);
			verbose_ = false;
            hscale_ = 1.0;
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
			cleanup();
			return path;
		}

		double
		get_length(uint32_t startid, uint32_t goalid)
		{
			warthog::search_node* goal = search(startid, goalid);
			warthog::cost_t len = warthog::INF;
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
			cleanup();
			return len / (double)warthog::ONE;
		}

		virtual inline size_t
		mem()
		{
			size_t bytes = 
				// memory for the priority quete
				open_->mem() + 
				// gridmap size and other stuff needed to expand nodes
				expander_->mem() +
				// misc
				sizeof(*this);
			return bytes;
		}

        inline double
        get_hscale() { return hscale_; } 

        inline void
        set_hscale(double hscale) { hscale_ = hscale; } 



	private:
		H* heuristic_;
		E* expander_;
		warthog::pqueue* open_;
        double hscale_; // heuristic scaling factor

		// no copy
		flexible_astar(const flexible_astar& other) { } 
		flexible_astar& 
		operator=(const flexible_astar& other) { return *this; }

		warthog::search_node*
		search(uint32_t startid, uint32_t goalid)
		{
			nodes_expanded_ = nodes_generated_ = nodes_touched_ = 0;
			search_time_ = 0;

			warthog::timer mytimer;
			mytimer.start();

			#ifndef NDEBUG
			if(verbose_)
			{
				std::cerr << "search: startid="<<startid<<" goalid=" <<goalid
					<< std::endl;
			}
			#endif

			warthog::problem_instance instance;
			instance.set_goal(goalid);
			instance.set_start(startid);
			instance.set_searchid(++warthog::search::searchid_);

			warthog::search_node* goal = 0;
			warthog::search_node* start = expander_->generate(startid);
			start->reset(instance.get_searchid());
			start->set_g(0);
			start->set_f(heuristic_->h(startid, goalid) * hscale_);
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
				nodes_expanded_++;

				warthog::search_node* current = open_->pop();
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
				warthog::cost_t cost_to_n = warthog::INF;
				for(expander_->first(n, cost_to_n); 
						n != 0;
					   	expander_->next(n, cost_to_n))
				{
					nodes_touched_++;
					if(n->get_expanded())
					{
						// skip neighbours already expanded
                        if(verbose_)
                        {
                            int32_t x, y;
                            expander_->get_xy(n, x, y);
                            std::cerr << "  closed; (edgecost=" << cost_to_n << ") ("<<x<<", "<<y<<")...";
                            n->print(std::cerr);
                            std::cerr << std::endl;

                            //warthog::cost_t gval = current->get_g() + cost_to_n;
                            //double hval = heuristic_->h(n->get_id(), goalid) * hscale_;
                            //warthog::cost_t fval = gval + hval;
                            //std::cerr << "  alt g: " << gval << " h: " << hval << " f: " << fval << std::endl;
                            //if(gval < n->get_g())
                            //{
                            //    assert(gval >= n->get_g());
                            //}

                        }
						continue;
					}

					if(open_->contains(n))
					{
						// update a node from the fringe
						warthog::cost_t gval = current->get_g() + cost_to_n;
						if(gval < n->get_g())
						{
							n->relax(gval, current);
							open_->decrease_key(n);
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
						warthog::cost_t gval = current->get_g() + cost_to_n;
						n->set_g(gval);
						n->set_f(gval + heuristic_->h(n->get_id(), goalid) * hscale_);
					   	n->set_parent(current);
						open_->push(n);
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
						nodes_generated_++;
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

