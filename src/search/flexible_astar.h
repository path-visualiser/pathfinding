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
#include "solution.h"
#include "timer.h"

#include <functional>
#include <iostream>
#include <memory>
#include <vector>

namespace warthog
{

// H is a heuristic function
// E is an expansion policy
// F is a node filtering (== pruning) policy
template< class H, 
          class E, 
          class F = warthog::dummy_filter>
class flexible_astar : public warthog::search
{
	public:
		flexible_astar(H* heuristic, E* expander, F* filter = 0) :
            heuristic_(heuristic), expander_(expander), filter_(filter)
		{
			open_ = new warthog::pqueue(1024, true);
            cost_cutoff_ = warthog::INF;
            exp_cutoff_ = warthog::INF;
            on_relax_fn_ = [](warthog::search_node*){ };
            on_norelax_fn_ = 
                [](warthog::search_node*, warthog::search_node*, double){ };
		}

		virtual ~flexible_astar()
		{
			cleanup();
			delete open_;
		}

        void
		get_path(warthog::problem_instance& instance, warthog::solution& sol)
		{
            assert(
                sol.path_.empty() && sol.sum_of_edge_costs_ == warthog::INF);
            pi_ = instance;

			warthog::search_node* target = search(sol);
			if(target)
			{
                sol.sum_of_edge_costs_ = target->get_g();

				// follow backpointers to extract the path
				assert(target->get_id() == pi_.target_id_);
                warthog::search_node* current = target;
				while(current)
                {
					sol.path_.push_back(current->get_id());
                    current = current->get_parent();
				}
				assert(sol.path_.back() == pi_.start_id_);

                #ifndef NDEBUG
                if(pi_.verbose_)
                {
                    for(auto& node_id : sol.path_)
                    {
                        int32_t x, y;
                        expander_->get_xy(node_id, x, y);
                        std::cerr 
                            << "final path: (" << x << ", " << y << ")...";
                        warthog::search_node* n = 
                            expander_->get_ptr(node_id, pi_.instance_id_);
                        assert(n);
                        n->print(std::cerr);
                        std::cerr << std::endl;
                    }
                }
                #endif
            }
		}
        
        // return a list of the nodes expanded during the last search
        // @param coll: an empty list
        void
        closed_list(std::vector<warthog::search_node*>& coll)
        {
            for(uint32_t i = 0; i < expander_->get_nodes_pool_size(); i++)
            {
                warthog::search_node* current = 
                    expander_->get_ptr(i, pi_.instance_id_);
                if(current) { coll.push_back(current); }
            }
        }

        // apply @param fn to every node on the closed list
        void
        apply_to_closed(std::function<void(warthog::search_node*)>& fn)
        {
            for(uint32_t i = 0; i < expander_->get_nodes_pool_size(); i++)
            {
                warthog::search_node* current = 
                    expander_->get_ptr(i, pi_.instance_id_);
                if(current) { fn(current); }
            }
        }

        // apply @param fn every time a node is successfully relaxed
        void
        apply_on_relax(std::function<void(warthog::search_node*)>& fn)
        {
            on_relax_fn_ = fn;
        }

        // apply @param fn every time a node is reached from a new parent
        // with the same cost as a previous parent
        void
        apply_on_norelax( std::function<void(
                    warthog::search_node* n, 
                    warthog::search_node* current, 
                    double edge_cost)>& fn)
        {
            on_norelax_fn_ = fn;
        }

        // set a cost-cutoff to run a bounded-cost A* search.
        // the search terminates when the target is found or the f-cost 
        // limit is reached.
        inline void
        set_cost_cutoff(double cutoff) { cost_cutoff_ = cutoff; }

        inline double
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
        F* filter_;
		warthog::pqueue* open_;
        warthog::problem_instance pi_;

        // early termination limits
        double cost_cutoff_; 
        uint32_t exp_cutoff_;

        // callback for when a node is relaxex
        std::function<void(warthog::search_node*)> on_relax_fn_;

        // callback for when a node is not relaxed
        std::function<void(
                warthog::search_node*, 
                warthog::search_node*, 
                double edge_cost)> on_norelax_fn_;

		// no copy ctor
		flexible_astar(const flexible_astar& other) { } 
		flexible_astar& 
		operator=(const flexible_astar& other) { return *this; }

		warthog::search_node*
		search(warthog::solution& sol)
		{
			warthog::timer mytimer;
			mytimer.start();

            cleanup();

            // generate the start and goal. then
            // update the instance with their internal ids 
            // (this is just to make debugging easier)
			warthog::search_node* start;
            if(pi_.start_id_ == warthog::INF) { return 0; }
            start = expander_->generate_start_node(&pi_);
            pi_.start_id_ = start->get_id();

			warthog::search_node* target = 0;
            if(pi_.target_id_ != warthog::INF)
            { 
                target = expander_->generate_target_node(&pi_); 
                pi_.target_id_ = target->get_id();
                target = 0; // just need the id; FIXME: hacky 
            }

			#ifndef NDEBUG
			if(pi_.verbose_) { pi_.print(std::cerr); std:: cerr << "\n";}
			#endif

            int32_t sx, sy, gx, gy;
            expander_->get_xy(pi_.start_id_, sx, sy);
            expander_->get_xy(pi_.target_id_, gx, gy);

			start->init(pi_.instance_id_, 0, 0, 
                    heuristic_->h(sx, sy, gx, gy));
			open_->push(start);

			while(open_->size())
			{
				if(open_->peek()->get_id() == pi_.target_id_)
				{
					#ifndef NDEBUG
					if(pi_.verbose_)
					{
						int32_t x, y;
						warthog::search_node* current = open_->peek();
                        expander_->get_xy(current->get_id(), x, y);
						std::cerr << "target found ("<<x<<", "<<y<<")...";
						current->print(std::cerr);
						std::cerr << std::endl;
					}
					#endif
					target = open_->peek();
					break;
				}

                // early termination tests (in case we want bounded-cost 
                // search or if we want to impose some memory limit)
                if(open_->peek()->get_f() > cost_cutoff_) { break; } 
                if(sol.nodes_expanded_ >= exp_cutoff_) { break; }

				warthog::search_node* current = open_->pop();
				sol.nodes_expanded_++;

				#ifndef NDEBUG
				if(pi_.verbose_)
				{
					int32_t x, y;
                    expander_->get_xy(current->get_id(), x, y);
					std::cerr 
                        << sol.nodes_expanded_
                        << ". expanding ("<<x<<", "<<y<<")...";
					current->print(std::cerr);
					std::cerr << std::endl;
				}
				#endif
				current->set_expanded(true); // NB: set before generating
				assert(current->get_expanded());
				expander_->expand(current, &pi_);

				warthog::search_node* n = 0;
				double cost_to_n = warthog::INF;
				for(expander_->first(n, cost_to_n); 
						n != 0;
					   	expander_->next(n, cost_to_n))
				{
                    sol.nodes_touched_++;
					if(n->get_expanded())
					{
						// skip neighbours already expanded
                        #ifndef NDEBUG
                        if(pi_.verbose_)
                        {
                            int32_t x, y;
                            expander_->get_xy(n->get_id(), x, y);
                            std::cerr << "  closed; (edgecost=" 
                                << cost_to_n << ") ("<<x<<", "<<y<<")...";
                            n->print(std::cerr);
                            std::cerr << std::endl;
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
                            sol.nodes_updated_++;

							#ifndef NDEBUG
							if(pi_.verbose_)
							{
								int32_t x, y;
                                expander_->get_xy(n->get_id(), x, y);
								std::cerr 
                                    << "  open; updating (edgecost="
                                    << cost_to_n<<") ("<<x<<", "<<y<<")...";
								n->print(std::cerr);
								std::cerr << std::endl;
							}
							#endif
                            on_relax_fn_(n);
						}
						else
						{
                            on_norelax_fn_(n, current, cost_to_n);

							#ifndef NDEBUG
							if(pi_.verbose_)
							{
								int32_t x, y;
                                expander_->get_xy(n->get_id(), x, y);
								std::cerr 
                                    << "  open; not updating (edgecost=" 
                                    << cost_to_n<< ") ("<<x<<", "<<y<<")...";
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
                        int32_t nx, ny;
                        expander_->get_xy(n->get_id(), nx, ny);
                        n->init(pi_.instance_id_, current, 
                            gval, gval + heuristic_->h(nx, ny, gx, gy));
                        
                        // but only if the node is not provably redundant
                        if(filter_ && filter_->filter(n))
                        {
                            #ifndef NDEBUG
                            if(pi_.verbose_)
                            {
                                std::cerr 
                                    << "  filtered-out (edgecost=" 
                                    << cost_to_n<<") ("<<nx<<", "<<ny<<")...";
                                n->print(std::cerr);
                                std::cerr << std::endl;
                            }
                            #endif
                            continue;
                        }

                        open_->push(n);
                        sol.nodes_inserted_++;

                        #ifndef NDEBUG
                        if(pi_.verbose_)
                        {
                            std::cerr 
                                << "  generating (edgecost=" 
                                << cost_to_n<<") ("<< nx <<", "<< ny <<")...";
                            n->print(std::cerr);
                            std::cerr << std::endl;
                        }
                        #endif

                        on_relax_fn_(n);
					}
				}
			}

            #ifndef NDEBUG
            if(pi_.verbose_)
            {
                if(target == 0) 
                {
                    std::cerr 
                        << "search failed; no solution exists " << std::endl;
                }
            }
            #endif

			mytimer.stop();

			sol.sum_of_edge_costs_ = mytimer.elapsed_time_micro();
            return target;
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

