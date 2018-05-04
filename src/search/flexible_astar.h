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
template< class H, 
          class E >
class flexible_astar : public warthog::search
{
	public:
		flexible_astar(H* heuristic, E* expander) :
            heuristic_(heuristic), expander_(expander)
		{
			open_ = new warthog::pqueue(1024, true);
            cost_cutoff_ = warthog::INF;
            exp_cutoff_ = warthog::INF;
            on_relax_fn_ = 0;
            on_generate_fn_ = 0;
            on_expand_fn_ = 0;
		}

		virtual ~flexible_astar()
		{
			delete open_;
		}

        virtual void
		get_distance(
                warthog::problem_instance& instance, warthog::solution& sol)
        {

            assert(
                sol.path_.empty() && sol.sum_of_edge_costs_ == warthog::INF);
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
                    current = expander_->generate(current->get_parent());
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
                            expander_->generate(node_id);
                        assert(n->get_search_id() == pi_.instance_id_);
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
                warthog::search_node* current = expander_->generate(i);
                if(current->get_search_id() == pi_.instance_id_) 
                { 
                    coll.push_back(current);
                }
            }
        }

        warthog::search_node* 
        get_search_node(uint32_t id)
        {
            return expander_->get_ptr(id, pi_.instance_id_);
        }

        // apply @param fn to every node on the closed list
        void
        apply_to_closed(std::function<void(warthog::search_node*)>& fn)
        {
            for(uint32_t i = 0; i < expander_->get_nodes_pool_size(); i++)
            {
                warthog::search_node* current = expander_->generate(i);
                if(current->get_search_id() == pi_.instance_id_)
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
                    double edge_cost, 
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
		warthog::pqueue* open_;
        warthog::problem_instance pi_;

        // early termination limits
        double cost_cutoff_; 
        uint32_t exp_cutoff_;

        // callback for when a node is relaxed
        std::function<void(warthog::search_node*)>* on_relax_fn_;

        // callback for when a node is reached / generated
        std::function<void(
                warthog::search_node*, 
                warthog::search_node*, 
                double edge_cost, 
                uint32_t edge_id)>* on_generate_fn_;

        // callback for when a node is expanded
        std::function<void(warthog::search_node*)>* on_expand_fn_;

		// no copy ctor
		flexible_astar(const flexible_astar& other) { } 
		flexible_astar& 
		operator=(const flexible_astar& other) { return *this; }

		warthog::search_node*
		search(warthog::solution& sol)
		{
			warthog::timer mytimer;
			mytimer.start();
			open_->clear();

            // generate the start node. 
            // also, update the instance with the start and target internal
            // ids (this is just to make debugging easier)
			warthog::search_node* start;
            if(pi_.start_id_ == warthog::INF) { return 0; }
            start = expander_->generate_start_node(&pi_);
            pi_.start_id_ = start->get_id();
            if(on_generate_fn_) { (*on_generate_fn_)(start, 0, 0, UINT32_MAX); }

			warthog::search_node* target = 0;
            if(pi_.target_id_ != warthog::INF)
            { 
                target = expander_->generate_target_node(&pi_); 
                pi_.target_id_ = target->get_id();
                target = 0; 
            }

            // initialise start node
            int32_t sx, sy, gx, gy;
            expander_->get_xy(pi_.start_id_, sx, sy);
            expander_->get_xy(pi_.target_id_, gx, gy);
			start->init(pi_.instance_id_, warthog::NODE_NONE, 
                    0, heuristic_->h(sx, sy, gx, gy));

			#ifndef NDEBUG
			if(pi_.verbose_) { pi_.print(std::cerr); std:: cerr << "\n";}
			#endif

            // begin expanding
			open_->push(start);
			while(open_->size())
			{
				warthog::search_node* current = open_->pop();
				current->set_expanded(true); // NB: set before generating
				assert(current->get_expanded());
				sol.nodes_expanded_++;
                if(on_expand_fn_) { (*on_expand_fn_)(current); }

                // goal test
                if(current->get_id() == pi_.target_id_) 
                {
                    target = current;
                    break;
                }
                    
                // early termination: in case we want bounded-cost 
                // search or if we want to impose some memory limit
                if(current->get_f() > cost_cutoff_) { break; } 
                if(sol.nodes_expanded_ >= exp_cutoff_) { break; }


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

                // generate successors
				expander_->expand(current, &pi_);
				warthog::search_node* n = 0;
				double cost_to_n = warthog::INF;
                uint32_t edge_id = 0;
				for(expander_->first(n, cost_to_n); 
						n != 0;
					   	expander_->next(n, cost_to_n))
				{
                    sol.nodes_touched_++;
                    if(on_generate_fn_) 
                    { (*on_generate_fn_)(n, current, cost_to_n, edge_id++); }
                    
                    // add new nodes to the fringe
                    if(n->get_search_id() != current->get_search_id())
                    {
						double gval = current->get_g() + cost_to_n;
                        int32_t nx, ny;
                        expander_->get_xy(n->get_id(), nx, ny);
                        n->init(current->get_search_id(), current->get_id(),
                            gval, gval + heuristic_->h(nx, ny, gx, gy));

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

                        if(on_relax_fn_) { (*on_relax_fn_)(n); }
                        continue;
                    }

                    // skip neighbours already expanded
					if(n->get_expanded())
					{
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

                    // update a node from the fringe
					if(open_->contains(n))
					{
						double gval = current->get_g() + cost_to_n;
						if(gval < n->get_g())
						{
							n->relax(gval, current->get_id());
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
                            if(on_relax_fn_) { (*on_relax_fn_)(n); }
						}
						else
						{
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
				}
			}

			mytimer.stop();
			sol.time_elapsed_micro_ = mytimer.elapsed_time_micro();

            #ifndef NDEBUG
            if(pi_.verbose_)
            {
                if(target == 0) 
                {
                    std::cerr 
                        << "search failed; no solution exists " << std::endl;
                }
                else
                {
                    int32_t x, y;
                    expander_->get_xy(target->get_id(), x, y);
                    std::cerr << "target found ("<<x<<", "<<y<<")...";
                    target->print(std::cerr);
                    std::cerr << std::endl;
                }
            }
            #endif

            return target;
		}
};

}

#endif

