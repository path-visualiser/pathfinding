#ifndef WARTHOG_GRIDMAP_TIME_EXPANSION_POLICY_H
#define WARTHOG_GRIDMAP_TIME_EXPANSION_POLICY_H

// search/gridmap_time_expansion_policy.h
//
// An time-based expansion policy for square uniform-cost grids.
// At each location there are three different types of actions:
// - diagonal moves (at a cost of sqrt(2))
// - cardinal moves (at a cost of 1)
// - wait moves (at a cost of 1)
//
// Each action, except wait, moves the agent from its current grid
// location to an adjacent grid location. Each action (including wait)
// advances time by one time-step.
//
// @author: dharabor
// @created: 2018-11-01
//

#include "expansion_policy.h"
#include "gridmap.h"
#include "search_node.h"

#include <memory>

namespace warthog
{

class problem_instance;
class gridmap_time_expansion_policy 
{
	public:
		gridmap_time_expansion_policy(warthog::gridmap* map, bool manhattan = false);
		~gridmap_time_expansion_policy();

		inline void
		reset()
		{
			current_ = 0;
            neis_->clear();
		}

		inline void
		first(warthog::search_node*& ret, double& cost)
		{
            current_ = 0;
            n(ret, cost);
		}

		inline void
		n(warthog::search_node*& ret, double& cost)
		{
            if(current_ < neis_->size())
            {
                ret = (*neis_)[current_].node_;
                cost = (*neis_)[current_].cost_;
            }
            else
            {
                ret = 0;
                cost = 0;
            }
		}

		inline void
		next(warthog::search_node*& ret, double& cost)
		{
            current_++;
            n(ret, cost);
		}

		void 
		expand(warthog::search_node*, warthog::problem_instance*);

        void
        get_xy(uint32_t node_id, int32_t& x, int32_t& y);

        warthog::search_node* 
        generate_start_node(warthog::problem_instance* pi);

        warthog::search_node*
        generate_target_node(warthog::problem_instance* pi);
        
        warthog::search_node*
        generate(uint32_t node_id)
        {
            uint32_t xy_id = node_id & id_mask_;
            uint32_t timestep = node_id >> bitwidth_map_;
            return __generate(xy_id, timestep);
        }

        inline bool
        is_target(warthog::search_node* n, warthog::problem_instance* pi)
        {
            return ((n->get_id() & id_mask_) == pi->target_id_);
        }

		size_t 
        mem();

	
	private:
		warthog::gridmap* map_;
        uint32_t id_mask_;
        uint32_t bitwidth_map_;
        bool manhattan_;
        std::vector<warthog::node_pool*>* time_map_;

        struct neighbour_record
        {
            neighbour_record(warthog::search_node* node, double cost)
            {
                node_ = node;
                cost_ = cost;
            }
            warthog::search_node* node_;
            double cost_;
        };

        arraylist<neighbour_record>* neis_;
        uint32_t current_;

        inline warthog::search_node* 
        __generate(uint32_t xy_id, uint32_t timestep)
        {
            while(timestep >= time_map_->size())
            {
                time_map_->push_back(
                            new warthog::node_pool(
                                map_->height() * map_->width()));
            }
            warthog::search_node* nei = time_map_->at(timestep)->generate(xy_id);
            nei->set_id((timestep << bitwidth_map_) | xy_id);
            return nei;
        }


        inline void 
        add_neighbour(warthog::search_node* nei, double cost)
        {
            neis_->push_back(neighbour_record(nei, cost));
            //std::cout << " neis_.size() == " << neis_->size() << std::endl;
        }

};

}

#endif

