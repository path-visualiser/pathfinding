#ifndef WARTHOG_MANHATTAN_TIME_EXPANSION_POLICY_H
#define WARTHOG_MANHATTAN_TIME_EXPANSION_POLICY_H

// search/manhattan_time_expansion_policy.h
//
// An time-based expansion policy for uniform-cost manhattan grids.
// At each location there are two different types of actions:
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

// this data structure describes the constraints imposed on a single 
// grid cell. a brief description of the members:
//
// timestep_ is the time at which the constraints apply
// v_ is true if the cell is constrained / blocked
// e_ indicates which of the 8 outgoing directions are constrained / blocked
struct cell_constraints
{
    cell_constraints() 
        : timestep_(0), v_(0), e_(0) { }

    cell_constraints(uint16_t timestep, uint8_t v, uint8_t e) 
        : timestep_(timestep), v_(v), e_(e) { }

    uint16_t timestep_; 
    uint8_t v_;
    uint8_t e_;
};

// this data structure describes the all constraints that currently apply
class time_constraints
{
    public:

       time_constraints(warthog::gridmap* map)
       {
           uint32_t map_sz = map->height() * map->width();
           cons_ = new std::vector< std::vector<cell_constraints> >(map_sz);
       } 

       ~time_constraints()
       {
           delete cons_;
       }

       inline void
       add_constraint(uint32_t padded_id, cell_constraints con)
       {
           assert(padded_id < cons_->size());
           cons_->at(padded_id).push_back(con);
       }

       // return all constraints associated with the xy location
       // @param padded_id
       inline std::vector<cell_constraints>& 
       get_constraints(uint32_t padded_id)
       {
           return cons_->at(padded_id);
       }

       // return any constraints associated with the xy location
       // @param padded_id at time @param timestep
       inline cell_constraints
       get_constraints(uint32_t padded_id, uint32_t timestep)
       {
            cell_constraints retval;
            std::vector<cell_constraints>::iterator con_iter = 
                std::find_if(
                        cons_->at(padded_id).begin(), 
                        cons_->at(padded_id).end(),
                    [timestep](warthog::cell_constraints& tmp)
                    -> bool
                    {
                        return tmp.timestep_  == timestep;
                    });
            if(con_iter != cons_->at(padded_id).end())
            {
                retval = *con_iter;
            }
            return retval;
       }

       void
       clear_constraints(uint32_t padded_id)
       {
           cons_->at(padded_id).clear();
       }

       void
       clear_constraints()
       {
           for(uint32_t i = 0; i < cons_->size(); i++)
           {
               cons_->at(i).clear();
           }   
       }

    private:
        std::vector< std::vector<cell_constraints> >* cons_;
};

class manhattan_time_expansion_policy 
{
	public:
		manhattan_time_expansion_policy(warthog::gridmap* map);
		~manhattan_time_expansion_policy();

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

