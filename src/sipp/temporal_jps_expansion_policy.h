#ifndef WARTHOG_TEMPORAL_JPS_EXPANSION_POLICY_H
#define WARTHOG_TEMPORAL_JPS_EXPANSION_POLICY_H

// sipp/temporal_jps_expansion_policy.h
// 
// An implementation of Temporal Jump Point Search (JPST):
// an algorithm for time-expanded gridmaps which 
// combines Jump Point Search (JPS) with Safe Interval Path Planning (SIPP)
// In the time dimension, JPST proceeds like SIPP
// In the xy dimensions, JPST proceeds like JPS
//
// @author: dharabor
// @created: 2019-10-28
//

#include "memory/node_pool.h"
#include "jps/online_jump_point_locator.h"
#include "search/problem_instance.h"
#include "sipp/jpst_locator.h"
#include "sipp/jpst_gridmap.h"
#include "util/arraylist.h"

namespace warthog
{

class temporal_jps_expansion_policy
{
    public:
        temporal_jps_expansion_policy(warthog::jpst_gridmap* gm);
        virtual ~temporal_jps_expansion_policy(); 

        // get a search_node memory pointer associated with @param node_id. 
        // (value is null if @param node_id is bigger than nodes_pool_size_)
		inline warthog::search_node*
		generate(warthog::sn_id_t node_id)
		{
            warthog::sn_id_t xy_id = node_id & UINT32_MAX;
            warthog::sn_id_t index = node_id >> 32;

            while(pool_.size() <= index)
            { pool_.push_back(new warthog::mem::node_pool(sz_xy)); }
            warthog::search_node* node = pool_.at(index)->generate(xy_id);
            node->set_id(node_id);
            return node;
		}

        // the start node is the first safe interval for the 
        // associated grid location. We require this interval
        // to start at time 0 
        inline warthog::search_node* 
        generate_start_node(warthog::problem_instance* pi)
        {
            warthog::sipp::safe_interval start_si = 
                jpst_gm_->get_safe_interval(
                    (uint32_t)pi->start_id_, 0);
            if(start_si.s_time_ != 0) { return 0; }
            return generate(pi->start_id_);
        }

        // the target node in SIPP is an xy location which 
        // the agent must reach. This function can be adjusted.
        // e.g. for some types of problems the target may be a 
        // at a specific time index or part of a specific 
        // safe interval. In MAPF the agent needs to wait 
        // indefinitely once it reaches its target 
        // location. This can only be satisfied if the safe
        // interval ends at time infinity.
        // 
        // NB: See also the function ::is_target where we need to
        // detect whether a node being expanded is the target
        warthog::search_node*
        generate_target_node(warthog::problem_instance* pi)
        {
            uint32_t xy_id = (pi->target_id_ & UINT32_MAX);
            uint32_t gm_id = jpst_gm_->gm_->to_padded_id(xy_id);
            if(!jpst_gm_->gm_->get_label(gm_id))
            {
                return 0; // target is an obstacle
            }
            target_gm_id_ = gm_id;
            return generate(pi->target_id_);
        }

        // returns true if @param n is the target node and false otherwise.
        // Currently the function assumes the objective is to reach a 
        // specific xy location as early as possible. 
        // 
        // NB: This function may need to be adjusted depending on the
        // problem at hand. In some settings the target can be time indexed.
        // In others, such as MAPF, the target must satisfy certain constraints.
        // e.g. find the earliest time at which the target's xy location can
        // be reached such that the agent can wait at that location
        // indefinitely (i.e. the safe interval must end at time INF)
        bool
        is_target(warthog::search_node* n, warthog::problem_instance* pi)
        {
            return (n->get_id() & UINT32_MAX) == (pi->target_id_ & UINT32_MAX);
        }


        // return the xy coordinates corresponding to the search node 
        // identifier @param node_id
        void
        get_xy(sn_id_t node_id, int32_t& x, int32_t& y) 
        {
            uint32_t xy_id = (uint32_t)(node_id & INT32_MAX);
            y = (int32_t)(xy_id / jpst_gm_->gm_->header_width());
            x = (int32_t)(xy_id % jpst_gm_->gm_->header_width());
        }

        // SIPP looks for successors among the set of safe intervals stored
        // with each xy location adjacent to that of @param current
        // Generated are all safe intervals which begin before the end (<=)
        // of the safe interval associated with @param current
		inline void 
		expand(warthog::search_node* c_node, warthog::problem_instance* problem)
        {
            reset();

            c_node_ = c_node;
            c_xy_id_ = (uint32_t)(c_node->get_id() & UINT32_MAX);
            c_index_ = (uint32_t)(c_node->get_id() >> 32);
            c_gm_id_ = jpst_gm_->gm_->to_padded_id(c_xy_id_);
            c_si_ = &jpst_gm_->get_safe_interval(c_xy_id_, c_index_);
            problem_ = problem;
            uint32_t p_xy_id = (uint32_t)(c_node_->get_parent() & UINT32_MAX);

            // get the parent direction
            warthog::jps::direction lastmove;
            if(p_xy_id == warthog::INF32)
            {
                lastmove = warthog::jps::NONE;
            }
            else
            {
                lastmove = warthog::jps::compute_direction_4c(c_xy_id_, 
                               p_xy_id, map_width_);
            }

            // get the tiles around the current node c and determine
            // which of the available moves are forced and which are natural
            uint32_t c_tiles;
            jpst_gm_->gm_->get_neighbours(c_gm_id_, (uint8_t*)&c_tiles);
            uint32_t succ_dirs = warthog::jps::compute_successors_4c(lastmove, c_tiles);

             // the parent can be forced if there exist temporal obstacles 
             // at that location (reasoning about time is weird, man)
             if(lastmove && jpst_gm_->get_all_intervals(p_xy_id).size() > 1)
             {
                succ_dirs |= opposite_dir[__builtin_ffs(lastmove)];
             }

            // generate successors
            for(uint32_t i = 0; i < 4; i++)
            {
                warthog::jps::direction d = (warthog::jps::direction) (1 << i);
                if(succ_dirs & d)
                {
                    double jumpcost;
                    uint32_t succ_id;
                    jpl_t_->jump(d, c_gm_id_, target_gm_id_, succ_id, jumpcost);

                    if(succ_id != warthog::INF32)
                    {
                        int32_t succ_xy_id = 
                            (int32_t)c_xy_id_ + (int32_t)(xy_id_offsets_[(i+1)]*jumpcost);
                        assert(succ_xy_id >= 0);
                        warthog::cbs::move ec_move = ec_moves_[(i+1)];
                        generate_successors((uint32_t)succ_xy_id, ec_move, jumpcost);
                    }
                }
            }

            //if(p_index < (jpst_gm_->get_all_intervals(p_xy_id).size()-1))
            //{ jump_east( c_xy_id_, c_gm_id_ ); }
        }

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

        size_t
        mem()
        {
            size_t retval = 0;
            for(uint32_t i = 0; i < pool_.size(); i++)
            {
                retval += pool_.at(i)->mem();
            }
            retval += sizeof(this);
            retval += sizeof(neighbour_record)*neis_->size();
            return retval;
        }

    private:
        struct neighbour_record
        {
            neighbour_record(warthog::search_node* node, warthog::cost_t cost)
            {
                node_ = node;
                cost_ = cost;
            }
            warthog::search_node* node_;
            warthog::cost_t cost_;
        };


        uint32_t sz_xy;
        uint32_t current_;
        arraylist<neighbour_record>* neis_;

        warthog::jpst_gridmap* jpst_gm_;
        warthog::jpst_locator* jpl_t_;
        int32_t xy_id_offsets_[5];
        warthog::cbs::move ec_moves_[5];
        warthog::jps::direction opposite_dir[5];

        // convenience variables
        warthog::search_node* c_node_;
        uint32_t target_gm_id_;
        uint32_t c_xy_id_;
        uint32_t c_index_;
        uint32_t c_gm_id_;
        warthog::sipp::safe_interval* c_si_;

        // convenience variables
        uint32_t map_width_;
        uint32_t gm_map_width_;
        warthog::problem_instance* problem_;

        // we need to hash each sipp::safe_interval to a warthog::search_node 
        // and each warthog::search_node to a sipp::safe_interval
        // we use several pools of preallocated memory for this:
        // one pool for the set of sipp::safe_intervals which have index 0
        // one pool for the set of sipp::safe_intervals which have index 1
        // one pool for the set of sipp::safe_intervals which have index 2
        // ...
        // and so on
        // the upper 4 bytes of each warthog::search_node::id specifies
        // the (x, y) location of the node containing the sipp::safe_interval
        // the lower 4 bytes specifies the index of the sipp::safe_interval
        std::vector<warthog::mem::node_pool*> pool_;

        // generates all safe-interval successors at location @param succ_xy_id
        // the location may be adjacent to the node being currently expanded or
        // it can be a location reached through jumping. If the location is
        // reached through jumping, the code expects that all ``jumped over''
        // locations have a single safe interval: [0, INF)
        inline void 
        generate_successors( uint32_t succ_xy_id, warthog::cbs::move ec_direction, 
                             warthog::cost_t action_cost )
        {
            // iterate over adjacent safe intervals
            std::vector<warthog::sipp::safe_interval>& neis 
                = jpst_gm_->get_all_intervals(succ_xy_id);
            for(uint32_t i = 0; i < neis.size(); i++)
            {
                // we generate safe intervals for adjacent cells but:
                // (i) only if the successor safe interval begins before 
                // (i.e. <) the end of the current safe interval and; 
                // (ii) only if the current safe interval is safe for
                // the duration of the action that moves the agent
                // (iii) only the successor is safe at the time the 
                // agent finishes moving.
                warthog::sipp::safe_interval* succ_si = &(neis.at(i));

                if( succ_si->s_time_ < c_si_->e_time_ && 
                    c_node_->get_g() < succ_si->e_time_)
                {
                    // if the adjacent safe interval begins at some time
                    // in the future then we wait at the current
                    // safe interval and move away at the earliest time
                    if(succ_si->s_time_ >= (c_node_->get_g() + action_cost))
                    {
                        // the wait+move action completes exactly as 
                        // the successor interval begins
                        action_cost = succ_si->s_time_ - c_node_->get_g();

                        // prune: avoid edge colllisions when moving
                        // into the successor interval
                        if(succ_si->action_ == ec_direction)
                        { continue; }
                    }

                    // prune: not enough time to execute the action
                    //if((c_node->get_g() + action_cost) > c_si.e_time_) 
                    //{ continue;  }

                    // prune: moved to successor but the location is not the
                    // target and not enough time remains to move again
                    //if((current->get_g() + action_cost) == succ_si.e_time_
                    //   && !(succ_xy_id == pi->target_id_))
                    //{ continue; }

                    // generate successor
                    warthog::sn_id_t succ_node_id = succ_xy_id;
                    succ_node_id = ((warthog::sn_id_t)i << 32) + succ_node_id;
                    add_neighbour(generate(succ_node_id), action_cost);
                }
            }
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
