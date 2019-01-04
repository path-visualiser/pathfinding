#ifndef WARTHOG_CBS_LL_HEURISTIC_H
#define WARTHOG_CBS_LL_HEURISTIC_H

// mapf/cbs_ll_heuristic.h
//
// The low-level (i.e. single-agent) heuristic function used in 
// Conflict-based Search. Pre-computes all distances, from every 
// target node to every other node in the input graph. 
//
// This implementation assumes the input graph is a 4-connected 
// uniform-cost grid. 
//
// For more details see: 
// Sharon, Guni, et al. "Conflict-based search for optimal multi-agent pathfinding." 
// Artificial Intelligence 219 (2015): 40-66.
//
// @author: dharabor
// @created: 2018-11-04
//

#include "forward.h"
#include "labelled_gridmap.h"
#include <map>

namespace warthog
{

class cbs_ll_heuristic
{
    public:
        cbs_ll_heuristic() 
        {
        }

        ~cbs_ll_heuristic() 
        {
        }

        inline warthog::cost_t
        h(warthog::sn_id_t id, warthog::sn_id_t id2)
        {
            assert(t_map_.find((uint32_t)id2) != t_map_.end());
            assert(t_map_.find((uint32_t)id2)->second == t_index_);
            return h_[t_index_][(uint32_t)id];
        }

        // run a Dijkstra grid search (not time expanded!) from each 
        // target whose identifier is  in @param target_nodes. 
        // We store the g-value from each target to every node in the grid.
        // These precomputed g-values are a lower-bound on the
        // distance from any location@time to each target.
        void
        compute_h_values(
                std::vector<uint32_t>& target_nodes,
                warthog::gridmap* gm);

        // the current target specifies which set of g-values to
        // refer to when answering ::h queries
        // this function returns true if @param target_id (unpadded)
        // corresponds to one of the nodes for which a set
        // of pre-computed g-values exist. 
        // Otherwise the function returns false.
        //
        // @param target_id: unpadded xy index specifying the current target
        // @return true indicates success, false indicates failure
        bool
        set_current_target(warthog::sn_id_t target_id)
        {
            uint32_t target_xy_id = (uint32_t)target_id;
            std::map<uint32_t, uint32_t>::iterator it = t_map_.find(target_xy_id);
            if(it == t_map_.end()) { return false; }
            t_index_ = it->second;
            return true;
        }

        size_t
        mem() 
        { 
            size_t sz = 0;
            for(uint32_t i = 0; i < h_.size(); i++)
            {
                sz += sizeof(warthog::cost_t) * h_[i].size();
            }
            sz += sizeof(void*)*h_.size();
            sz += sizeof(uint32_t)*2*t_map_.size();
            sz += sizeof(this);
            return sz;
        }

    private:
        std::vector<std::vector<warthog::cost_t>> h_;
        std::map<uint32_t, uint32_t> t_map_;
        uint32_t t_index_;
};

}

#endif
