#ifndef WARTHOG_CBS_SA_HEURISTIC_H
#define WARTHOG_CBS_SA_HEURISTIC_H

// heuristics/cbs_sa_heuristic.h
//
// The single-agent heuristic function used in Conflict-based Search.
// Pre-computes all distances, from every target node to every other
// node in the input graph. 
//
// This implementation assumes the input graph is a 4-connected 
// uniform-cost grid. 
//
// @author: dharabor
// @created: 2018-11-04
//

#include "gridmap.h"
#include <map>

namespace warthog
{

class cbs_sa_heuristic
{
    public:
        cbs_sa_heuristic(warthog::gridmap* gm);
        ~cbs_sa_heuristic() { } 

        inline double
        h(int32_t x, int32_t y, int32_t x2, int32_t y2)
        {
            return h(y * gm_->width() + x, y2 * gm_->width() + x2);
        }

        inline double
        h(unsigned int id, unsigned int id2)
        {
            return h_[t_index_][id & id_mask_];
        }

        // run a Dijkstra grid search (not time expanded!) from each 
        // target whose identifier is  in @param target_nodes. 
        // We store the g-value from each target to every node in the grid.
        // These precomputed g-values are a lower-bound on the
        // distance from any location@time to each target.
        void
        compute_h_values(std::vector<uint32_t>& target_nodes);

        // the current target specifies which set of g-values to
        // refer to when answering ::h queries
        // this function returns true if @param padded_id
        // corresponds to one of the nodes for which a set
        // of pre-computed g-values exist. 
        // Otherwise the function returns false.
        bool
        set_current_target(uint32_t padded_id)
        {
            std::map<uint32_t, uint32_t>::iterator it = t_map_.find(padded_id);
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
                sz += sizeof(uint32_t) * h_[i].size();
            }
            sz += sizeof(void*)*h_.size();
            sz += sizeof(uint32_t)*2*t_map_.size();
            sz += sizeof(this);
            return sz;
        }

    private:
        warthog::gridmap* gm_;
        std::vector<std::vector<double>> h_;
        std::map<uint32_t, uint32_t> t_map_;
        uint32_t t_index_;
        uint32_t id_mask_;
};

}

#endif
