#ifndef WARTHOG_AFHD_FILTER_H
#define WARTHOG_AFHD_FILTER_H

// afhd_filter.h
// 
// Hierarchy-based arc flags using dijkstra 
// to pre-compute arcflags on down edges
// and local neighbour analysis to label 
// up edges.
//
// @author: dharabor
// @created: 2016-09-12
//

#include <cassert>
#include <cstdint>
#include <iostream>
#include <set>
#include <vector>

#include "planar_graph.h"

namespace warthog
{

class problem_instance;
class search_node;

class afhd_filter
{
    public:
        // @param g: the search graph
        // @param part: a partitioning of the graph nodes
        afhd_filter(
                warthog::graph::planar_graph* g, 
                std::vector<uint32_t>* part);

        afhd_filter(
                warthog::graph::planar_graph* g, 
                std::vector<uint32_t>* part,
                const char* afhd_file);

        ~afhd_filter();

        // return true if the ith edge of @param node_id
        // (as specified by @param edge_index) cannot possibly 
        // appear on any optimal path to the current target;
        // return false otherwise.
        inline bool 
        filter(uint32_t node_id, uint32_t edge_idx)
        {
            uint8_t* label = flags_.at(node_id).at(edge_idx);
            return (label[t_byte_] & t_bitmask_) == 0;
        }

        void
        set_goal(uint32_t goal_id) 
        { 
            uint32_t t_part = part_->at(goal_id);
            t_byte_ = t_part >> 3;
            t_bitmask_ = 1 << (t_part & 7);
        }

        void
        compute(std::vector<uint32_t>* rank);

        void
        compute(uint32_t firstid, uint32_t lastid, 
                std::vector<uint32_t>* labels);

        void
        print(std::ostream& out);

        bool
        load_labels(const char* filename);

    private:    
        uint32_t nparts_;
        warthog::graph::planar_graph* g_;
        std::vector<uint32_t>* part_;
        std::vector<std::vector<uint8_t*>> flags_;

        // sometimes we want to compute arcflags for just 
        // a subset of nodes; those in the range [firstid, lastid)
        uint32_t firstid_;
        uint32_t lastid_;

        // values to quickly extract the flag bit for the target at hand
        uint32_t bytes_per_label_;
        uint32_t t_byte_;
        uint32_t t_bitmask_;
        
        void 
        init(warthog::graph::planar_graph* g, 
                std::vector<uint32_t>* part);

        void
        compute_down_flags(
                std::vector<uint32_t>* ids_by_rank, 
                std::vector<uint32_t>* rank);

        void
        compute_up_flags(
                std::vector<uint32_t>* ids_by_rank,
                std::vector<uint32_t>* rank);

};

}

#endif

