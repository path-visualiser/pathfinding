#ifndef WARTHOG_AF_FILTER_H
#define WARTHOG_AF_FILTER_H

// af_filter.h
//
// An arcflags filter. Works on regular graphs or contraction hierarchies.
// The contraction hierarchy specialisation allows faster preprocessing by 
// limiting the scope of each search; e.g. once a node is reached by following 
// a down-arc the only successors which are generated and relaxed must also 
// be reached by a down-arc.
//
// For theoretical details relating to the basic arcflags technique see:
//
// [Fast Point-to-Point Shortest Path Computations with Arc-Flags,
//  Ekkehard Ko ̈hler, Rolf H. Mo ̈hring, and Heiko Schilling,
//  In The Shortest Path Problem: Ninth DIMACS Implementation Challenge, 
//  Edited by Demetrescu, Camil and Goldberg, Andrew V. and Johnson, David S 
//  pp 41-72, American Mathematical Society, 2009]
//
// For theoretical details on contraction hierarchies see:
//
// [Geisbergerger, Sanders, Schultes and Delling. 
// Contraction Hierarchies: Faster and Simpler Hierarchical 
// Routing in Road Networks. In Proceedings of the 2008
// Workshop on Experimental Algorithms (WEA)]
//
// @author: dharabor
// @created: 2016-08-16
//

#include "arclabels.h"

#include <vector>
#include <cassert>
#include <cstdint>
#include <iostream>

namespace warthog
{

class problem_instance;
class search_node;

namespace graph
{

class planar_graph;

}

class af_filter
{
    public:
        // @param g: the search graph
        // @param part: a partitioning of the graph nodes
        af_filter(std::vector<uint32_t>* part);

        ~af_filter();

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

        inline void
        set_target(uint32_t goal_id) 
        { 
            uint32_t t_part = part_->at(goal_id);
            t_byte_ = t_part >> 3;
            t_bitmask_ = 1 << (t_part & 7);
        }

        inline uint8_t*
        get_label(uint32_t node_id, uint32_t edge_index)
        { return flags_.at(node_id).at(edge_index); }
        
        inline uint32_t
        get_partition(uint32_t node_id)
        {
            return part_->at(node_id);
        }

        inline void
        print(std::ostream& out);

        bool
        load_labels(const char* filename, warthog::graph::planar_graph* g);

    private:    
        std::vector<uint32_t>* part_;
        std::vector<std::vector<uint8_t*>> flags_;
        warthog::arclabels::af_params par_;
        warthog::graph::planar_graph* g_;

        // we cache which bit corresponds to the target partition
        uint32_t t_byte_;
        uint32_t t_bitmask_;
};

}

#endif

