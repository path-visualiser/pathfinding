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

#include <vector>
#include <cassert>
#include <cstdint>
#include <iostream>

namespace warthog
{

namespace graph
{

class planar_graph;

}

class problem_instance;
class search_node;

class af_filter
{
    public:
        // @param g: the search graph
        // @param part: a partitioning of the graph nodes
        af_filter(
                warthog::graph::planar_graph* g, 
                std::vector<uint32_t>* part);

        af_filter(
                warthog::graph::planar_graph* g, 
                std::vector<uint32_t>* part,
                const char* arcflags_file);

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

        void
        set_target(uint32_t goal_id) 
        { 
            uint32_t t_part = part_->at(goal_id);
            t_byte_ = t_part >> 3;
            t_bitmask_ = 1 << (t_part & 7);
        }

        void
        compute();

        void
        compute(uint32_t firstid, uint32_t lastid);

        void
        compute_ch(std::vector<uint32_t>* rank);

        void
        compute_ch(uint32_t firstid, uint32_t lastid, std::vector<uint32_t>* rank);

        void
        print(std::ostream& out);

        bool
        load_arcflags_file(const char* filename);

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

};

}

#endif

