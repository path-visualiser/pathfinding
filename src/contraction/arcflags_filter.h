#ifndef WARTHOG_ARCFLAGS_FILTER_H
#define WARTHOG_ARCFLAGS_FILTER_H

// arcflags_filter.h
//
// Given a partitioning of the graph, this class adds labels to every edge
// indicating whether the edge appears on any optimal path from its tail node
// to a any node in each partition. Each edge label is stored as a bitfield.
// One bit per partition. The labels are used to filter nodes during search
// as follows: if the edge used to reach a candidate node is not labeled
// with a 1-value for the partition containing the goal node then the edge 
// cannot possibly appear on any optimal path to the goal.
//
// For theoretical details see:
//
// [Fast Point-to-Point Shortest Path Computations with Arc-Flags,
//  Ekkehard Ko ̈hler, Rolf H. Mo ̈hring, and Heiko Schilling,
//  In The Shortest Path Problem: Ninth DIMACS Implementation Challenge, 
//  Edited by Demetrescu, Camil and Goldberg, Andrew V. and Johnson, David S 
//  pp 41-72, American Mathematical Society, 2009]
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

class arcflags_filter
{
    public:
        // @param g: the search graph
        // @param part: a partitioning of the graph nodes
        arcflags_filter(
                warthog::graph::planar_graph* g, 
                std::vector<uint32_t>* rank,
                std::vector<uint32_t>* part);

        arcflags_filter(
                warthog::graph::planar_graph* g, 
                std::vector<uint32_t>* rank,
                std::vector<uint32_t>* part,
                char* arcflags_file);

        ~arcflags_filter();

        inline bool 
        filter(uint32_t node_id, uint32_t edge_idx)
        {
            uint8_t* label = flags_.at(node_id).at(edge_idx);
            return label[t_byte_] & t_bitmask_;
        }

        void
        set_instance(warthog::problem_instance* instance);

        void
        compute();

        void
        compute(uint32_t firstid, uint32_t lastid);

        void
        print(std::ostream& out);

    private:    
        uint32_t nparts_;
        warthog::graph::planar_graph* g_;
        std::vector<uint32_t>* part_;
        std::vector<uint32_t>* rank_;
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
                std::vector<uint32_t>* rank,
                std::vector<uint32_t>* part);

        void
        load_arcflags_file(const char* filename);
};

}

#endif

