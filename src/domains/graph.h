#ifndef WARTHOG_GRAPH_H

// graph.h
// 
// A simple general purpose data structure for directed 
// weighted graphs. Supported input types are:
//  - warthog::gridmap objects
//  - road network graphs in the format of the 9th DIMACS competition
//
// This implemention stores all nodes and edges in two separate flat arrays
// and uses one to index the other. The graph can contain a maximum of
// 2^32 nodes and edges.
//
// @author: dharabor
// @created: 2015-01-07
//

#include <ostream>
#include <unordered_map>
#include <vector>

namespace warthog
{


class gridmap;
class graph
{
    struct node
    {
        int32_t x_;
        int32_t y_;
        uint32_t begin_;    // index of first edge
        uint32_t degree_;   // number of (incoming and outgoing) edges
    } ;

    struct edge
    {
        uint32_t head_idx_; // the node pointed to
        uint32_t wt_;       // the edge weight
    };

    public:
        graph();
        ~graph();

        // read in a grid map in the format used at the 
        // international Grid-based Path Planning Competition
        bool 
        load_grid(char* filename);

        // read in a map in the format of the 9th DIMACS
        // competition. In this format graphs are specified
        // using two files: (i) a gr file which defines edge
        // weights and endpoints and; (ii) a co file which defines 
        // node ids and planar coordinates
        bool
        load_dimacs(char*, char*);

        void
        print_dimacs(std::ostream& oss);


    private:
        void init();

        std::vector<warthog::graph::node>* nodes_;
        std::vector<warthog::graph::edge>* edges_;
        std::vector<uint32_t>* lat_lng_;
        std::string* filename_;

};

}

#endif

