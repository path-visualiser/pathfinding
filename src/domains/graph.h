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

#include <unordered_map>
#include <vector>

namespace warthog
{

// nodes store only enough information to identify
// where the list of incoming and outgoing edges
// begin in the edges array. the index of the last
// edge is not kept and we rely instead on detecting a 
// terminator character.
struct node
{
    uint32_t edge_idx_; // first edge
    uint32_t degree_;   // number of (incoming and outgoing) edges
} ;

struct edge
{
    uint32_t node_idx_; // the (other) endpoint of this edge
    uint32_t wt_;       // the edge weight
};

class gridmap;
class graph
{
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

    private:
        void init();

        std::vector<warthog::node>* nodes_;
        std::vector<warthog::edge>* edges_;
        std::vector<uint32_t>* lat_lng_;
        std::string* filename_;

};

}

#endif


