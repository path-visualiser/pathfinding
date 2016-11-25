#ifndef WARTHOG_PLANAR_GRAPH_H
#define WARTHOG_PLANAR_GRAPH_H

// planar_graph.h
// 
// A simple general purpose data structure for directed 
// weighted planar graphs. Supported input types are:
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

#include "constants.h"
#include "graph.h"

#include <ostream>
#include <unordered_map>
#include <vector>

namespace warthog
{
namespace graph
{

class gridmap;
class planar_graph
{

    public:
        planar_graph();
        planar_graph(warthog::graph::planar_graph& other);
        ~planar_graph();

        // read in a grid map in the format used at the 
        // international Grid-based Path Planning Competition
        bool 
        load_grid(const char* filename, bool store_incoming=true);

        // read in a map in the format of the 9th DIMACS
        // competition. In this format graphs are specified
        // using two files: (i) a gr file which defines edge
        // weights and endpoints and; (ii) a co file which defines 
        // node ids and planar coordinates
        //
        // @param reverse_arcs: reverses the direction of each edge
        // @param duplicate_edges: store edges with both head and tail node
        // @param enforce_euclidean: arc lengths must be >= euclidean distance
        bool
        load_dimacs(const char*, const char*, 
                bool reverse_arcs =false, 
                bool store_incoming_edges = false,
                bool enforce_euclidean=true);

        void
        print_dimacs_gr(std::ostream& oss);

        void
        print_dimacs_co(std::ostream& oss);

        // returns number of nodes + ID_OFFSET number of
        // dummy elements located at the start of the nodes array
        inline uint32_t
        get_num_nodes()
        {
            return nodes_sz_;
        }

        // returns the number of dummy nodes at the start 
        // of the nodes array
        inline uint32_t
        id_offset() { return ID_OFFSET; }

        inline uint32_t 
        get_num_edges()
        {
            uint32_t num_edges = 0;
            for(uint32_t i = ID_OFFSET; i < nodes_sz_; i++)
            {
                num_edges += nodes_[i].out_degree();
            }
            return num_edges;
        }

        inline void
        get_xy(uint32_t id, int32_t& x, int32_t& y)
        {
            if(id < nodes_sz_)
            {
                x = xy_[id*2];
                y = xy_[id*2+1];
            }
            else
            {
                x = y = warthog::INF;
            }
        }

        inline void
        set_xy(uint32_t id, int32_t x, int32_t y)
        {
            if(id < nodes_sz_)
            {
                xy_[id*2] = x;
                xy_[id*2+1] = y;
            }
        }

        inline warthog::graph::node* 
        get_node(uint32_t id)
        {
            //if(id >= ID_OFFSET && id < nodes_sz_)
            if(id < nodes_sz_)
            {
                return &nodes_[id];
            }
            return 0;
        }

        inline uint32_t
        add_node(int32_t x, int32_t y)
        {
            if(nodes_cap_ == nodes_sz_)
            {
                resize(nodes_cap_ == 0 ? 1 : (nodes_cap_*2));
            }
            uint32_t index = nodes_sz_;
            nodes_sz_++;
            xy_[index*2] = x;
            xy_[index*2+1] = y;
            return index;
        }

        inline void 
        set_verbose(bool verbose) { verbose_ = verbose; } 

        inline bool
        get_verbose() { return verbose_; }

        inline const char* 
        get_filename() { return filename_.c_str(); }

        // increase the node capacity of the graph
        // (NB: capacity is not the same as size!)
        bool
        reserve(uint32_t new_cap);

        inline size_t
        mem()
        {
            size_t mem = 0;
            for(uint32_t i = 0; i < nodes_sz_; i++)
            {
                mem += nodes_[i].mem();
            }
            mem += sizeof(*xy_) * nodes_sz_ * 2;
            mem += sizeof(char)*filename_.length() +
                sizeof(*this);
            return mem;
        }


    private:
        std::string filename_;
        
        // the set of nodes that comprise the graph
        uint32_t nodes_sz_;
        uint32_t nodes_cap_;
        warthog::graph::node* nodes_;

        // planar coordinates stored as adjacent pairs (x, then y)
        int32_t* xy_;

        bool verbose_;

        // DIMACS-fomat node ids are specified as 1-indexed. Our internal
        // representation is zero indexed. We convert between them by 
        // adding an offset
        uint32_t ID_OFFSET; 

        size_t
        resize(uint32_t new_cap);

};

}
}

#endif

