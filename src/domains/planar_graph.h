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

#include <ostream>
#include <unordered_map>
#include <vector>

namespace warthog
{

class gridmap;
class planar_graph
{
    public:
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

        planar_graph();
        ~planar_graph();

        // read in a grid map in the format used at the 
        // international Grid-based Path Planning Competition
        bool 
        load_grid(const char* filename);

        // read in a map in the format of the 9th DIMACS
        // competition. In this format graphs are specified
        // using two files: (i) a gr file which defines edge
        // weights and endpoints and; (ii) a co file which defines 
        // node ids and planar coordinates
        // @param backward is used to control arc direction.
        // setting this to true swaps the head and tail of every
        // arc (the default is false)
        bool
        load_dimacs(const char*, const char*, bool backward=false);

        void
        print_dimacs(std::ostream& oss);

        inline uint32_t
        get_num_nodes()
        {
            return nodes_->size();
        }

        inline uint32_t 
        get_num_edges()
        {
            return edges_->size();
        }

        inline void
        get_xy(uint32_t id, int32_t& x, int32_t& y)
        {
            if(id < nodes_->size()) 
            {
                warthog::planar_graph::node n = nodes_->at(id);
                x = n.x_;
                y = n.y_; 
            }
        }

        inline warthog::planar_graph::node
        get_node(uint32_t id)
        {
           return nodes_->at(id) ;
        }

        inline warthog::planar_graph::edge
        get_edge(uint32_t id)
        {
            return edges_->at(id);
        }

        size_t
        mem()
        {
            return 
                sizeof(warthog::planar_graph::node)*nodes_->size() +
                sizeof(warthog::planar_graph::edge)*edges_->size() +
                sizeof(char)*filename_->size() +
                sizeof(*this);
        }


    private:
        void init();

        std::vector<warthog::planar_graph::node>* nodes_;
        std::vector<warthog::planar_graph::edge>* edges_;
        std::string* filename_;
};

}

#endif

