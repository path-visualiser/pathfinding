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
// @created: 2016-01-07
//

#include "constants.h"
#include "graph.h"

#include <ostream>
#include <unordered_map>
#include <vector>

namespace warthog
{
class gridmap;

namespace graph
{
class planar_graph
{

    public:
        planar_graph();
        planar_graph(warthog::gridmap* gm, bool store_incoming=true);
        planar_graph(warthog::graph::planar_graph& other);
        ~planar_graph();

        warthog::graph::planar_graph&
        operator=(const warthog::graph::planar_graph& other);

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

        // print text descriptions of the set of arcs associated with 
        // all nodes in the range [firstid, lastid)
        // NB: if the range is invalid, nothing is printed
        void
        print_dimacs_gr(std::ostream& oss, uint32_t first_id, uint32_t last_id);

        // print text descriptions of the set of nodes 
        // in the range [firstid, lastid)
        // NB: if the range is invalid, nothing is printed
        void
        print_dimacs_co(std::ostream& oss, uint32_t first_id, uint32_t last_id);

        inline uint32_t
        get_num_nodes()
        {
            return nodes_sz_;
        }

        inline uint32_t 
        get_num_edges()
        {
            uint32_t num_edges = 0;
            for(uint32_t i = 0; i < nodes_sz_; i++)
            {
                num_edges += nodes_[i].out_degree();
            }
            return num_edges;
        }

        // Fetch the planar coordinates of a node
        //
        // @param id: an internal graph id
        // @return x: the x coordinate of node @param id
        // @return y: the y coordinate of node @param id
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

        // Set the planar coordinates of a node
        //
        // @param id: an internal graph id
        // @param x: the x coordinate of node @param id
        // @param y: the y coordinate of node @param id
        inline void
        set_xy(uint32_t id, int32_t x, int32_t y)
        {
            if(id < nodes_sz_)
            {
                xy_[id*2] = x;
                xy_[id*2+1] = y;
            }
        }

        // Fetch a node
        //
        // @param id: an internal graph id
        // @return: the node object associated with @param id
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

        // Add a new node into the graph. If a node already exists in the 
        // graph with the same external id as @param ext_id then then nothing 
        // is added.
        // NB: a new node is always added if @param ext_id is equal to the 
        // value warthog::INF
        //
        // @param x: the x-coordinate of the new node
        // @param y: the y-coordinate of the new node
        // @param ext_id: an (optional) external id for this node
        // @return: the internal graph id of the new node or the id of the
        // existing node whose graph id is equal to @param ext_id
        uint32_t
        add_node(int32_t x, int32_t y, uint32_t ext_id = warthog::INF);

        // print extra stuff to std::err 
        inline void 
        set_verbose(bool verbose) { verbose_ = verbose; } 

        inline bool
        get_verbose() { return verbose_; }

        // @return the name of the file from which the 
        // current graph object was constructed
        inline const char* 
        get_filename() { return filename_.c_str(); }

        // increase the node capacity of the graph
        // (NB: capacity is not the same as size!)
        inline size_t 
        reserve(uint32_t new_cap) { return resize(new_cap); }

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
        
        // convert an external node id (e.g. as it appears in an input file)
        // to the equivalent internal id used by the current graph
        //
        // @param ex_id: the external id of the node
        // @return: the corresponding internal id for @param ex_id
        // if ex_id did not appear in the input file (or if the graph
        // was not created from an input file) the function returns
        // the value warthog::INF
        inline uint32_t 
        to_graph_id(uint32_t ext_id) 
        { 
            if(id_map_.size() == 0) { return warthog::INF; }

        //    uint32_t min = 0;
        //    uint32_t max = id_map_.size()-1;
        //    uint32_t pos = max >> 1;
        //    while(min <= max)
        //    {
        //        uint32_t current = id_map_.at(pos);
        //        if(current == ex_id) { return pos; }
        //        if(ex_id < current) { max = pos-1; }
        //        else { min = pos+1; }
        //        pos = (max+min) >> 1;
        //    } 
        //    return warthog::INF;

            std::unordered_map<uint32_t, uint32_t>::iterator it 
                    = ext_id_map_.find(ext_id);

            if(it == ext_id_map_.end()) { return warthog::INF; }
            return (*it).second;
        }

        // convert an internal node id (i.e. as used by the current graph
        // to the equivalent external id (e.g. as appears in an input file)
        //
        // @param ex_id: the external id of the node
        // @return: the corresponding internal id for @param ex_id
        // if ex_id did not appear in the input file (or if the graph
        // was not created from an input file) the function returns
        // the value warthog::INF
        inline uint32_t 
        to_external_id(uint32_t in_id) 
        { 
            if(in_id > nodes_sz_) { return warthog::INF; }
            return id_map_.at(in_id);
        }


    private:
        std::string filename_;

        // these containers serve to map from external graph ids to 
        // internal graph ids
        std::vector<uint32_t> id_map_; 
        std::unordered_map<uint32_t, uint32_t> ext_id_map_; 
        
        // the set of nodes that comprise the graph
        uint32_t nodes_sz_;
        uint32_t nodes_cap_;
        warthog::graph::node* nodes_;

        // planar coordinates stored as adjacent pairs (x, then y)
        int32_t* xy_;

        bool verbose_;

        size_t
        resize(uint32_t new_cap);

        bool
        grid2graph(warthog::gridmap*, bool);


};

}
}

#endif

