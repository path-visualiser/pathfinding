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
        planar_graph(uint32_t num_nodes=0);
        planar_graph(warthog::gridmap* gm, bool store_incoming=true);
        ~planar_graph() { }

        warthog::graph::planar_graph&
        operator=(const warthog::graph::planar_graph&& other);
        
        bool
        operator==(const warthog::graph::planar_graph& other);

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
        print_dimacs_gr(std::ostream& oss, 
                uint32_t first_id, uint32_t last_id);

        // print text descriptions of the set of nodes 
        // in the range [firstid, lastid)
        // NB: if the range is invalid, nothing is printed
        void
        print_dimacs_co(std::ostream& oss, 
                uint32_t first_id, uint32_t last_id);

        // write a binary blob description of the graph
        void
        save(std::ostream& oss);
        
        // load graph from a binary blob description
        void
        load( std::istream& oss, 
                    bool store_incoming_edges = false, 
                    bool enforce_euclidean=true);

        // set the number of nodes and edges in the graph to zero.
        // any currently assigned nodes and edges have their destructors 
        // called
        void
        clear();

        // grow the graph so that the number of vertices is equal to 
        // @param num_nodes. if @param num_nodes is less thna the current
        // number of nodes, this function does nothing.
        void
        grow(uint32_t num_nodes);

        // allocate capacity for at least @param num_nodes
        // when @param num_nodes <= the number of nodes in the graph, this
        // function does nothing
        void
        capacity(uint32_t num_nodes);

        // check if arc weights are Euclidean and (optionally) fix if not.
        // "fixing" means arc weights must be at least as large as the
        // Euclidean distance between the arc's head and tail vertex.
        bool
        is_euclidean(bool fix_if_not=true);

        inline uint32_t
        get_num_nodes() const
        {
            return nodes_.size();
        }

        inline uint32_t 
        get_num_edges_out() const
        {
            uint32_t num_edges = 0;
            for(auto& node : nodes_)
            {
                num_edges += node.out_degree();
            }
            return num_edges;
        }

        inline uint32_t
        get_num_edges_in() const
        {
            uint32_t num_edges = 0;
            for(auto& node : nodes_)
            {
                num_edges += node.in_degree();
            }
            return num_edges;
        }

        // Fetch the planar coordinates of a node
        //
        // @param id: an internal graph id
        // @return x: the x coordinate of node @param id
        // @return y: the y coordinate of node @param id
        inline void
        get_xy(uint32_t id, int32_t& x, int32_t& y) const
        {
            if(id < get_num_nodes())
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
            uint32_t nodes_sz = get_num_nodes();
            if(id >= nodes_sz) { return; }
            if(xy_.size() != nodes_sz) { xy_.resize(nodes_sz); }
            xy_[id*2] = x;
            xy_[id*2+1] = y;
        }

        // Fetch a node
        //
        // @param id: an internal graph id
        // @return: the node object associated with @param id
        inline warthog::graph::node* 
        get_node(uint32_t id) 
        {
            //if(id >= ID_OFFSET && id < nodes_sz_)
            if(id < get_num_nodes())
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
        add_node(int32_t x, int32_t y, uint32_t ext_id);

        uint32_t
        add_node(int32_t x, int32_t y);

        uint32_t
        add_node();

        // print extra stuff to std::err 
        inline void 
        set_verbose(bool verbose) { verbose_ = verbose; } 

        inline bool
        get_verbose() const { return verbose_; }

        // @return the name of the file from which the 
        // current graph object was constructed
        inline const char* 
        get_filename() const { return filename_.c_str(); }

        inline size_t
        mem()
        {
            size_t mem = 0;
            for(uint32_t i = 0; i < get_num_nodes(); i++)
            {
                mem += nodes_[i].mem();
            }
            mem += sizeof(int32_t) * xy_.size() * 2;
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
        to_external_id(uint32_t in_id)  const
        { 
            if(in_id > get_num_nodes()) { return warthog::INF; }
            return id_map_.at(in_id);
        }

        // compute the proportion of bytes allocated to edges with respect
        // to the size of the address space those bytes span.
        // a value of 1 using this metric indicates that the edges perfectly 
        // fit into the allocated address space
        inline double
        edge_mem_frag() 
        {
            warthog::graph::edge *min_addr, *max_addr;
            min_addr = this->get_node(0)->outgoing_begin();
            max_addr = this->get_node(0)->outgoing_end();
            for(uint32_t i = 0; i < this->get_num_nodes(); i++)
            {
                warthog::graph::node* n = this->get_node(i);
                warthog::graph::edge* out_begin = n->outgoing_begin();
                warthog::graph::edge* out_end = n->outgoing_end();
                warthog::graph::edge* in_begin = n->incoming_begin();
                warthog::graph::edge* in_end = n->incoming_end();

                min_addr = out_begin ? 
                            (out_begin < min_addr ? out_begin : min_addr) :
                            min_addr;
                max_addr = out_end ? 
                            (out_end > max_addr ? out_end : max_addr) :
                            max_addr;
                min_addr = in_begin ? 
                            (in_begin < min_addr ? in_begin : min_addr) :
                            min_addr;
                max_addr = in_end ? 
                            (in_end > max_addr ? in_end : max_addr) :
                            max_addr;
            }
            
            uint64_t mem_lb = sizeof(warthog::graph::edge) *
                (this->get_num_edges_out() + this->get_num_edges_in());
            uint64_t mem_actual =
                (max_addr-min_addr) * sizeof(warthog::graph::edge);
            return mem_actual / (double)mem_lb;
        }

        //inline void
        //shrink_to_fit()
        //{
        //    edge* tmp = new edge[get_num_edges()];
        //    uint32_t e_index = 0;
        //    for(uint32_t i = 0; i < this->get_num_nodes(); i++)
        //    { 
        //        uint32_t in_deg = nodes_[i].in_degree();
        //        uint32_t out_deg = nodes_[i].out_degree();
        //        nodes_[i].relocate(&tmp[e_index], &tmp[e_index + in_deg]);
        //        e_index += (in_deg + out_deg);
        //    }
        //}

    private:
        // the set of nodes that comprise the graph
        std::vector<warthog::graph::node> nodes_;
        
        // planar coordinates stored as adjacent pairs (x, then y)
        std::vector<int32_t> xy_;

        // these containers serve to map from external graph ids to 
        // internal graph ids
        std::vector<uint32_t> id_map_; 
        std::unordered_map<uint32_t, uint32_t> ext_id_map_; 

        std::string filename_;
        bool verbose_;

        bool
        grid2graph(warthog::gridmap*, bool);

};

}
}

#endif

