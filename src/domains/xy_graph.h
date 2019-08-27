#ifndef WARTHOG_XY_GRAPH_H
#define WARTHOG_XY_GRAPH_H

// xy_graph.h
// 
// A simple general purpose data structure for directed weighted xy graphs 
// (i.e. spatial networks embedded in two dimensions). 
// Supported input types are:
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


#include "dimacs_parser.h"
#include "euclidean_heuristic.h"
#include "forward.h"
#include "graph.h"
#include "gridmap_expansion_policy.h"
#include "constants.h"

#include <ostream>
#include <unordered_map>
#include <vector>

namespace warthog
{

namespace graph
{

template<class T_NODE, class T_EDGE>
class xy_graph_base
{
    public:
        // create an empty graph
        xy_graph_base(uint32_t num_nodes=0) : verbose_(false)
        {
            grow(num_nodes);
        }

        ~xy_graph_base() 
        { }

        // create a graph from a grid map. 
        bool
        load_from_grid(warthog::gridmap* gm, bool store_incoming=true)
        {
            std::vector<uint32_t> id_map(gm->header_height() * gm->header_width());

            // add each traversable tile as a node in the graph
            uint32_t num_traversable_tiles = 0;
            warthog::gridmap_expansion_policy exp(gm);
            for(uint32_t y = 0; y < gm->header_height(); y++)
            {
                for(uint32_t x = 0; x < gm->header_width(); x++)
                {
                    // we differentiate between external grid ids
                    // (unpadded) and internal grid ids (with padding)
                    uint32_t from_gm_id = y * gm->header_width() + x;

                    // skip obstacles
                    if(!gm->get_label(gm->to_padded_id(from_gm_id))) 
                    { continue; }
                    num_traversable_tiles++;

                    // add graph node (we scale up all costs and coordinates)
                    id_map[from_gm_id] = add_node(
                        (int32_t)(x * warthog::graph::GRID_TO_GRAPH_SCALE_FACTOR),
                        (int32_t)(y * warthog::graph::GRID_TO_GRAPH_SCALE_FACTOR),
                        from_gm_id);
                }
            }
            assert(num_traversable_tiles == get_num_nodes());

            // add edges
            for(uint32_t y = 0; y < gm->header_height(); y++)
            {
                for(uint32_t x = 0; x < gm->header_width(); x++)
                {
                    // we use a grid-expansion policy here which enforces
                    // standard constraints on grid moves re; e.g. corner cutting
                    uint32_t from_graph_id;
                    {
                        from_graph_id = id_map[y*gm->header_width() + x];
                        assert(from_graph_id < num_traversable_tiles);
                    }

                    warthog::search_node* nei = 0;
                    warthog::search_node* n = 0;
                    edge_cost_t edge_cost = 0;

                    n = exp.generate(gm->to_padded_id(y*gm->header_width() + x));
                    exp.expand(n, 0);
                    for(exp.first(nei, edge_cost); nei != 0; exp.next(nei, edge_cost))
                    {
                        uint32_t to_graph_id;
                        {
                            uint32_t nei_x, nei_y;
                            gm->to_unpadded_xy((uint32_t)nei->get_id(), nei_x, nei_y);
                            uint32_t to_gm_id = nei_y * gm->header_width() + nei_x;
                            to_graph_id = id_map[to_gm_id];
                            assert(from_graph_id != to_graph_id);
                            assert(to_graph_id < num_traversable_tiles);
                        }

                        T_NODE* gr_from = get_node(from_graph_id);
                        T_NODE* gr_to = get_node(to_graph_id);
                        edge_cost_t gr_weight = edge_cost * warthog::graph::GRID_TO_GRAPH_SCALE_FACTOR;

                        gr_from->add_outgoing(T_EDGE(to_graph_id, gr_weight));
                        if(store_incoming)
                        {
                            gr_to->add_incoming(T_EDGE(from_graph_id, gr_weight));
                        }
                    }
                }
            }
            return true;
        }

        // read in a map in the format of the 9th DIMACS
        // competition. In this format graphs are specified
        // using two files: (i) a gr file which defines edge
        // weights and endpoints and; (ii) a co file which defines 
        // node ids and xy coordinates
        //
        // @param reverse_arcs: reverses the direction of each edge
        // @param duplicate_edges: store edges with both head and tail node
        // @param enforce_euclidean: arc lengths must be >= euclidean distance
        bool
        load_from_dimacs(const char* gr_file, const char* co_file,
                bool reverse_arcs =false,
                bool store_incoming_edges = false,
                bool enforce_euclidean=true)
        {
            clear();

            //warthog::dimacs_parser dimacs(gr_file, co_file);
            warthog::dimacs_parser dimacs;
            if(!( dimacs.load_graph(gr_file) && dimacs.load_graph(co_file) ))
            {
                return false;
            }
            //std::cout << "loaded" << std::endl;
            
            filename_.assign(gr_file);

            // sort nodes by their ids
            bool (*node_comparator) 
                 (warthog::dimacs_parser::node, warthog::dimacs_parser::node) =
                    [](warthog::dimacs_parser::node n1,
                       warthog::dimacs_parser::node n2) -> bool
                        {
                            return n1.id_ < n2.id_;
                        };
            std::sort(dimacs.nodes_begin(), dimacs.nodes_end(), node_comparator);
            if(verbose_) { std::cerr << "nodes, sorted" << std::endl; }
            
            // allocate memory for nodes
            uint32_t num_nodes_dimacs = dimacs.get_num_nodes();
            capacity(num_nodes_dimacs);
            for(warthog::dimacs_parser::node_iterator it = dimacs.nodes_begin();
                    it != dimacs.nodes_end(); it++)
            {
               // NB: xy_graph uses 0-indexed ids, DIMACS uses 1-indexed ids
               // here, we update the id as we add it to the graph
               add_node((*it).x_, (*it).y_, (*it).id_ - 1);
            }
            if(verbose_) { std::cerr << "nodes, converted" << std::endl; }

            // scan the list of edges so as to know the in and out degree of each node
            std::vector<warthog::graph::ECAP_T> in_deg(dimacs.get_num_nodes(), 0);
            std::vector<warthog::graph::ECAP_T> out_deg(dimacs.get_num_nodes(), 0);
            for(warthog::dimacs_parser::edge_iterator it = dimacs.edges_begin();
                    it != dimacs.edges_end(); it++)
            {
               // NB: xy_graph uses 0-indexed ids, DIMACS uses 1-indexed ids
               // here, we convert to 0-indexing on-the-fly
                in_deg[it->head_id_ - 1]++;
                out_deg[it->tail_id_ - 1]++;
            }

            // allocate memory for edges
            for(warthog::dimacs_parser::node_iterator it = dimacs.nodes_begin();
                    it != dimacs.nodes_end(); it++)
            {
               // NB: xy_graph uses 0-indexed ids, DIMACS uses 1-indexed ids
               // here, we convert to 0-indexing on-the-fly
                uint32_t nid = (*it).id_ - 1;
                nodes_[nid].capacity( 
                    store_incoming_edges ? in_deg[nid] : 0, out_deg[nid]);
            }

            // convert edges to graph format
            for(warthog::dimacs_parser::edge_iterator it = dimacs.edges_begin();
                    it != dimacs.edges_end(); it++)
            {
               // NB: xy_graph uses 0-indexed ids, DIMACS uses 1-indexed ids
               // here, we convert to 0-indexing on-the-fly
                uint32_t hid = (*it).head_id_ - 1;
                uint32_t tid = (*it).tail_id_ - 1;

                // in a reverse graph the head and tail of every 
                // edge are swapped
                if(reverse_arcs)
                {
                    uint32_t tmp = hid;
                    hid = tid;
                    tid = tmp;
                }
                warthog::graph::edge e;
                e.node_id_ = hid;
                e.wt_ = (*it).weight_;

#ifndef NDEBUG
                warthog::graph::ECAP_T deg_before = nodes_[tid].out_degree();
#endif

                nodes_[tid].add_outgoing(e);

#ifndef NDEBUG
                warthog::graph::ECAP_T deg_after = nodes_[tid].out_degree();
#endif

                // edges can be stored twice: once as an incoming edge 
                // and once as an outgoing edge
                if(store_incoming_edges)
                {
                    e.node_id_ = tid;
                    nodes_[hid].add_incoming(e);

#ifndef NDEBUG
                    // sanity check: either the edge was a duplicate (and nothing
                    // was added) or the added was added in which case it should
                    // be the last element in the outgoing list of the tail node
                    warthog::graph::edge sanity = 
                        *(nodes_[tid].outgoing_end()-1);
                    assert(deg_before == deg_after || 
                            sanity.node_id_ == hid);
#endif
                }
                if(verbose_ && ((it - dimacs.edges_begin()) % 1000) == 0)
                {
                    std::cerr 
                        << "\rconverted K edges " 
                        << (it - dimacs.edges_begin()) / 1000;
                }
            }
            is_euclidean(enforce_euclidean);
            if(verbose_) { std::cout << "edges, converted" << std::endl; }

            std::cerr << "edge memory fragmentation: (1=none): " 
                << this->edge_mem_frag() << std::endl;
            //print_dimacs_co(std::cerr, 0, (uint32_t)nodes_.size()-1);
            //print_dimacs_gr(std::cerr, 0, (uint32_t)nodes_.size()-1);
            return true;
        }

        // load graph from a binary blob description
        bool
        load_from_blob(  std::istream& oss,
                           bool store_incoming_edges = false,
                           bool enforce_euclidean=true)
        {
            // read number of nodes 
            uint32_t num_nodes;
            oss.read((char*)&num_nodes, sizeof(num_nodes));
            if(num_nodes > 0)
            {
                nodes_.resize(num_nodes);
                if(!oss.good()) { std::cerr << "err; reading #nodes data\n"; return false; }
            }

            // ready xy coordinates
            uint32_t xy_sz;
            oss.read((char*)&xy_sz, sizeof(xy_sz));
            if(!oss.good())
            { std::cerr << "err; reading #xy data\n"; return false;}
            if(xy_sz > 0)
            {
                xy_.resize(xy_sz);
                oss.read((char*)&xy_[0], sizeof(int32_t)*xy_sz);
                if(!oss.good()) { std::cerr << "err; reading xy data\n"; return false; }
            }

            // read external ids
            uint32_t id_map_sz;
            oss.read((char*)&id_map_sz, sizeof(id_map_sz));
            if(!oss.good()) { std::cerr << "err; reading #external-ids\n"; return false; }
            if(id_map_sz == get_num_nodes())
            {
                id_map_.resize(id_map_sz);
                oss.read((char*)&id_map_[0], sizeof(int32_t)*id_map_sz);
                if(!oss.good()) 
                { std::cerr << "err; reading external id data\n"; return false;}

                // build hashtable for reverse lookups
                ext_id_map_.reserve((size_t)get_num_nodes());
                for(uint32_t i = 0; i < id_map_.size(); i++)
                {
                    ext_id_map_.insert(std::pair<uint32_t, uint32_t>(id_map_[i], i));
                }
            }
            else
            {
                std::cerr 
                    << "warn; ignoring external id values (read " << id_map_sz 
                    << " ids but graph has " << get_num_nodes() << " nodes)" 
                    << std::endl;
            }

            // read node and edge data
            uint32_t nid = 0;
            for( ; nid < num_nodes; nid++)
            {
                nodes_[nid].load(oss);
                if(!oss.good()) 
                { std::cerr 
                    << "err; reading outgoing adj list for node " << nid
                    << std::endl;
                    break;
                }
            }
            if(nid != num_nodes)
            {
                std::cerr 
                    << "err; read " << nid << " of " << get_num_nodes()
                    << " adjacency lists. missing data? "
                    << std::endl;
            }

            if(store_incoming_edges)
            {
                std::cerr << "#outgoing " << this->get_num_edges_out() << std::endl;

                // collect all edges and compute in degree of every node
                std::vector<warthog::graph::ECAP_T> in_deg_all(get_num_nodes(), 0);
                std::vector<uint32_t> from_id;
                std::vector<uint32_t> to_id;
                std::vector<edge_cost_t> weight;

                for( nid = 0; nid < num_nodes; nid++)
                {
                    for( warthog::graph::edge_iter it = nodes_[nid].outgoing_begin();
                            it != nodes_[nid].outgoing_end(); 
                            it++ ) 
                    { 
                        assert(in_deg_all[it->node_id_] < warthog::graph::ECAP_MAX);
                        in_deg_all[it->node_id_]++;
                        from_id.push_back(nid);
                        to_id.push_back(it->node_id_);
                        weight.push_back(it->wt_);
                    }
                }

                // allocate memory
                for( nid = 0; nid < get_num_nodes(); nid++)
                { nodes_[nid].capacity(in_deg_all[nid], nodes_[nid].out_degree()); }

                // add incoming edges 
                for(uint32_t i = 0; i < from_id.size(); i++)
                {
                    nodes_[to_id[i]].add_incoming(
                            warthog::graph::edge(from_id[i], weight[i]));
                }
                std::cerr << "#incoming " << this->get_num_edges_in() << std::endl;
            }

            return is_euclidean(enforce_euclidean);
        }

        warthog::graph::xy_graph_base<T_NODE, T_EDGE>&
        operator=(const warthog::graph::xy_graph_base<T_NODE, T_EDGE>&& other)
        {
            // my memory, not want
            clear();

            // your memory, can has?
            verbose_ = other.verbose_;
            filename_ = other.filename_;
            nodes_ = std::move(other.nodes_);
            xy_ = std::move(other.xy_);
            id_map_ = std::move(other.id_map_);
            ext_id_map_ = std::move(other.ext_id_map_);

            return *this;
        }
        
        bool
        operator==(const warthog::graph::xy_graph_base<T_NODE, T_EDGE>& other)
        {
            if(xy_.size() == other.xy_.size())
            {
                for(uint32_t i = 0; i < xy_.size(); i++)
                {
                    if(xy_[i] != other.xy_[i]) { return false; }
                }
            } 
            else { return false; }


            if(id_map_.size() == other.id_map_.size())
            {
                for(uint32_t i = 0; i < id_map_.size(); i++)
                {
                    if(id_map_[i] != other.id_map_[i]) { return false; }
                }
            }
            else { return false; }
            
            if(get_num_nodes() == other.get_num_nodes())
            {
                for(uint32_t i = 0; i < get_num_nodes(); i++)
                {
                    if(!(nodes_[i] == other.nodes_[i])) { return false; }
                }
                return true;
            }
            return false;
        }



        // set the number of nodes and edges in the graph to zero.
        // any currently assigned nodes and edges have their destructors 
        // called
        void
        clear()
        {
            nodes_.clear();
            xy_.clear();
            id_map_.clear();
            ext_id_map_.clear();
        }

        // grow the graph so that the number of vertices is equal to 
        // @param num_nodes. if @param num_nodes is less thna the current
        // number of nodes, this function does nothing.
        void
        grow(size_t num_nodes)
        {
            if(num_nodes <= nodes_.size()) { return; }
            nodes_.resize(num_nodes);
            xy_.resize(num_nodes*2);
            id_map_.resize(num_nodes);
            ext_id_map_.reserve(num_nodes);
        }

        // allocate capacity for at least @param num_nodes
        // when @param num_nodes <= the number of nodes in the graph, this
        // function does nothing
        void
        capacity(size_t num_nodes)
        {
            nodes_.reserve(num_nodes);
            xy_.reserve(num_nodes*2);
            id_map_.reserve(num_nodes);
            ext_id_map_.reserve(num_nodes);
        }


        inline uint32_t
        get_num_nodes() const
        {
            return (uint32_t)nodes_.size();
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

        // Fetch the xy coordinates of a node
        //
        // @param id: an internal graph id
        // @return x: the x coordinate of node @param id
        // @return y: the y coordinate of node @param id
        inline void
        get_xy(uint32_t id, int32_t& x, int32_t& y) const
        {
            assert(id < (xy_.size()>>1) );
            x = xy_[id*2];
            y = xy_[id*2+1];
        }

        // Set the xy coordinates of a node
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
        inline T_NODE*
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
        add_node(int32_t x, int32_t y, uint32_t ext_id)
        {
            // check if a node with the same external id already exists; if so 
            // return the id of the existing node. otherwise, add a new node
            //uint32_t graph_id = to_graph_id(ext_id);
            //if(graph_id != warthog::INF32) { return graph_id; } 

            uint32_t graph_id = (uint32_t)get_num_nodes();
            nodes_.push_back(warthog::graph::node());
            xy_.push_back(x);
            xy_.push_back(y);
            //id_map_.push_back(ext_id);
            //if(ext_id != warthog::INF32)
            //{
            //    ext_id_map_.insert(std::pair<uint32_t, uint32_t>(ext_id, graph_id));
            //}
            return graph_id;
        }

        uint32_t
        add_node(int32_t x, int32_t y)
        {
            return add_node(x, y, warthog::INF32);
        }

        uint32_t
        add_node()
        {
            return add_node(warthog::INF32, warthog::INF32, warthog::INF32);
        }

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
        // the value warthog::INF32
        inline uint32_t 
        to_graph_id(uint32_t ext_id) 
        { 
//            assert(ext_id < nodes_.size());
            return ext_id;
//            if(id_map_.size() == 0) { return warthog::INF32; }
//
//            std::unordered_map<uint32_t, uint32_t>::iterator it 
//                    = ext_id_map_.find(ext_id);
//
//            if(it == ext_id_map_.end()) { return warthog::INF32; }
//            return (*it).second;
      }

        // convert an internal node id (i.e. as used by the current graph
        // to the equivalent external id (e.g. as appears in an input file)
        //
        // @param ex_id: the external id of the node
        // @return: the corresponding internal id for @param ex_id
        // if ex_id did not appear in the input file (or if the graph
        // was not created from an input file) the function returns
        // the value warthog::INF32
        inline uint32_t 
        to_external_id(uint32_t in_id)  const
        { 
            //assert(in_id < nodes_.size());
            return in_id;
            //if(in_id > get_num_nodes()) { return warthog::INF32; }
            //return id_map_.at(in_id);
        }

        // compute the proportion of bytes allocated to edges with respect
        // to the size of the address space those bytes span.
        // a value of 1 using this metric indicates that the edges perfectly 
        // fit into the allocated address space
        inline double
        edge_mem_frag() 
        {
            T_EDGE *min_addr, *max_addr;
            min_addr = this->get_node(0)->outgoing_begin();
            max_addr = this->get_node(0)->outgoing_end();
            for(uint32_t i = 0; i < this->get_num_nodes(); i++)
            {
                T_NODE* n = this->get_node(i);
                T_EDGE* out_begin = n->outgoing_begin();
                T_EDGE* out_end = n->outgoing_end();
                T_EDGE* in_begin = n->incoming_begin();
                T_EDGE* in_end = n->incoming_end();

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
            
            size_t mem_lb = sizeof(warthog::graph::edge) *
                (this->get_num_edges_out() + this->get_num_edges_in());
            size_t mem_actual =
                ((size_t)(max_addr-min_addr)) * sizeof(warthog::graph::edge);
            return mem_actual / (double)mem_lb;
        }

        // write a binary blob description of the graph
        void
        save(std::ostream& oss)
        {
            // write number of nodes
            size_t nodes_sz = get_num_nodes();
            oss.write((char*)(&nodes_sz), sizeof(nodes_sz));
            if(!oss.good()) 
            { std::cerr << "err writing #nodes\n"; return; }

            // write xy values
            size_t xy_sz = xy_.size();
            oss.write((char*)&xy_sz, sizeof(xy_sz));
            if(!oss.good()) 
            { std::cerr << "err writing #xy\n"; return; }
            oss.write((char*)&xy_[0], (long)(sizeof(int32_t)*xy_sz));
            if(!oss.good()) 
            { std::cerr << "err writing xy coordinates data\n"; return; }

            // write external ids
            size_t id_map_sz = id_map_.size();
            oss.write((char*)&id_map_sz, sizeof(int32_t));
            if(!oss.good()) 
            { std::cerr << "err writing #external ids\n"; return; }
            oss.write((char*)&id_map_[0], (long)(sizeof(int32_t)*id_map_sz));
            if(!oss.good()) 
            { std::cerr << "err writing external ids\n"; return; }
            
            // write edge data
            size_t nid = 0;
            for( ; nid < get_num_nodes(); nid++)
            {
                nodes_[nid].save(oss);
                if(!oss.good()) 
                { std::cerr 
                    << "err; writing outgoing adj list for node " << nid
                    << std::endl;
                    break;
                }
            }
        }

        // print text descriptions of the set of nodes 
        // in the range [firstid, lastid)
        // NB: if the range is invalid, nothing is printed
        void
        print_dimacs_co(
                std::ostream& oss, uint32_t first_id, uint32_t last_id)
        {
            if(first_id > last_id || last_id > get_num_nodes())
            { return; }

            if(get_num_nodes() > 0)
            {
                oss << "p aux sp co " << (last_id - first_id) << std::endl;
                for(uint32_t i = first_id; i < last_id; i++)
                {
                    oss << "v " << to_external_id(i) 
                        << " " << xy_[i*2] << " " << xy_[i*2+1] << std::endl;
                }
            }
        }

        // print text descriptions of the set of arcs associated with 
        // all nodes in the range [firstid, lastid)
        // NB: if the range is invalid, nothing is printed
        void
        print_dimacs_gr(
                std::ostream& oss,
                uint32_t first_id, uint32_t last_id)
        {
            if(first_id > last_id || last_id > get_num_nodes())
            { return; }

            if(get_num_nodes() > 0)
            {
                oss << "p sp " << (last_id - first_id) << " " << 
                    this->get_num_edges_out() << std::endl;
                for(uint32_t i = first_id; i < last_id; i++)
                {
                    warthog::graph::node* n = get_node(i);
                    for(warthog::graph::edge_iter it = n->outgoing_begin(); 
                            it != n->outgoing_end(); it++)
                    {
                        oss << "a " << to_external_id(i) << " " 
                            << to_external_id((*it).node_id_)
                            << " " << (uint32_t)((*it).wt_) << std::endl;
                    }
                }
            }
        }

        // check if arc weights are Euclidean and (optionally) fix if not.
        // "fixing" means arc weights must be at least as large as the
        // Euclidean distance between the arc's head and tail vertex.
        bool
        is_euclidean(bool fix_if_not=true)
        {
            warthog::euclidean_heuristic h_euc(this);
            for(uint32_t t_id= 0; t_id < get_num_nodes(); t_id++)
            {
                int32_t tx, ty, hx, hy;
                get_xy(t_id, tx, ty);

                for(warthog::graph::edge_iter it = nodes_[t_id].outgoing_begin();
                    it != nodes_[t_id].outgoing_end();
                    it++)
                {
                    uint32_t h_id = (*it).node_id_;
                    get_xy(h_id, hx, hy);

                    warthog::cost_t hdist = h_euc.h(tx, ty, hx, hy);
                    if((*it).wt_ < hdist)
                    {
                        if(!fix_if_not) { return false; }
                        (*it).wt_ = static_cast<uint32_t>(ceil(hdist));
                    } 
                }
            }
            return true;
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
        std::vector<T_NODE> nodes_;
        
        // xy coordinates stored as adjacent pairs (x, then y)
        std::vector<int32_t> xy_;

        // these containers serve to map from external graph ids to 
        // internal graph ids
        std::vector<uint32_t> id_map_; 
        std::unordered_map<uint32_t, uint32_t> ext_id_map_; 

        std::string filename_;
        bool verbose_;
};
typedef xy_graph_base<warthog::graph::node, warthog::graph::edge> xy_graph;

}
}

#endif

