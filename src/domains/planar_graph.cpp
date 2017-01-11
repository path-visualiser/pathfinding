#include "constants.h"
#include "dimacs_parser.h"
#include "euclidean_heuristic.h"
#include "gridmap.h"
#include "gridmap_expansion_policy.h"
#include "planar_graph.h"

#include <algorithm>
#include <cassert>
#include <functional>
#include <iostream>
#include <string>

warthog::graph::planar_graph::planar_graph()
{
    xy_ = 0;
    nodes_ = 0;
    nodes_sz_ = 0;
    nodes_cap_ = 0;
    verbose_ = false;
}

warthog::graph::planar_graph::planar_graph(
        warthog::gridmap* map, bool store_incoming)
{
    xy_ = 0;
    nodes_ = 0;
    nodes_sz_ = 0;
    nodes_cap_ = 0;
    verbose_ = false;

    grid2graph(map, store_incoming);
}

warthog::graph::planar_graph::planar_graph(
        warthog::graph::planar_graph& other)
{
    (*this) = other;
}

warthog::graph::planar_graph&
warthog::graph::planar_graph::operator=(const warthog::graph::planar_graph& other)
{
    verbose_ = other.verbose_;
    filename_ = other.filename_;
    resize(other.nodes_sz_); // NB: minimum memory for new graph

    for(uint32_t i = 0; i < nodes_sz_; i++)
    {
        nodes_[i] = other.nodes_[i];
        xy_[2*i] = other.xy_[2*i];
        xy_[2*i+1] = other.xy_[2*i+1];
    }

    id_map_.clear();
    for(uint32_t i = 0; i < other.id_map_.size(); i++)
    {
        id_map_.push_back(other.id_map_.at(i));
    }
    return *this;
}

warthog::graph::planar_graph::~planar_graph()
{
    delete [] xy_;
    delete [] nodes_;
}


bool
warthog::graph::planar_graph::load_dimacs(const char* gr_file, const char* co_file, 
        bool reverse_arcs, 
        bool store_incoming_edges, 
        bool enforce_euclidean)
{
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
    
    uint32_t num_nodes_dimacs = dimacs.get_num_nodes();
    resize(num_nodes_dimacs);

    for(warthog::dimacs_parser::node_iterator it = dimacs.nodes_begin();
            it != dimacs.nodes_end(); it++)
    {
       add_node((*it).x_, (*it).y_, (*it).id_);
    }
    if(verbose_) { std::cerr << "nodes, converted" << std::endl; }

    // convert edges to graph format
    warthog::euclidean_heuristic h(0);
    for(warthog::dimacs_parser::edge_iterator it = dimacs.edges_begin();
            it != dimacs.edges_end(); it++)
    {
        uint32_t hid = to_graph_id( (*it).head_id_ );
        uint32_t tid = to_graph_id( (*it).tail_id_ );

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

        // in a euclidean consistent graph the edge costs are at least
        // as large as the euclidean norm between the edge endpoints
        if(enforce_euclidean)
        {
            int32_t tx, ty, hx, hy;
            get_xy(tid, tx, ty);
            get_xy(hid, hx, hy);

            double hdist = h.h(tx, ty, hx, hy);
            if(e.wt_ < hdist)
            {
                e.wt_ = ceil(hdist);
            } 
        }


        // add edge to the graph
#ifndef NDEBUG
        uint32_t deg_before = nodes_[tid].out_degree();
#endif

        nodes_[tid].add_outgoing(e);

#ifndef NDEBUG
        uint32_t deg_after = nodes_[tid].out_degree();
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
    }

    if(verbose_) { std::cout << "edges, converted" << std::endl; }
    return true;
}

// creates a graph from a gridmap input file. edge transition
// costs are fixed to 1 for straight moves and sqrt(2) for diagonal.
bool
warthog::graph::planar_graph::load_grid(
        const char* file, bool store_incoming)
{
    delete [] nodes_;
    delete [] xy_;

    warthog::gridmap gm(file);
    return grid2graph(&gm, store_incoming);
}

bool
warthog::graph::planar_graph::grid2graph(warthog::gridmap* gm, bool store_incoming)
{
    // enumerate the traversable nodes with degree > 0
    struct edge_tuple 
    {
        uint32_t from;
        uint32_t to;
        uint32_t wt;
    };

    filename_ = gm->filename();
    warthog::gridmap_expansion_policy exp(gm);
    std::vector<uint32_t> nodes;
    std::vector<edge_tuple> edges;
    for(uint32_t y = 0; y < gm->header_height(); y++)
    {
        for(uint32_t x = 0; x < gm->header_width(); x++)
        {
            uint32_t node_id = y * gm->header_width() + x;
            uint32_t gm_id = gm->to_padded_id(node_id);

            // skip obstacles
            if(!gm->get_label(gm_id)) { continue; }

            // each traversable node appears in the graph
            add_node(x * warthog::graph::GRID_TO_GRAPH_SCALE_FACTOR, 
                     y * warthog::graph::GRID_TO_GRAPH_SCALE_FACTOR,
                     node_id);

            warthog::search_node* n = exp.generate(gm_id);
            warthog::search_node* nei = 0;
            double edge_cost = 0;
            exp.expand(n, 0);
            for(exp.first(nei, edge_cost); nei != 0; exp.next(nei, edge_cost))
            {
                //// each traversable node with degree > 0 
                //// appears in the graph
                //if(id_map_.size() == 0 || id_map_.back() != node_id) 
                //{ 
                //    add_node(x * warthog::graph::GRID_TO_GRAPH_SCALE_FACTOR, 
                //             y * warthog::graph::GRID_TO_GRAPH_SCALE_FACTOR,
                //             node_id);
                //}

                // track the list of traversable edges
                uint32_t nei_x, nei_y;
                gm->to_unpadded_xy(nei->get_id(), nei_x, nei_y);
                uint32_t nei_id = nei_y * gm->header_width() + nei_x;

                assert(edge_cost > 0);

                edge_tuple e;
                e.from = node_id, 
                e.to = nei_id, 
                e.wt = edge_cost * warthog::graph::GRID_TO_GRAPH_SCALE_FACTOR;
                edges.push_back(e);
            }
        }
    }

    for(uint32_t i = 0; i < edges.size(); i++)
    {
        edge_tuple tuple = edges.at(i);    
        uint32_t from_id = to_graph_id(tuple.from);
        uint32_t to_id = to_graph_id(tuple.to);
        assert(from_id != warthog::INF && to_id != INF);

        warthog::graph::node* from = this->get_node(from_id);
        warthog::graph::node* to = this->get_node(to_id);
        assert(from && to);

        warthog::graph::edge e(to_id, tuple.wt);
        from->add_outgoing(warthog::graph::edge(to_id, tuple.wt));
        if(store_incoming)
        {
            to->add_incoming(warthog::graph::edge(from_id, tuple.wt));
        }
    }
    return true;
}

void
warthog::graph::planar_graph::print_dimacs_gr(std::ostream& oss,
                uint32_t first_id, uint32_t last_id)
{
    if(first_id > last_id || last_id > nodes_sz_)
    { return; }

    if(nodes_sz_ > 0)
    {
        oss << "p sp " << nodes_sz_ << " " << 
            this->get_num_edges() << std::endl;
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

void
warthog::graph::planar_graph::print_dimacs_co(std::ostream& oss,
                uint32_t first_id, uint32_t last_id)
{
    if(first_id > last_id || last_id > nodes_sz_)
    { return; }

    if(nodes_sz_ > 0)
    {
        // dimacs vertex ids should always begin from index 1
        // here we apply a delta to make sure that's always so
        oss << "p aux sp co " << nodes_sz_ << std::endl;
        for(uint32_t i = first_id; i < last_id; i++)
        {
            oss << "v " << to_external_id(i) 
                << " " << xy_[i*2] << " " << xy_[i*2+1] << std::endl;
        }
    }
}

size_t
warthog::graph::planar_graph::resize(uint32_t new_cap)
{
    if(new_cap <= nodes_cap_) { return false; }
    warthog::graph::node* big_nodes = new warthog::graph::node[new_cap];
    int32_t* big_xy = new int32_t[new_cap*2];
    for(uint32_t i = 0; i < nodes_sz_; i++)
    {
        big_nodes[i] = std::move(nodes_[i]);
        big_xy[i*2] = xy_[i*2];
        big_xy[i*2+1] = xy_[i*2+1];
    }
    delete [] nodes_;
    delete [] xy_;

    nodes_ = big_nodes;
    xy_ = big_xy;
    nodes_cap_ = new_cap;

    return nodes_cap_;
}


uint32_t
warthog::graph::planar_graph::add_node(int32_t x, int32_t y, uint32_t ext_id)
{
    if(nodes_cap_ == nodes_sz_)
    {
        resize(nodes_cap_ == 0 ? 1 : (nodes_cap_*2));
    }
    uint32_t index = nodes_sz_;
    nodes_sz_++;
    xy_[index*2] = x;
    xy_[index*2+1] = y;
    
    if(ext_id != warthog::INF)
    {
        id_map_.push_back(ext_id);
    }
    return index;
}
