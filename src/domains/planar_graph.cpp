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
    ID_OFFSET=0;
}

warthog::graph::planar_graph::planar_graph(
        warthog::graph::planar_graph& other)
{
    verbose_ = other.verbose_;
    filename_ = other.filename_;
    nodes_sz_ = other.nodes_sz_;
    nodes_ = new warthog::graph::node[nodes_sz_];
    xy_ = new int32_t[nodes_sz_*2];
    ID_OFFSET = other.ID_OFFSET;

    for(uint32_t i = 0; i < nodes_sz_; i++)
    {
        nodes_[i] = other.nodes_[i];
        xy_[2*i] = other.xy_[2*i];
        xy_[2*i+1] = other.xy_[2*i+1];
    }
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
    dimacs.load_graph(gr_file);
    dimacs.load_graph(co_file);
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
    
    // convert nodes. NB:
    // DIMACS format specifies that node ids should be 1-indexed but
    // our internal representation is 0-indexed. We add dummy elements
    // to the front of the node arrays so we can use the 1-indexed ids
    // without any internal conversion
    ID_OFFSET = (*dimacs.nodes_begin()).id_;
    uint32_t num_nodes_dimacs = dimacs.get_num_nodes() + ID_OFFSET;
    resize(num_nodes_dimacs);
    nodes_sz_ = num_nodes_dimacs;

    for(uint32_t i = 0; i < ID_OFFSET*2; i++)
    {
        xy_[i] = INT32_MAX;
    }

    for(warthog::dimacs_parser::node_iterator it = dimacs.nodes_begin();
            it != dimacs.nodes_end(); it++)
    {
       uint32_t index = (*it).id_;
       xy_[index*2] = (*it).x_;
       xy_[index*2+1] = (*it).y_;
    }
    if(verbose_) { std::cout << "nodes, converted" << std::endl; }

    // convert edges to graph format
    warthog::euclidean_heuristic h(0);
    for(warthog::dimacs_parser::edge_iterator it = dimacs.edges_begin();
            it != dimacs.edges_end(); it++)
    {
        uint32_t hid = (*it).head_id_;
        uint32_t tid = (*it).tail_id_;

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
            double hdist = 
                h.h(xy_[tid*2], xy_[tid*2+1], 
                    xy_[hid*2], xy_[hid*2+1]);
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

    // grid costs have double precision edge weights 
    // but our graphs use integer costs; we scale
    // up all the grid costs to equivalent integers
    // (NB: 10^5 precision, then floor)
    uint32_t EDGE_COST_SCALE_FACTOR = warthog::ONE;

    warthog::gridmap gm(file);
    warthog::gridmap_expansion_policy exp(&gm);
    this->nodes_sz_ = gm.header_width() * gm.header_height();
    this->nodes_ = new warthog::graph::node[nodes_sz_];
    this->xy_ = new int32_t[nodes_sz_*2];

    for(uint32_t y = 0; y < gm.header_height(); y++)
    {
        for(uint32_t x = 0; x < gm.header_width(); x++)
        {
            uint32_t node_id = y * gm.header_width() + x;
            uint32_t gm_id = gm.to_padded_id(node_id);
            warthog::search_node* n = exp.generate(gm_id);
            warthog::search_node* nei = 0;
            double edge_cost = 0;
            exp.expand(n, 0);

            xy_[node_id*2] = x;
            xy_[node_id*2+1] = y;

            // iterate over all the neighbours of n
            for(exp.first(nei, edge_cost); nei != 0; exp.next(nei, edge_cost))
            {
                uint32_t nei_x, nei_y;
                gm.to_unpadded_xy(nei->get_id(), nei_x, nei_y);
                uint32_t nei_id = nei_y * gm.header_width() + nei_x;
                uint32_t edge_cost_int = edge_cost * EDGE_COST_SCALE_FACTOR;

                nodes_[node_id].add_outgoing(
                        warthog::graph::edge(nei_id, edge_cost_int));

                if(store_incoming)
                {
                    nodes_[nei_id].add_incoming(
                            warthog::graph::edge(node_id, edge_cost_int));
                }
            }
        }
    }
    return true;
}

void
warthog::graph::planar_graph::print_dimacs_gr(std::ostream& oss)
{
    if(nodes_sz_ > 0)
    {
        // dimacs vertex ids should always begin from index 1
        // here we apply a delta to make sure that's always so
        int32_t delta;
        if(ID_OFFSET >= 1) { delta = ID_OFFSET - 1;  }
        else { delta = -1; }

        oss << "p sp " << (nodes_sz_ - ID_OFFSET) << " " << 
            this->get_num_edges() << std::endl;
        for(uint32_t i = ID_OFFSET; i < nodes_sz_; i++)
        {
            warthog::graph::node* n = get_node(i);
            for(warthog::graph::edge_iter it = n->outgoing_begin(); 
                    it != n->outgoing_end(); it++)
            {
                oss << "a " << i - delta << " " << (*it).node_id_  - delta << " " 
                    << (uint32_t)((*it).wt_) << std::endl;
            }
        }
    }
}

void
warthog::graph::planar_graph::print_dimacs_co(std::ostream& oss)
{
    if(nodes_sz_ > 0)
    {
        // dimacs vertex ids should always begin from index 1
        // here we apply a delta to make sure that's always so
        int32_t delta;
        if(ID_OFFSET >= 1) { delta = ID_OFFSET - 1;  }
        else { delta = -1; }

        oss << "p aux sp co " << (nodes_sz_ - ID_OFFSET) << std::endl;
        for(uint32_t i = ID_OFFSET; i < nodes_sz_; i++)
        {
            oss << "v " << i - delta << " " << xy_[i*2] << " " << xy_[i*2+1] << std::endl;
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
