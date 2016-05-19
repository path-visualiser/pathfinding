#include "constants.h"
#include "dimacs_parser.h"
#include "euclidean_heuristic.h"
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
        bool duplicate_edges, 
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
    
    delete [] xy_;
    delete [] nodes_;

    // convert nodes. NB:
    // DIMACS format specifies that node ids should be 1-indexed but
    // our internal representation is 0-indexed. We add dummy elements
    // to the front of the node arrays so we can use the 1-indexed ids
    // without any internal conversion
    ID_OFFSET = 1;
    nodes_sz_ = dimacs.get_num_nodes() + ID_OFFSET;
    xy_ = new int32_t[nodes_sz_*2];
    nodes_ = new node[nodes_sz_];
    xy_[0] = xy_[1] = INT32_MAX;

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
        nodes_[tid].add_outgoing(e);

        // duplicate edges are stored twice: once as an incoming edge 
        // and once as an outgoing edge
        if(duplicate_edges)
        {
            e.node_id_ = tid;
            nodes_[hid].add_incoming(e);

            // sanity check
#ifndef NDEBUG
            warthog::graph::edge sanity = 
                *(nodes_[tid].outgoing_end()-1);
            assert(sanity.node_id_ == hid);
#endif
        }
    }

    if(verbose_) { std::cout << "edges, converted" << std::endl; }
    return true;
}

bool
warthog::graph::planar_graph::load_grid(const char* file)
{
    std::cerr << "not yet!";
    return 1;
}

void
warthog::graph::planar_graph::print_dimacs(std::ostream& oss)
{
    if(nodes_sz_ > 0)
    {
        // we add +1 to all ids because dimacs ids are 1-indexed and our
        // internal representation is 0-indexed
        oss << "p aux sp co " << (nodes_sz_) << std::endl;
        oss << "p sp " << (nodes_sz_) << " " << 
            this->get_num_edges() << std::endl;
        for(uint32_t i = ID_OFFSET; i < nodes_sz_; i++)
        {
            warthog::graph::node n = nodes_[i];
            oss << "v " << i << " " << xy_[i*2] << " " << xy_[i*2+1] << std::endl;
            for(warthog::graph::edge_iter it = n.outgoing_begin(); it != n.outgoing_end(); it++)
            {
                oss << "a " << i << " " << (*it).node_id_ << " " 
                    << (*it).wt_ << std::endl;
            }
        }
    }
}

void
warthog::graph::planar_graph::print_dimacs_gr(std::ostream& oss)
{
    if(nodes_sz_ > 0)
    {
        // -1 because dimacs ids are 1-indexed and we use 0-indexed arrays
        oss << "p sp " << (nodes_sz_ - ID_OFFSET) << " " << 
            this->get_num_edges() << std::endl;
        for(uint32_t i = ID_OFFSET; i < nodes_sz_; i++)
        {
            warthog::graph::node n = nodes_[i];
            for(warthog::graph::edge_iter it = n.outgoing_begin(); it != n.outgoing_end(); it++)
            {
                oss << "a " << i << " " << (*it).node_id_ << " " 
                    << (*it).wt_ << std::endl;
            }
        }
    }
}

void
warthog::graph::planar_graph::print_dimacs_co(std::ostream& oss)
{
    if(nodes_sz_ > 0)
    {
        // -1 because dimacs ids are 1-indexed and we use 0-indexed arrays
        oss << "p aux sp co " << (nodes_sz_ - ID_OFFSET) << std::endl;
        for(uint32_t i = ID_OFFSET; i < nodes_sz_; i++)
        {
            warthog::graph::node n = nodes_[i];
            oss << "v " << i << " " << xy_[i*2] << " " << xy_[i*2+1] << std::endl;
        }
    }
}
