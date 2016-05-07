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
}

warthog::graph::planar_graph::~planar_graph()
{
    delete [] xy_;
    delete [] nodes_;
}


bool
warthog::graph::planar_graph::load_dimacs(const char* gr_file, const char* co_file, 
        bool reverse_arcs, 
        //bool duplicate_edges, 
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
    
    // convert nodes (NB: dimacs format specifies that node ids should be
    // 1-indexed; we keep an extra dummy node to maintain compatibility)
    delete [] xy_;
    delete [] nodes_;
    nodes_sz_ = dimacs.get_num_nodes() + 1;
    xy_ = new int32_t[nodes_sz_*2];
    nodes_ = new node[nodes_sz_];

    for(warthog::dimacs_parser::node_iterator it = dimacs.nodes_begin();
            it != dimacs.nodes_end(); it++)
    {
       uint32_t index = it - dimacs.nodes_begin() + 1;
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
        if(tid == 260284)
        {
            uint32_t tmp = 32;
            tmp = tmp;
            tmp = hid;
        }
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
//        if(duplicate_edges)
//        {
//            e.node_id_ = tid;
//            nodes_[hid].add_incoming(e);
//
//            // sanity check
//            warthog::graph::edge sanity = 
//                *(nodes_[tid].outgoing_end()-1);
//            assert(sanity.node_id_ == hid);
//        }
    }
    //std::cout << "edges, converted" << std::endl;
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
        // -1 because dimacs ids are 1-indexed and we use 0-indexed arrays
        oss << "p aux sp co " << (nodes_sz_-1) << std::endl;
        oss << "p sp " << (nodes_sz_-1) << " " << 
            this->get_num_edges() << std::endl;
        for(uint32_t i = 1; i < nodes_sz_; i++)
        {
            warthog::graph::node n = nodes_[i];
            oss << "v " << i << " " << xy_[i*2] << " " << xy_[i*2+1] << std::endl;
            for(warthog::graph::edge_iter it = n.outgoing_begin(); it != n.outgoing_end(); it++)
            {
                oss << "a " << i << " " << (*it).node_id_<< " " 
                    << (*it).wt_ << std::endl;
            }
        }
    }
}
