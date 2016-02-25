#include "constants.h"
#include "dimacs_parser.h"
#include "euclidean_heuristic.h"
#include "planar_graph.h"

#include <algorithm>
#include <cassert>
#include <functional>
#include <iostream>
#include <string>

warthog::planar_graph::planar_graph()
{
    filename_ = new std::string("");
    nodes_ = new std::vector<warthog::planar_graph::node>();
    edges_ = new std::vector<warthog::planar_graph::edge>();
}

warthog::planar_graph::~planar_graph()
{
    nodes_->clear();
    edges_->clear();
    delete nodes_;
    delete edges_;
    delete filename_;
}

bool
warthog::planar_graph::load_dimacs(const char* gr_file, const char* co_file, 
        bool backward_graph, bool enforce_euclidean)
{
    //warthog::dimacs_parser dimacs(gr_file, co_file);
    warthog::dimacs_parser dimacs;
    dimacs.load_graph(gr_file);
    dimacs.load_graph(co_file);
    //std::cout << "loaded" << std::endl;

    nodes_->reserve(dimacs.get_num_nodes());
    edges_->reserve(dimacs.get_num_edges());


    // sort the edges by their tail indexes
    // (or head indexes, for a backward graph)
    std::sort(dimacs.edges_begin(), dimacs.edges_end(), 
            [=](warthog::dimacs_parser::edge e1,
               warthog::dimacs_parser::edge e2) -> bool 
                {
                    if(backward_graph)
                    {
                        if(e1.head_id_ < e2.head_id_)
                        {
                            return true;
                        }
                        if(e1.head_id_ == e2.head_id_)
                        {
                            return e1.head_id_ <= e2.head_id_;
                        }
                        return false;
                    }

                    if(e1.tail_id_ < e2.tail_id_)
                    {
                        return true;
                    }
                    if(e1.tail_id_ == e2.tail_id_)
                    {
                        return e1.head_id_ <= e2.head_id_;
                    }
                    return false;
                }
    );
    //std::cout << "edges, sorted" << std::endl;
  
    // sort nodes by their ids
    bool (*node_comparator) 
         (warthog::dimacs_parser::node, warthog::dimacs_parser::node) =
            [](warthog::dimacs_parser::node n1,
               warthog::dimacs_parser::node n2) -> bool
                {
                    return n1.id_ < n2.id_;
                };
    std::sort(dimacs.nodes_begin(), dimacs.nodes_end(), node_comparator);
    //std::cout << "nodes, sorted" << std::endl;

    // convert nodes to graph format
    // first, insert a dummy element at index0. that way all the id
    // of each node is equal to its index in nodes_
    warthog::planar_graph::node n;
    n.x_ = warthog::INF;
    n.y_ = warthog::INF;
    n.begin_ = n.degree_ = 0;
    nodes_->push_back(n);

    for(warthog::dimacs_parser::node_iterator it = dimacs.nodes_begin();
            it != dimacs.nodes_end(); it++)
    {
       warthog::planar_graph::node n;
       n.degree_ = 0;
       n.begin_ = 0;
       n.x_ = (*it).x_;
       n.y_ = (*it).y_;
       assert((*it).id_ == nodes_->size());
       nodes_->push_back(n);
    }
    //std::cout << "nodes, converted" << std::endl;

    warthog::euclidean_heuristic h(0);

    // convert edges to graph format
    for(warthog::dimacs_parser::edge_iterator it = dimacs.edges_begin();
            it != dimacs.edges_end(); it++)
    {
        //std::cout << "debug; tail_id = " << (*it).tail_id_ << 
        //    " head_id = " << (*it).head_id_ << std::endl;
        warthog::planar_graph::edge e;
        e.head_idx_ = backward_graph ? 
            ((*it).tail_id_) :
            ((*it).head_id_);
        e.wt_ = (*it).weight_;

        // in a euclidean consistent graph the edge costs are at least
        // as large as the euclidean norm between the edge endpoints
        if(enforce_euclidean)
        {
            warthog::planar_graph::node head = nodes_->at((*it).head_id_);
            warthog::planar_graph::node tail = nodes_->at((*it).tail_id_);
            double hdist = h.h(tail.x_, tail.y_, head.x_, head.y_);
            if(e.wt_ < hdist)
            {
                e.wt_ = ceil(hdist);
            } 
        }
        edges_->push_back(e);

        // while we're at it, update node degree info
        // (NB: nodes_ is 0-indexed; DIMACS format is 1-indexed)
        int node_id = backward_graph ? 
                        ((*it).head_id_) :
                        ((*it).tail_id_);
        if((*nodes_)[node_id].degree_ == 0)
        {
            (*nodes_)[node_id].begin_ = it - dimacs.edges_begin();
        }
        (*nodes_)[node_id].degree_++;
    }
    //std::cout << "edges, converted" << std::endl;
    return true;
}

bool
warthog::planar_graph::load_grid(const char* file)
{
    std::cerr << "not yet!";
    return 1;
}


void
warthog::planar_graph::print_dimacs(std::ostream& oss)
{
    uint32_t nnodes = nodes_->size();
    if(nnodes > 0)
    {
        oss << "p aux sp co " << nodes_->size() << std::endl;
        oss << "p sp " << nodes_->size() << " " << edges_->size() 
            << std::endl;
        for(uint32_t i = 0; i < nnodes; i++)
        {
            warthog::planar_graph::node n = nodes_->at(i);
            oss << "v " << (i+1) << " " << n.x_ << " " << n.y_ << std::endl;
            for(uint32_t j = n.begin_; j < n.begin_ + n.degree_; j++)
            {
                warthog::planar_graph::edge e = edges_->at(j);
                oss << "a " << (i+1) << " " << e.head_idx_ << " " << e.wt_
                    << std::endl;
            }
        }
    }

}


