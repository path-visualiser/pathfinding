#include "dimacs_parser.h"
#include "graph.h"

#include <algorithm>
#include <functional>
#include <iostream>
#include <string>

warthog::graph::graph()
{
    filename_ = new std::string("");
    nodes_ = new std::vector<warthog::graph::node>();
    edges_ = new std::vector<warthog::graph::edge>();
    lat_lng_ = new std::vector<uint32_t>();
}

warthog::graph::~graph()
{
    delete nodes_;
    delete edges_;
    delete lat_lng_;
    delete filename_;
}

bool
warthog::graph::load_dimacs(char* gr_file, char* co_file)
{
    //warthog::dimacs_parser dimacs(gr_file, co_file);
    warthog::dimacs_parser dimacs;
    dimacs.load(gr_file);
    dimacs.load(co_file);
    std::cout << "loaded" << std::endl;

    nodes_->reserve(dimacs.get_num_nodes());
    edges_->reserve(dimacs.get_num_edges());


    // sort the edges by their tail indexes
    bool (*edge_comparator)
         (warthog::dimacs_parser::edge, warthog::dimacs_parser::edge) =
            [](warthog::dimacs_parser::edge e1,
               warthog::dimacs_parser::edge e2) -> bool 
                {
                    if(e1.tail_id_ < e2.tail_id_)
                    {
                        return true;
                    }
                    if(e1.tail_id_ == e2.tail_id_)
                    {
                        return e1.head_id_ <= e2.head_id_;
                    }
                    return false;
                };
    std::sort(dimacs.edges_begin(), dimacs.edges_end(), edge_comparator); 
    std::cout << "edges, sorted" << std::endl;
  
    // sort nodes by their ids
    bool (*node_comparator) 
         (warthog::dimacs_parser::node, warthog::dimacs_parser::node) =
            [](warthog::dimacs_parser::node n1,
               warthog::dimacs_parser::node n2) -> bool
                {
                    return n1.id_ < n2.id_;
                };
    std::sort(dimacs.nodes_begin(), dimacs.nodes_end(), node_comparator);
    std::cout << "nodes, sorted" << std::endl;

    // convert nodes to graph format
    for(warthog::dimacs_parser::node_iterator it = dimacs.nodes_begin();
            it != dimacs.nodes_end(); it++)
    {
       warthog::graph::node n;
       n.degree_ = 0;
       n.begin_ = 0;
       n.x_ = (*it).x_;
       n.y_ = (*it).y_;
       nodes_->push_back(n);
    }
    std::cout << "nodes, converted" << std::endl;

    // convert edges to graph format
    for(warthog::dimacs_parser::edge_iterator it = dimacs.edges_begin();
            it != dimacs.edges_end(); it++)
    {
        //std::cout << "debug; tail_id = " << (*it).tail_id_ << 
        //    " head_id = " << (*it).head_id_ << std::endl;
        warthog::graph::edge e;
        e.head_idx_ =  (*it).head_id_;
        e.wt_ = (*it).weight_;
        edges_->push_back(e);

        // while we're at it, update node degree info
        int node_id = (*it).tail_id_ - 1; // because nodes_ is 0-indexed
        if((*nodes_)[node_id].degree_ == 0)
        {
            (*nodes_)[node_id].begin_ = it - dimacs.edges_begin();
        }
        (*nodes_)[node_id].degree_++;
    }
    std::cout << "edges, converted" << std::endl;
    return true;
}

bool
warthog::graph::load_grid(char* file)
{
    std::cerr << "not yet!";
    return 1;
}


void
warthog::graph::print_dimacs(std::ostream& oss)
{
    uint32_t nnodes = nodes_->size();
    if(nnodes > 0)
    {
        oss << "p aux sp co " << nodes_->size() << std::endl;
        oss << "p sp " << nodes_->size() << " " << edges_->size() 
            << std::endl;
        for(uint32_t i = 0; i < nnodes; i++)
        {
            warthog::graph::node n = nodes_->at(i);
            oss << "v " << (i+1) << " " << n.x_ << " " << n.y_ << std::endl;
            for(uint32_t j = n.begin_; j < n.begin_ + n.degree_; j++)
            {
                warthog::graph::edge e = edges_->at(j);
                oss << "a " << (i+1) << " " << e.head_idx_ << " " << e.wt_
                    << std::endl;
            }
        }
    }

}


