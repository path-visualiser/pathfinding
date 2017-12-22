#include "firstmove_labelling.h"
#include "planar_graph.h"

#include "search_node.h"
#include "zero_heuristic.h"
#include "flexible_astar.h"
#include "problem_instance.h"
#include "solution.h"
#include "graph_expansion_policy.h"

#include <algorithm>

istream&
operator>>(fm_run& the_run, std::istream& in)
{
    in.read(&the_run.head_, 4);
    in.read(&the_run.label_, 1);
    return in;
}

ostream&
operator<<(fm_run& the_run, std::ostream& out)
{
    out.write(&the_run.head_, 4);
    out.write(&the_run.label_, 1);
    return out;
}

std::istream&
warthog::label::operator>>(std::istream& in, warthog::label::dfs_label& label)
{
    in >> label.rank_ >> label.ids_ >> label.bbox_;
    for(uint32_t i = 0; i < label.flags_.size(); i++)
    {   
        in.read( (char*)(&label.flags_.at(i)), 
                  sizeof(warthog::label::dfs_label::T_FLAG) );
    }
    return in;
}

std::ostream&
warthog::label::operator<<(std::ostream& out, 
        warthog::label::dfs_label& label)
{
    out << label.rank_ << label.ids_ << label.bbox_;
    for(uint32_t i = 0; i < label.flags_.size(); i++)
    { 
        out.write( (char*)(&label.flags_.at(i)), 
                    sizeof(warthog::label::dfs_label::T_FLAG) );
    }
    return out;
}

std::istream&
warthog::label::operator>>(std::istream& in, 
        warthog::label::firstmove_labelling& lab)
{
    for(uint32_t n_id = 0; n_id < lab.g_->get_num_nodes(); n_id++)
    {
        warthog::graph::node* n = lab.g_->get_node(n_id);
        warthog::graph::edge_iter begin = n->outgoing_begin();
        warthog::graph::edge_iter end = n->outgoing_end();
        for(warthog::graph::edge_iter it = begin; it != end; it++)
        {
            in >> lab.lab_->at(n_id).at(it - begin);
            assert(in.good());

            if(!in.good())
            {
                std::cerr << "unexpected error while reading labels\n";
                std::cerr 
                    << "[debug info] node: " << n_id
                    << " out-edge-index: " << (it - begin) << "\n";
                return in;
            }
        }
    }
    return in;
}

std::ostream&
warthog::label::operator<<(std::ostream& out,
        warthog::label::firstmove_labelling& lab)
{
    for(uint32_t n_id = 0; n_id < lab.g_->get_num_nodes(); n_id++)
    {
        warthog::graph::node* n = lab.g_->get_node(n_id);
        warthog::graph::edge_iter begin = n->outgoing_begin();
        warthog::graph::edge_iter end = n->outgoing_end();
        for(warthog::graph::edge_iter it = begin; it != end; it++)
        {
            out << lab.lab_->at(n_id).at(it - begin);
            if(!out.good())
            {
                std::cerr << "unexpected error while writing labels\n";
                std::cerr 
                    << "[debug info] node: " << n_id
                    << " out-edge-index: " << (it-begin) << "\n";
                return out;
            }
        }
    }
    return out;
}

warthog::label::firstmove_labelling::firstmove_labelling(
        warthog::graph::planar_graph* g, 
        std::vector<uint32_t>* rank,
        std::vector<uint32_t>* partitioning)
    : g_(g), rank_(rank), part_(partitioning)
{
    // find the apex of the hierarchy
    apex_id_ = 0;
    for(uint32_t i = 0; i < rank_->size(); i++)
    { 
        if(rank_->at(i) > rank_->at(apex_id_)) 
        { apex_id_ = i; } 
    }

    dfs_order_ = new std::vector< int32_t >();
    compute_dfs_ids(apex_id_);
    
    // figure out how many bytes are required per label
    uint32_t max_id = *(std::max_element(part_->begin(), part_->end()));
    bytes_per_af_label_ = (max_id / 8) + !!(max_id % 8);

    // allocate memory for edge labels
    lab_ = new std::vector< std::vector < dfs_label > >();
    lab_->resize(g_->get_num_nodes());
    dfs_label dummy(bytes_per_af_label_);
    for(uint32_t n_id = 0; n_id < g_->get_num_nodes(); n_id++)
    {
        warthog::graph::node* n = this->g_->get_node(n_id);
        lab_->at(n_id).resize(n->out_degree(), dummy);
    }   
}

warthog::label::firstmove_labelling::~firstmove_labelling()
{ 
    delete lab_;
    delete dfs_order_;
}

void
warthog::label::firstmove_labelling::compute_dfs_ids(uint32_t source_id)
{
    uint32_t dfs_id = 0;
    dfs_order_->resize(g_->get_num_nodes(), INT32_MAX);
    std::function<void(uint32_t)> dfs_id_fn =
    [this, &dfs_id, &dfs_id_fn] (uint32_t source_id) -> void
    {
        warthog::graph::node* source = this->g_->get_node(source_id);
        warthog::graph::edge_iter begin = source->outgoing_begin();
        for( warthog::graph::edge_iter it = begin; 
                it != source->outgoing_end();
                it++)
        {
            // skip up edges
            if(this->rank_->at(it->node_id_) > this->rank_->at(source_id)) 
            { continue; }

            // recurse
            if(this->dfs_order_->at(it->node_id_) == INT32_MAX)
            { dfs_id_fn(it->node_id_); }
        }
        if(this->dfs_order_->at(source_id) == INT32_MAX)
        { this->dfs_order_->at(source_id) = dfs_id++; }
    };

    // gogogo
    dfs_id_fn(source_id);
}

void
warthog::label::firstmove_labelling::compress(
        std::vector<std::vector< uint32_t >>& moves, 
        std::vector<uint32_t>& column_order)
{
    for(uint32_t i = 0; i < moves.size(); i++)
    {

         
    }
}
