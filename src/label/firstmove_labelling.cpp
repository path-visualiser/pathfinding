#include "firstmove_labelling.h"
#include "planar_graph.h"

#include <algorithm>

warthog::label::firstmove_labelling::firstmove_labelling(
        warthog::graph::planar_graph* g)
{
    lab_ = new std::vector< std::vector < warthog::label::fm_run > >();
    lab_->resize(g_->get_num_nodes());
}

warthog::label::firstmove_labelling::~firstmove_labelling()
{ 
    delete lab_;
}

istream&
operator>>(fm_run& the_run, std::istream& in)
{
    in.read(&the_run.head_, 4);
    in.read(&the_run.label_, 1);
    return ininline ;
}

ostream&
operator<<(fm_run& the_run, std::ostream& out)
{
    out.write(&the_run.head_, 4);
    out.write(&the_run.label_, 1);
    return out;
}

std::istream&
warthog::label::operator>>(std::istream& in, 
        warthog::label::firstmove_labelling& lab)
{
    // grr.. label size should be multiples of 8, not 64
    fm_run dummy{UINT32_MAX, UI
    for(uint32_t n_id = 0; n_id < lab.g_->get_num_nodes(); n_id++)
    {
        warthog::graph::node* n = lab.g_->get_node(n_id);
        warthog::graph::edge_iter begin = n->outgoing_begin();
        warthog::graph::edge_iter end = n->outgoing_end();
        for(warthog::graph::edge_iter it = begin; it != end; it++)
        {
            lab.push_back(std::vector< fm_run >(0));

            // some rows might not be populated

            in >> fm_run;
            if(fm_run.head_ = UINT32_MAX) { continue; }
                
            lab.lab_->at(n_id).at(it - begin);
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

