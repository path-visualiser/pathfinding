#include "planar_graph.h"
#include "graph_expansion_policy.h"
#include "problem_instance.h"
#include "search_node.h"

warthog::graph_expansion_policy::graph_expansion_policy(warthog::graph::planar_graph* g)
    : expansion_policy(g->get_num_nodes())
{
    g_ = g;
    filter_sz_ = (g->get_num_nodes() >> warthog::LOG2_DBWORD_BITS)+1;
    filter_ = new warthog::dbword[filter_sz_];
}

warthog::graph_expansion_policy::~graph_expansion_policy()
{
    delete [] filter_;
}

void
warthog::graph_expansion_policy::expand(
        warthog::search_node* current, warthog::problem_instance* problem)
{
    reset();

    uint32_t current_id = current->get_id();
    warthog::graph::node* n = g_->get_node(current_id);
    
    for(warthog::graph::edge_iter it = n->outgoing_begin();
            it != n->outgoing_end(); it++)
    {
        warthog::graph::edge& e = *it;
        assert(e.node_id_ < g_->get_num_nodes());
        
        uint32_t index = e.node_id_ >> warthog::LOG2_DBWORD_BITS;
        warthog::dbword mask = 1 << (e.node_id_ & warthog::DBWORD_BITS_MASK);
        if(! (filter_[index] & mask))
        {
             this->add_neighbour(this->generate(e.node_id_), e.wt_);
        }
    }
}

void
warthog::graph_expansion_policy::get_xy(warthog::search_node* n, 
        int32_t& x, int32_t& y)
{
    g_->get_xy(n->get_id(), x, y);
}

void 
warthog::graph_expansion_policy::filter(uint32_t node_id)
{
    int index = node_id >> warthog::LOG2_DBWORD_BITS;
    int pos = node_id & DBWORD_BITS_MASK;
    filter_[index] |= 1 << pos;
}

void
warthog::graph_expansion_policy::unfilter(uint32_t node_id)
{
    int index = node_id >> warthog::LOG2_DBWORD_BITS;
    int pos = node_id & DBWORD_BITS_MASK;
    filter_[index] &= ~(1 << pos);
}

void 
warthog::graph_expansion_policy::reset_filter()
{
    for(uint32_t i = 0; i < filter_sz_; i++)
    {
        filter_[i] = 0;
    }
}

uint32_t
warthog::graph_expansion_policy::mem()
{
    return 
        expansion_policy::mem() + 
        sizeof(*filter_)*filter_sz_ +
        g_->mem() +
        sizeof(this);
}
