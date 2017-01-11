#include "planar_graph.h"
#include "graph_expansion_policy.h"
//#include "node_filter.h"
#include "problem_instance.h"
#include "search_node.h"

warthog::graph_expansion_policy::graph_expansion_policy(warthog::graph::planar_graph* g)
    : expansion_policy(g->get_num_nodes())
{
    g_ = g;
    nf_ = 0;
}

warthog::graph_expansion_policy::~graph_expansion_policy()
{
}

void
warthog::graph_expansion_policy::expand(
        warthog::search_node* current, warthog::problem_instance* problem)
{
    reset();

    uint32_t current_id = current->get_id();
    warthog::graph::node* n = g_->get_node(current_id) ;
    
    for(warthog::graph::edge_iter it = n->outgoing_begin();
            it != n->outgoing_end(); it++)
    {
        warthog::graph::edge& e = *it;
        assert(e.node_id_ < g_->get_num_nodes());
        
//        if(!nf_ || !(nf_->filter(e.node_id_)))
//        {
             this->add_neighbour(this->generate(e.node_id_), e.wt_);
//        }
    }
}

void
warthog::graph_expansion_policy::get_xy(uint32_t nid, int32_t& x, int32_t& y)
{
    g_->get_xy(nid, x, y);
}


size_t
warthog::graph_expansion_policy::mem()
{
    return 
        expansion_policy::mem() + 
        //nf_->mem() +
        g_->mem() +
        sizeof(this);
}

warthog::search_node* 
warthog::graph_expansion_policy::generate_start_node(
        warthog::problem_instance* pi)
{
    uint32_t s_graph_id = g_->to_graph_id(pi->get_start_id());
    if(s_graph_id == warthog::INF) { return 0; }
    return generate(s_graph_id);
}

warthog::search_node*
warthog::graph_expansion_policy::generate_target_node(
        warthog::problem_instance* pi)
{
    uint32_t t_graph_id = g_->to_graph_id(pi->get_target_id());
    if(t_graph_id == warthog::INF) { return 0; }
    return generate(t_graph_id);
}
