#include "contraction.h"
#include "fch_expansion_policy.h"
#include "planar_graph.h"
#include "search_node.h"

warthog::fch_expansion_policy::fch_expansion_policy(
        warthog::graph::planar_graph* g, std::vector<uint32_t>* rank)
    : expansion_policy(g->get_num_nodes()), g_(g) 
{
    rank_ = rank;
}

warthog::fch_expansion_policy::~fch_expansion_policy()
{
}

void
warthog::fch_expansion_policy::expand(
        warthog::search_node* current, warthog::problem_instance*)
{
    reset();

    warthog::search_node* pn = current->get_parent();
    uint32_t current_id = current->get_id();
    uint32_t current_rank = get_rank(current_id);

    warthog::graph::node* n = g_->get_node(current_id);
    warthog::graph::edge_iter begin, end;
    begin = n->outgoing_begin();
    end = n->outgoing_end();

    // traveling up the hierarchy we generate all neighbours;
    // traveling down, we generate only "down" neighbours
    bool up_travel = !pn || (current_rank > get_rank(pn->get_id()));
    for(warthog::graph::edge_iter it = begin; it != end; it++)
    {
        warthog::graph::edge& e = *it;
        assert(e.node_id_ < g_->get_num_nodes());
        
        if(up_travel || (get_rank(e.node_id_) < current_rank))
        {
            this->add_neighbour(this->generate(e.node_id_), e.wt_);
        }
    }
}

void
warthog::fch_expansion_policy::get_xy(uint32_t nid, int32_t& x, int32_t& y)
{
    g_->get_xy(nid, x, y);
}

warthog::search_node* 
warthog::fch_expansion_policy::generate_start_node(
        warthog::problem_instance* pi)
{
    uint32_t s_graph_id = g_->to_graph_id(pi->start_id_);
    if(s_graph_id == warthog::INF) { return 0; }
    return generate(s_graph_id);
}

warthog::search_node*
warthog::fch_expansion_policy::generate_target_node(
        warthog::problem_instance* pi)
{
    uint32_t t_graph_id = g_->to_graph_id(pi->target_id_);
    if(t_graph_id == warthog::INF) { return 0; }
    return generate(t_graph_id);
}
