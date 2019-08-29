#include "contraction.h"
#include "fch_expansion_policy.h"
#include "xy_graph.h"
#include "search_node.h"

warthog::fch_expansion_policy::fch_expansion_policy(
        warthog::graph::xy_graph* g, 
        std::vector<uint32_t>* rank)
    : expansion_policy(g->get_num_nodes()), g_(g) 
{
    rank_ = rank;
}

void
warthog::fch_expansion_policy::expand(
        warthog::search_node* current, warthog::problem_instance*)
{
    reset();

    warthog::search_node* pn = generate(current->get_parent());
    warthog::sn_id_t current_id = current->get_id();
    uint32_t current_rank = get_rank(current_id);
    warthog::graph::node* n = g_->get_node((uint32_t)current_id);


    // traveling up the hierarchy we generate all neighbours;
    // traveling down, we generate only "down" neighbours
    bool up_travel = !pn || (current_rank > get_rank(pn->get_id()));
    if(up_travel)
    {
        // generate up edges
        for(uint32_t i = 0; i < n->out_degree(); i++)
        {
            warthog::graph::edge& e = *(n->outgoing_begin() + i);
            assert(e.node_id_ < g_->get_num_nodes());
            this->add_neighbour(this->generate(e.node_id_), e.wt_);
        }
        // generate down edges
        // NB: we abuse the incoming list and store here the 
        // set of outgoing successors which go down the hierarchy
        for(uint32_t i = 0; i < n->in_degree(); i++)
        {
            warthog::graph::edge& e = *(n->incoming_begin() + i);
            assert(e.node_id_ < g_->get_num_nodes());
            this->add_neighbour(this->generate(e.node_id_), e.wt_);
        }
    }
    else
    {
        for(uint32_t i = 0; i < n->in_degree(); i++)
        {
            warthog::graph::edge& e = *(n->incoming_begin() + i);
            assert(e.node_id_ < g_->get_num_nodes());
            this->add_neighbour(this->generate(e.node_id_), e.wt_);
        }
    }
}

void
warthog::fch_expansion_policy::get_xy(warthog::sn_id_t nid, int32_t& x, int32_t& y)
{
    g_->get_xy((uint32_t)nid, x, y);
}

warthog::search_node* 
warthog::fch_expansion_policy::generate_start_node(
        warthog::problem_instance* pi)
{
    uint32_t s_graph_id = g_->to_graph_id((uint32_t)pi->start_id_);
    if(s_graph_id == warthog::INF32) { return 0; }
    return generate(s_graph_id);
}

warthog::search_node*
warthog::fch_expansion_policy::generate_target_node(
        warthog::problem_instance* pi)
{
    uint32_t t_graph_id = g_->to_graph_id((uint32_t)pi->target_id_);
    if(t_graph_id == warthog::INF32) { return 0; }
    return generate(t_graph_id);
}
