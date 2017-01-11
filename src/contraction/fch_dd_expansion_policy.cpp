#include "contraction.h"
#include "down_distance_filter.h"
#include "euclidean_heuristic.h"
#include "fch_dd_expansion_policy.h"
#include "planar_graph.h"
#include "search_node.h"
#include "node_filter.h"

warthog::fch_dd_expansion_policy::fch_dd_expansion_policy(
        warthog::graph::planar_graph* g, 
        std::vector<uint32_t>* rank,
        warthog::down_distance_filter* nf,
        warthog::euclidean_heuristic* heuristic)
    : expansion_policy(g->get_num_nodes()), g_(g) 
{
    rank_ = rank;
    nf_ = nf;
    heuristic_ = heuristic;
}

warthog::fch_dd_expansion_policy::~fch_dd_expansion_policy()
{
}

void
warthog::fch_dd_expansion_policy::expand(
        warthog::search_node* current, warthog::problem_instance* instance)
{
    reset();

    warthog::search_node* pn = current->get_parent();
    uint32_t current_id = current->get_id();
    uint32_t current_rank = get_rank(current_id);

    warthog::graph::node* n = g_->get_node(current_id);
    warthog::graph::edge_iter begin, end;
    begin = n->outgoing_begin();
    end = n->outgoing_end();

    // determine whether current was reached via an up edge or a down edge
    bool up_travel = !pn || (current_rank > get_rank(pn->get_id()));
    for(warthog::graph::edge_iter it = begin; it != end; it++)
    {
        warthog::graph::edge& e = *it;
        assert(e.node_id_ < g_->get_num_nodes());

        // try to prune every down successor
        if(nf_ && (get_rank(e.node_id_) < current_rank))
        {
            double ddist = nf_->get_down_distance(e.node_id_);
            double hval = heuristic_->h(e.node_id_, instance->get_target_id());
            if(ddist >= hval)
            {
                warthog::search_node* tmp = this->generate(e.node_id_);
                this->add_neighbour(tmp, e.wt_);
            }

            else
            {
// enabling below code will add neighbours that can be relaxed with
// @param current as a down-parent. the total number of nodes expanded and
// generated can be smaller with this code enabled but the time overhead
// from calling ::generate more often is usually not worth the tradeoff.
// when the code is disabled nodes with suboptimal up-parents are expanded
// and possibly some of their successors. these nodes do not affect optimality
//
//                warthog::search_node* tmp = this->generate(e.node_id_);
//                if(tmp->get_g() != warthog::INF && gval < tmp->get_g())
//                {
//                    this->add_neighbour(tmp, e.wt_);
//                }
            }
        }

        // up successors pruned only when traveling down
        else if(up_travel)
        {
            this->add_neighbour(this->generate(e.node_id_), e.wt_);
        }
    }
}

void
warthog::fch_dd_expansion_policy::get_xy(
        uint32_t id, int32_t& x, int32_t& y)
{
    g_->get_xy(id, x, y);
}

warthog::search_node* 
warthog::fch_dd_expansion_policy::generate_start_node(
        warthog::problem_instance* pi)
{
    uint32_t s_graph_id = g_->to_graph_id(pi->get_start_id());
    if(s_graph_id == warthog::INF) { return 0; }
    return generate(s_graph_id);
}

warthog::search_node*
warthog::fch_dd_expansion_policy::generate_target_node(
        warthog::problem_instance* pi)
{
    uint32_t t_graph_id = g_->to_graph_id(pi->get_target_id());
    if(t_graph_id == warthog::INF) { return 0; }
    return generate(t_graph_id);
}
