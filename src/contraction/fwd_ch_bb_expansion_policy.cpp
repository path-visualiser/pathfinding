#include "contraction.h"
#include "bbox_filter.h"
#include "euclidean_heuristic.h"
#include "fwd_ch_bb_expansion_policy.h"
#include "planar_graph.h"
#include "search_node.h"
#include "node_filter.h"

warthog::fwd_ch_bb_expansion_policy::fwd_ch_bb_expansion_policy(
        warthog::graph::planar_graph* g, 
        std::vector<uint32_t>* rank,
        warthog::bbox_filter* nf)
    : expansion_policy(g->get_num_nodes()), g_(g) 
{
    rank_ = rank;
    nf_ = nf;
}

warthog::fwd_ch_bb_expansion_policy::~fwd_ch_bb_expansion_policy()
{
}

void
warthog::fwd_ch_bb_expansion_policy::expand(
        warthog::search_node* current, warthog::problem_instance* instance)
{
    reset();
    nf_->init(instance);

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

        // try to prune every down successor, regardless of 
        // wheter the parent was reached by an up edge or a
        // down edge
        bool down_succ = get_rank(e.node_id_) < current_rank;
        if(down_succ && !nf_->filter__(current_id, (it - begin)))
        {
            warthog::search_node* tmp = this->generate(e.node_id_);
            this->add_neighbour(tmp, e.wt_);
            continue;
        }

        // generate up successors only when traveling up
        // (the rest are implicitly pruned)
        else if(up_travel && !down_succ)
        {
            this->add_neighbour(this->generate(e.node_id_), e.wt_);
            continue;
        }
    }
}

void
warthog::fwd_ch_bb_expansion_policy::get_xy(
        warthog::search_node* n, int32_t& x, int32_t& y)
{
    g_->get_xy(n->get_id(), x, y);
}

