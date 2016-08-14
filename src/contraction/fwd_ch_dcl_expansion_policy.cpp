#include "contraction.h"
#include "dcl_filter.h"
#include "euclidean_heuristic.h"
#include "fwd_ch_dcl_expansion_policy.h"
#include "planar_graph.h"
#include "search_node.h"
#include "node_filter.h"

warthog::fwd_ch_dcl_expansion_policy::fwd_ch_dcl_expansion_policy(
        warthog::graph::planar_graph* g, 
        std::vector<uint32_t>* rank,
        warthog::dcl_filter* nf)
    : expansion_policy(g->get_num_nodes()), g_(g) 
{
    rank_ = rank;
    nf_ = nf;
    apex_ = 0;
    apex_reached_ = false;
}

warthog::fwd_ch_dcl_expansion_policy::~fwd_ch_dcl_expansion_policy()
{
}

void
warthog::fwd_ch_dcl_expansion_policy::expand(
        warthog::search_node* current, warthog::problem_instance* instance)
{
    reset();
    nf_->init(instance);

    warthog::search_node* pn = current->get_parent();
    uint32_t current_id = current->get_id();
    uint32_t current_rank = get_rank(current_id);

    if(rank_->at(current_id) == apex_)
    {
        apex_reached_ = true;
    }

    warthog::graph::node* n = g_->get_node(current_id);
    warthog::graph::edge_iter begin, end;
    begin = n->outgoing_begin();
    end = n->outgoing_end();

    // determine whether current was reached via an up edge or a down edge
    bool up_travel = !pn || (current_rank > get_rank(pn->get_id()));
    //std::cerr << (up_travel ? "(UPTRAVEL) " : "(DNTRAVEL) ");
    //int32_t num_down_succ = -1;
    for(warthog::graph::edge_iter it = begin; it != end; it++)
    {
        warthog::graph::edge& e = *it;
        assert(e.node_id_ < g_->get_num_nodes());

        // try to prune every down successor, regardless of 
        // wheter the parent was reached by an up edge or a
        // down edge
        bool down_succ = get_rank(e.node_id_) < current_rank;
        //num_down_succ += down_succ;
        //if(down_succ && !nf_->filter(current_id, num_down_succ, current->get_g()))
        if(down_succ && !nf_->filter(current_id, &(*it), current->get_g() + e.wt_))
        {
            // prune down successors before the apex is reached
            if(apex_ != warthog::INF && !apex_reached_) { continue; }
            // prune down successors below the goal
            if(rank_->at(e.node_id_) < rank_->at(instance->get_goal())) { continue; }

            //std::cerr << " (D) ";

            warthog::search_node* tmp = this->generate(e.node_id_);
            this->add_neighbour(tmp, e.wt_);
            continue;
        }

        // generate up successors only when traveling up
        // (the rest are implicitly pruned)
        else if(up_travel && !down_succ)
        {
            // prune up successors after the apex is reached
            if(apex_ != warthog::INF && apex_reached_) { continue; }
            // prune up successors above the apex
            if(rank_->at(e.node_id_) > apex_) { continue; }

            //std::cerr << " (U) ";

            this->add_neighbour(this->generate(e.node_id_), e.wt_);
            continue;
        }
    }
    //std::cerr << "\n";
}

void
warthog::fwd_ch_dcl_expansion_policy::get_xy(
        warthog::search_node* n, int32_t& x, int32_t& y)
{
    g_->get_xy(n->get_id(), x, y);
}

