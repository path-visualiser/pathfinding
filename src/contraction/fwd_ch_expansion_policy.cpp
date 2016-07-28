#include "contraction.h"
#include "down_distance_filter.h"
#include "euclidean_heuristic.h"
#include "fwd_ch_expansion_policy.h"
#include "planar_graph.h"
#include "search_node.h"
#include "node_filter.h"

warthog::fwd_ch_expansion_policy::fwd_ch_expansion_policy(
        warthog::graph::planar_graph* g, std::vector<uint32_t>* rank,
        warthog::down_distance_filter* nf, warthog::euclidean_heuristic* heuristic)

//        warthog::node_filter* filter)
    : expansion_policy(g->get_num_nodes()), g_(g) // , nf_(filter)
{
    rank_ = rank;
//    nf_ = filter;
    heuristic_ = 0;
    nf_ = 0;
    if(nf && heuristic)
    {
        nf_ = nf;
        heuristic_ = heuristic;
    }

    //stalled_ = new warthog::apriori_filter(g_->get_num_nodes());
}

warthog::fwd_ch_expansion_policy::~fwd_ch_expansion_policy()
{
    //delete stalled_;
}


//void
//warthog::fwd_ch_expansion_policy::expand(
//        warthog::search_node* current, warthog::problem_instance*)
//{
//    reset();
//
//    warthog::search_node* pn = current->get_parent();
//    uint32_t current_id = current->get_id();
//    uint32_t current_rank = get_rank(current_id);
//
//    warthog::graph::node* n = g_->get_node(current_id);
//    warthog::graph::edge_iter begin, end;
//    begin = n->outgoing_begin();
//    end = n->outgoing_end();
//
//    // traveling up the hierarchy we generate all neighbours;
//    // traveling down, we generate only "down" neighbours
//    bool up_travel = !pn || (current_rank > get_rank(pn->get_id()));
//    for(warthog::graph::edge_iter it = begin; it != end; it++)
//    {
//        warthog::graph::edge& e = *it;
//        assert(e.node_id_ < g_->get_num_nodes());
//        
//        if((up_travel || (get_rank(e.node_id_) < current_rank)))
//            // further prune remaining neighbours according
//            // to some strategy (usually: nodes which cannot 
//            // possibly be on an optimal path to the goal)
////            && !nf_->filter(e.node_id_))
//        {
//            this->add_neighbour(this->generate(e.node_id_), e.wt_);
//        }
//    }
//}


void
warthog::fwd_ch_expansion_policy::expand(
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
            double hval = heuristic_->h(e.node_id_, instance->get_goal());
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
warthog::fwd_ch_expansion_policy::get_xy(
        warthog::search_node* n, int32_t& x, int32_t& y)
{
    g_->get_xy(n->get_id(), x, y);
}

