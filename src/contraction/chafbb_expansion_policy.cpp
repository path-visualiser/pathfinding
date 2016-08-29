#include "bb_af_filter.h"
#include "chafbb_expansion_policy.h"
#include "contraction.h"
#include "problem_instance.h"
#include "search_node.h"

warthog::chafbb_expansion_policy::chafbb_expansion_policy(
        warthog::graph::planar_graph* g, 
        std::vector<uint32_t>* rank, 
        warthog::bb_af_filter* filter,
        bool backward,
        warthog::ch::search_direction sd)
    : expansion_policy(g->get_num_nodes())
{
    g_ = g;
    rank_ = rank;
    filter_ = filter;
    backward_ = backward;
    sd_ = sd;

    assert(g->get_num_nodes() == rank->size());
}

void
warthog::chafbb_expansion_policy::expand(warthog::search_node* current,
        warthog::problem_instance* problem)
{
    reset();
    if(problem->get_searchid() != search_id_)
    {
        search_id_ = problem->get_searchid();
        filter_->set_goal(problem->get_goal());
    }

    uint32_t current_id = current->get_id();
    uint32_t current_rank = get_rank(current_id);
    warthog::graph::node* n = g_->get_node(current_id);
   
    warthog::graph::edge_iter begin, end;
    if(backward_)
    {
        begin = n->incoming_begin();
        end = n->incoming_end();
    }
    else
    {
        begin = n->outgoing_begin();
        end = n->outgoing_end();
    }

    for(warthog::graph::edge_iter it = begin; it != end; it++)
    {
        warthog::graph::edge& e = *it;
        assert(e.node_id_ < g_->get_num_nodes());

        if((sd_ & warthog::ch::UP) && get_rank(e.node_id_) > current_rank)
        {
            uint32_t from_id = current_id;
            int32_t edge_index = it - begin;
            if(backward_)
            {
                from_id = e.node_id_;
                edge_index = -1;
                warthog::graph::node* m = g_->get_node(from_id);
                for(warthog::graph::edge_iter mit = m->outgoing_begin();
                        mit != m->outgoing_end(); mit++)
                {
                    if((*mit).node_id_ == current_id)
                    {
                        edge_index = mit - m->outgoing_begin();
                        break;
                    }
                }
                assert(edge_index > -1 && edge_index < m->out_degree());
            }

            if(!filter_->filter(from_id, edge_index, backward_))
            {
                this->add_neighbour(this->generate(e.node_id_), e.wt_);
                continue;
            }
        }

        if((sd_ & warthog::ch::DOWN) && get_rank(e.node_id_) < current_rank)
        {
             this->add_neighbour(this->generate(e.node_id_), e.wt_);
        }
    }
}

uint32_t
warthog::chafbb_expansion_policy::mem()
{
    return 
        expansion_policy::mem() + 
        sizeof(this);
}
