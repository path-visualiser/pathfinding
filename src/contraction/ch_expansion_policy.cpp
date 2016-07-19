#include "ch_expansion_policy.h"
#include "contraction.h"
#include "problem_instance.h"
#include "search_node.h"

warthog::ch_expansion_policy::ch_expansion_policy(
        warthog::graph::planar_graph* g, 
        std::vector<uint32_t>* ooc, bool backward) 
    : expansion_policy(g->get_num_nodes()), g_(g)
{
    backward_ = backward;

    assert(g->get_num_nodes() - g->id_offset() == ooc->size());

    // we convert the order of contraction into a ranked 
    // list s.t. rank[i] gives the contraction index of
    // node i
    warthog::ch::convert_order_of_contraction_to_ranked_list(
            *ooc, rank_);
}


void
warthog::ch_expansion_policy::expand(warthog::search_node* current,
        warthog::problem_instance* problem)
{
    reset();

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
        
        if(get_rank(e.node_id_) > current_rank)
        {
             this->add_neighbour(this->generate(e.node_id_), e.wt_);
        }
    }
}

uint32_t
warthog::ch_expansion_policy::mem()
{
    return 
        expansion_policy::mem() + 
        sizeof(this);
}
