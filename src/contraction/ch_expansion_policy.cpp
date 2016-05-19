#include "ch_expansion_policy.h"
#include "problem_instance.h"
#include "search_node.h"

void
warthog::ch_expansion_policy::expand(warthog::search_node* current,
        warthog::problem_instance* problem)
{
    reset();

    uint32_t current_id = current->get_id();
    warthog::graph::node* n = g_->get_node(current_id);
    
    for(warthog::graph::edge_iter it = n->outgoing_begin();
            it != n->outgoing_end(); it++)
    {
        warthog::graph::edge& e = *it;
        assert(e.node_id_ < g_->get_num_nodes());
        
        if(order_->at(e.node_id_) > order_->at(current_id))
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
