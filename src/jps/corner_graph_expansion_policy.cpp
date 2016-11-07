#include "corner_graph.h"
#include "corner_graph_expansion_policy.h"
#include "corner_point_locator.h"
#include "online_jump_point_locator2.h"
#include "planar_graph.h"

warthog::jps::corner_graph_expansion_policy::corner_graph_expansion_policy(
        warthog::graph::corner_graph* g) 
    : expansion_policy(g->get_num_nodes())
{
    g_ = g;
}

warthog::jps::corner_graph_expansion_policy::~corner_graph_expansion_policy()
{
}

void
warthog::jps::corner_graph_expansion_policy::get_xy(
        warthog::search_node* n, int32_t& x, int32_t& y)
{
    assert(n);
    g_->get_xy(n->get_id(), x, y);
}

void 
warthog::jps::corner_graph_expansion_policy::expand(
        warthog::search_node* current, warthog::problem_instance* problem)
{
    reset();

    uint32_t graph_id = current->get_id();
    warthog::graph::node* n = g_->get_node(graph_id);
    assert(n);
    for(warthog::graph::edge_iter eit = n->outgoing_begin(); 
            eit != n->outgoing_end(); eit++)
    {
        add_neighbour(generate(eit->node_id_), eit->wt_);
    }
}

size_t
warthog::jps::corner_graph_expansion_policy::mem()
{
    return 
        expansion_policy::mem() + sizeof(*this) + g_->mem();
}

