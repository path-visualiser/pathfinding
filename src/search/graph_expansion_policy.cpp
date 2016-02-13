#include "planar_graph.h"
#include "graph_expansion_policy.h"
#include "problem_instance.h"
#include "search_node.h"

warthog::graph_expansion_policy::graph_expansion_policy(warthog::planar_graph* g)
    : expansion_policy(g->get_num_nodes())
{
    this->g_ = g;
}

warthog::graph_expansion_policy::~graph_expansion_policy()
{
    delete g_;
}

void
warthog::graph_expansion_policy::expand(
        warthog::search_node* current, warthog::problem_instance* problem)
{
    reset();

    uint32_t current_id = current->get_id();
    warthog::planar_graph::node n = g_->get_node(current_id);
    
    for(uint32_t i = n.begin_; i < n.begin_ + n.degree_; i++)
    {
        warthog::planar_graph::edge e = g_->get_edge(i);
        assert(e.head_idx_ < g_->get_num_nodes());
        this->add_neighbour(this->generate(e.head_idx_), e.wt_);
    }
}

void
warthog::graph_expansion_policy::get_xy(warthog::search_node* n, 
        int32_t& x, int32_t& y)
{
    g_->get_xy(n->get_id(), x, y);
}
