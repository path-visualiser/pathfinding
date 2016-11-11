#include "jpg_expansion_policy.h"
#include "corner_graph.h"

warthog::jps::jpg_expansion_policy::jpg_expansion_policy(
        warthog::graph::corner_graph* g)
    : expansion_policy(g->get_num_nodes())
{
    g_ = g;
}

warthog::jps::jpg_expansion_policy::~jpg_expansion_policy()
{
}

void
warthog::jps::jpg_expansion_policy::get_xy(
        warthog::search_node* n, int32_t& x, int32_t& y)
{
    assert(n);
    return g_->get_xy(n->get_id(), x, y);
}

void 
warthog::jps::jpg_expansion_policy::expand(
        warthog::search_node* current, warthog::problem_instance* problem)
{
    reset();
    uint32_t graph_id = current->get_id();

    // special case for the start node (which has no parent)
    if(graph_id == problem->get_start())
    {
        warthog::graph::node* n = g_->get_node(graph_id);
        for(warthog::graph::edge_iter eit = n->outgoing_begin(); 
                eit != n->outgoing_end(); eit++)
        {
            add_neighbour(generate(eit->node_id_), eit->wt_);
        }
    }
    // all other nodes, all of which are in the graph
    else
    {
        // compute the travel direction from the parent
        int32_t px, py, x, y;
        this->get_xy(current->get_parent(), px, py);
        this->get_xy(current, x, y);
        warthog::jps::direction dir_c = 
            warthog::jps::compute_direction((uint32_t)px, (uint32_t)py, 
                    (uint32_t)x, (uint32_t)y);
        
        // compute the set of directions in which we should 
        // generate successors
        warthog::gridmap* gm = g_->get_gridmap();
        uint32_t grid_id = gm->to_padded_id(x, y);
        uint32_t tiles;
        gm->get_neighbours(grid_id, (uint8_t*)&tiles);
        uint32_t succ_dirs = warthog::jps::compute_successors(dir_c, tiles);

        warthog::graph::node* n = g_->get_node(graph_id);
        warthog::graph::edge_iter eit = n->outgoing_begin();
        while(succ_dirs)
        {
            uint32_t label = __builtin_ffs(succ_dirs)-1;
            uint32_t first = g_->labelled_edge_offset(graph_id, label);
            uint32_t last = g_->labelled_edge_offset(graph_id, label+1);

            for( uint32_t i = first; i < last; i++)
            {
                add_neighbour(generate((eit+i)->node_id_), (eit+i)->wt_);
            }
            succ_dirs = succ_dirs & ~(1 << label);
        }

        // goal test
        eit += n->out_degree() - 1;
        if(eit->node_id_ == g_->get_inserted_target_id())
        {
            add_neighbour(generate(eit->node_id_), eit->wt_);
        }
    }
}

size_t
warthog::jps::jpg_expansion_policy::mem()
{
    return 
        expansion_policy::mem() + sizeof(*this) + g_->mem();
}

