#include "cpd.h"
#include "xy_graph.h"


void
warthog::cpd::compute_dfs_preorder(
        warthog::graph::xy_graph* g,
        std::vector<uint32_t>* column_order)
{
    assert(g->get_num_nodes() < UINT32_MAX);
    column_order->assign(g->get_num_nodes(), UINT32_MAX);

    typedef std::pair<uint32_t, uint32_t> dfs_pair;
    std::vector<dfs_pair> dfs_stack;

    uint32_t postorder_id = 0;
    uint32_t dfs_seed = ((uint32_t)rand() % (uint32_t)g->get_num_nodes());
    dfs_stack.push_back(dfs_pair(dfs_seed, 0));
    dfs_stack.reserve(g->get_num_nodes());
    while(dfs_stack.size())
    {
        // current node in the dfs tree 
        dfs_pair& dfs_node = dfs_stack.back();

        warthog::graph::node* n = g->get_node(dfs_node.first);
        if(column_order->at(dfs_node.first) == UINT32_MAX)
        { 
            // first time we reach this node. assign it an id
            column_order->at(dfs_node.first) = postorder_id++;
        }

        while(true)
        {
            warthog::graph::edge_iter begin = n->outgoing_begin() + dfs_node.second++;
            if(begin >= n->outgoing_end())
            {
                // current node has no more unreached successors. backtrack.
                dfs_stack.pop_back();
                break;
            }

            assert(begin->node_id_ < g->get_num_nodes());
            if(column_order->at(begin->node_id_) == UINT32_MAX)
            {
                // otherwise, traverse down the next branch
                dfs_stack.push_back(dfs_pair(begin->node_id_, 0));
                break;
            }
        }
    }
}

