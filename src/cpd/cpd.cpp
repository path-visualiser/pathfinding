#include "cpd.h"
#include "xy_graph.h"


void
warthog::cpd::compute_dfs_preorder(
        warthog::graph::xy_graph* g, 
        std::vector<uint32_t>* column_order)
{
    std::cerr << "warthog::cpd::compute_dfs_preorder\n";
    typedef std::pair<uint32_t, uint32_t> dfs_pair;

    // dfs bookkeeping (we store some data to avoid cycles)
    std::vector<uint32_t> tmp(g->get_num_nodes(), INT32_MAX); 
    std::function<void(uint32_t)> dfs_preorder =
        [column_order, &tmp, g](uint32_t dfs_seed) -> void
    {
        std::vector<dfs_pair> dfs_stack;
        dfs_stack.push_back(dfs_pair(dfs_seed, 0));
        dfs_stack.reserve(g->get_num_nodes()*2);
        while(dfs_stack.size())
        {
            // current node in the dfs tree 
            dfs_pair dfs_node = dfs_stack.back();
            dfs_stack.pop_back();
            warthog::graph::node* n = g->get_node(dfs_node.first);

            //std::cerr << dfs_node.first << " " << dfs_node.second;

            // pre-order traversal
            assert(dfs_node.first != INT32_MAX);
            if(tmp.at(dfs_node.first) == INT32_MAX)
            { 
                tmp.at(dfs_node.first) = dfs_node.first;
                column_order->push_back(dfs_node.first); 
            }

            warthog::graph::edge_iter begin = n->outgoing_begin();
            warthog::graph::edge_iter left = begin + dfs_node.second;

            // right branch 
            for(uint32_t i = dfs_node.second+1; i < n->out_degree(); i++)
            {
                assert((begin+i)->node_id_ != INT32_MAX);
                if(tmp.at((begin+i)->node_id_) == INT32_MAX)
                {
                    //std::cerr << " right " << (begin+i)->node_id_;
                    dfs_stack.push_back(dfs_pair(dfs_node.first, i));
                    break;
                }
            }

            // left branch 
            if(left == n->outgoing_end()) { return; }
            assert(left->node_id_ != INT32_MAX);
            if(tmp.at(left->node_id_) == INT32_MAX)
            {
                //std::cerr << " left " << left->node_id_;
                dfs_stack.push_back(dfs_pair(left->node_id_, 0)); 
            }

            //std::cerr << "\n";
        }
    };

    // pick a random seed and run DFS; continue until every node is labelled
    while(true)
    {
        std::vector<uint32_t> seeds;
        for(uint32_t i = 0; i < tmp.size(); i++)
        {
            if(tmp.at(i) == INT32_MAX) { seeds.push_back(i); } 
        }
        if(!seeds.size()) { break; }

        uint32_t dfs_seed = seeds.at((uint32_t)rand() % (uint32_t)seeds.size());
        dfs_preorder(dfs_seed);
    }
}
