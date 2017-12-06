#include "down_labelling.h"
#include "planar_graph.h"

#include <algorithm>

warthog::label::down_labelling::down_labelling()
{
    dfs_order_ = new std::vector<int32_t>();
    lab_idr_ = new std::vector<std::vector< id_range_label >>();
}

warthog::label::down_labelling::~down_labelling()
{
    delete lab_bbaf_;
    delete lab_idr_;
    delete dfs_order_;
}

void
warthog::label::down_labelling::compute_id_range_labels(
        warthog::graph::planar_graph* g, 
        std::vector<uint32_t>* rank)
{
    // allocate some memory for id-range labels
    dfs_order_->resize(g->get_num_nodes(), INT32_MAX);
    lab_idr_->resize(g->get_num_nodes());

    // find the apex of the hierarchy
    uint32_t apex_id = 0;
    for(uint32_t i = 0; i < rank->size(); i++)
    { 
        if(rank->at(i) > rank->at(apex_id)) 
        { apex_id = i; } 
    }

    // traverse the graph and compute node and edge labels using DFS postorder
    std::vector<id_range_label> node_range(g->get_num_nodes());
    uint32_t dfs_id = 0;

    std::function<id_range_label(uint32_t)> label_fn = 
        [this, g, rank, &node_range, &apex_id, &dfs_id, &label_fn] 
        (uint32_t source_id) -> id_range_label
        {
            id_range_label source_range;

            warthog::graph::node* source = g->get_node(source_id);
            warthog::graph::edge_iter begin = source->outgoing_begin();
            for( warthog::graph::edge_iter it = begin; 
                    it != source->outgoing_end();
                    it++)
            {
                // skip up edges
                if(rank->at(it->node_id_) > rank->at(source_id)) 
                { continue; }

                // DFS
                id_range_label subtree_range = node_range.at((*it).node_id_);
                if(subtree_range.left == INT32_MAX)
                {
                    subtree_range = label_fn(it->node_id_);
                    node_range.at(it->node_id_) = subtree_range;
                }

                // label the edge
                lab_idr_->at(source_id).push_back(subtree_range);

                // update the range down-reachable from source_id
                source_range.merge(subtree_range);

                assert(source_range.contains(
                            this->dfs_order_->at(it->node_id_)));
            }
            if(this->dfs_order_->at(source_id) == INT32_MAX)
            {
               this->dfs_order_->at(source_id) = dfs_id++;
            }
            source_range.grow(this->dfs_order_->at(source_id));
            return source_range;
        };

    label_fn(apex_id);
}
