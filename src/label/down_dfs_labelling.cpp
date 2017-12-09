#include "down_dfs_labelling.h"
#include "planar_graph.h"

#include <algorithm>

warthog::label::down_dfs_labelling::down_dfs_labelling(
        warthog::graph::planar_graph* g, 
        std::vector<uint32_t>* partitioning)
{
    dfs_order_ = new std::vector< int32_t >();
    lab_ = new std::vector< std::vector < down_dfs_label > >();

    g_ = g;
    part_ = partitioning ;

    // figure out how many bytes are required per label
    uint32_t max_id = *(std::max_element(part_->begin(), part_->end()));
    bytes_per_af_label_ = (max_id / 8) + !!(max_id % 8);
}

warthog::label::down_dfs_labelling::~down_dfs_labelling()
{ 
    delete lab_;
    delete dfs_order_;
}

void
warthog::label::down_dfs_labelling::compute_id_range_labels(
        warthog::graph::planar_graph* g, 
        std::vector<uint32_t>* rank)
{
    // allocate some memory for id-range labels
    dfs_order_->resize(g->get_num_nodes(), INT32_MAX);
    lab_->resize(g->get_num_nodes());

    // find the apex of the hierarchy
    uint32_t apex_id = 0;
    for(uint32_t i = 0; i < rank->size(); i++)
    { 
        if(rank->at(i) > rank->at(apex_id)) 
        { apex_id = i; } 
    }

    // traverse the graph and compute node and edge labels using DFS postorder
    std::vector< down_dfs_label > node_labels(
            g->get_num_nodes(), down_dfs_label(bytes_per_af_label_));

    {
        down_dfs_label dummy(bytes_per_af_label_);
        for(uint32_t n_id = 0; n_id < g->get_num_nodes(); n_id++)
        {
            uint32_t num_down = 0;
            warthog::graph::node* n = g_->get_node(n_id);
            for(warthog::graph::edge_iter it = n->outgoing_begin();
                    it != n->outgoing_end(); 
                    it++)
            { if(rank->at(it->node_id_) < rank->at(n_id)) { num_down++; } }
            lab_->at(n_id).resize(num_down, dummy);
        }
    }

    uint32_t dfs_id = 0;
    std::function<void(uint32_t)> label_fn = 
        [this, g, rank, &node_labels, &apex_id, &dfs_id, &label_fn] 
        (uint32_t source_id) -> void
        {
            down_dfs_label& s_lab = node_labels.at(source_id);
            warthog::graph::node* source = g->get_node(source_id);
            warthog::graph::edge_iter begin = source->outgoing_begin();
            uint32_t down_idx = 0;
            for( warthog::graph::edge_iter it = begin; 
                    it != source->outgoing_end();
                    it++)
            {
                // skip up edges
                if(rank->at(it->node_id_) > rank->at(source_id)) 
                { continue; }

                // DFS
                if(dfs_order_->at(it->node_id_) == INT32_MAX)
                { label_fn(it->node_id_); }

                // label the edge
                lab_->at(source_id).at(down_idx++) = 
                    node_labels.at(it->node_id_);

                // update the range down-reachable from source_id
                s_lab.merge(node_labels.at(it->node_id_));
                assert(s_lab.ids_.contains(
                            this->dfs_order_->at(it->node_id_)));
            }

            if(dfs_order_->at(source_id) == INT32_MAX)
            { this->dfs_order_->at(source_id) = dfs_id++; }

            s_lab.rank_.grow(rank->at(source_id));

            s_lab.ids_.grow(this->dfs_order_->at(source_id));

            int32_t x, y;
            g->get_xy(source_id, x, y);
            s_lab.bbox_.grow(x, y);

            uint32_t s_part = this->part_->at(source_id);
            s_lab.flags_[s_part >> 3] |= (1 << (s_part & 7)); // div8, mod8
            assert(s_lab.flags_[s_part >> 3] & (1 << (s_part & 7)));
        };

    label_fn(apex_id);
}
