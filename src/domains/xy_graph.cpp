#include "xy_graph.h"

void
warthog::graph::gridmap_to_xy_graph(
    warthog::gridmap* gm, warthog::graph::xy_graph* g,
    bool store_incoming)
{
    assert(gm && g && g->get_num_nodes() == 0);
    g->set_filename(gm->filename());
    g->grow(gm->get_num_traversable_tiles());

    // add each traversable tile as a node in the graph
    std::vector<uint32_t> id_map(gm->header_height() * gm->header_width());
    warthog::gridmap_expansion_policy exp(gm);
    uint32_t next_graph_id = 0;
    for(uint32_t y = 0; y < gm->header_height(); y++)
    {
        for(uint32_t x = 0; x < gm->header_width(); x++)
        {
            // we differentiate between external grid ids
            // (unpadded) and internal grid ids (with padding)
            uint32_t from_gm_id = y * gm->header_width() + x;

            // skip obstacles
            if(!gm->get_label(gm->to_padded_id(from_gm_id))) 
            { continue; }

            // add graph node (we scale up all costs and coordinates)
            id_map[from_gm_id] = next_graph_id;
            g->set_xy(next_graph_id, 
                (int32_t)(x * warthog::graph::GRID_TO_GRAPH_SCALE_FACTOR),
                (int32_t)(y * warthog::graph::GRID_TO_GRAPH_SCALE_FACTOR));
            next_graph_id++;
        }
    }

    // add edges
    for(uint32_t y = 0; y < gm->header_height(); y++)
    {
        for(uint32_t x = 0; x < gm->header_width(); x++)
        {
            // we use a grid-expansion policy here which enforces
            // standard constraints on grid moves re; e.g. corner cutting
            uint32_t from_graph_id;
            {
                from_graph_id = id_map[y*gm->header_width() + x];
                assert(from_graph_id < next_graph_id);
            }

            warthog::search_node* nei = 0;
            warthog::search_node* n = 0;
            edge_cost_t edge_cost = 0;

            n = exp.generate(gm->to_padded_id(y*gm->header_width() + x));
            exp.expand(n, 0);
            for(exp.first(nei, edge_cost); nei != 0; exp.next(nei, edge_cost))
            {
                uint32_t to_graph_id;
                {
                    uint32_t nei_x, nei_y;
                    gm->to_unpadded_xy((uint32_t)nei->get_id(), nei_x, nei_y);
                    uint32_t to_gm_id = nei_y * gm->header_width() + nei_x;
                    to_graph_id = id_map[to_gm_id];
                    assert(from_graph_id != to_graph_id);
                    assert(to_graph_id < next_graph_id);
                }

                warthog::graph::node* gr_from = g->get_node(from_graph_id);
                edge_cost_t gr_weight = edge_cost * warthog::graph::GRID_TO_GRAPH_SCALE_FACTOR;
                gr_from->add_outgoing(warthog::graph::edge(to_graph_id, gr_weight));
                if(store_incoming)
                {
                    warthog::graph::node* gr_to = g->get_node(to_graph_id);
                    gr_to->add_incoming(warthog::graph::edge(from_graph_id, gr_weight));
                }
            }
        }
    }
}
