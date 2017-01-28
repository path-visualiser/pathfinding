#include "arclabels_bb.h"
#include "ch_expansion_policy.h"
#include "corner_point_graph.h"
#include "flexible_astar.h"
#include "graph_expansion_policy.h"
#include "planar_graph.h"
#include "problem_instance.h"
#include "zero_heuristic.h"

#include <ostream>

void
warthog::arclabels_bb_print(t_arclabels_bb& labels_, 
        warthog::graph::planar_graph* g_, std::ostream& out)
{
    out << "# Bounding box edge labels for the graph file \n"
        << "# " << g_->get_filename() << "\n" 
        << "# Before printing labels are sorted by tail index. \n"
        << "# Ties are broken using the original input ordering \n"
        << "# (lower is better)\n";
    out << "firstnode " << 0 << " lastnode " << (labels_.size()-1) << "\n";
    warthog::geom::rectangle dummy;
    for(uint32_t i = 0; i < labels_.size(); i++)
    {
        std::vector<warthog::geom::rectangle> vrect = labels_.at(i);
        for(uint32_t j = 0; j < vrect.size(); j++)
        {
            warthog::geom::rectangle rect = vrect.at(j);
            if(rect != dummy) { assert(rect.is_valid()); }
            out 
                //<< i << " " << j << " " << 
                << rect.x1 << "\t" << rect.y1 << "\t"
                << rect.x2 << "\t" << rect.y2 << "\t"
                << std::endl;
        }
    }
}

void
warthog::arclabels_bb_compute(
        warthog::graph::planar_graph* g_, std::ostream& out)
{
    if(!g_) { return; } 

    uint32_t start_id_ = 0;
    uint32_t last_id_ = g_->get_num_nodes()-1;
    t_arclabels_bb labels_;

    std::cerr << "computing bb labels\n";
    warthog::zero_heuristic heuristic;
    warthog::graph_expansion_policy expander(g_);

    warthog::flexible_astar<
        warthog::zero_heuristic,
        warthog::graph_expansion_policy> dijkstra(&heuristic, &expander);

    // need to keep track of the first edge on the way to the current node
    // (the solution is a bit hacky as we break the chain of backpointers 
    // to achieve this; it doesn't matter, we don't care about the path)
    std::function<void(warthog::search_node*)> relax_fn = 
            [] (warthog::search_node* n) -> void
            {
                // the start node and its children don't need their 
                // parent pointers updated. for all other nodes the
                // grandparent becomes the parent
                if(n->get_parent()->get_parent() != 0)
                {
                    if(n->get_parent()->get_parent()->get_parent() != 0)
                    {
                        n->set_parent(n->get_parent()->get_parent());
                    }
                }
            };
    dijkstra.apply_on_relax(relax_fn);

    for(uint32_t i = start_id_; i <= last_id_; i++)
    {
        std::cerr << "\rprocessing node : " << i 
                  << ". continue to " << last_id_;

        uint32_t source_id = i;
        uint32_t ext_source_id = g_->to_external_id(source_id);
        warthog::problem_instance pi(ext_source_id, warthog::INF);
        dijkstra.get_length(pi);
        
        // create a rectangle label for every outgoing edge
        warthog::graph::node* source = g_->get_node(source_id);
        labels_.push_back(
                std::vector<warthog::geom::rectangle>(
                    source->out_degree()));

        // we need an easy way to convert between the ids of nodes
        // adjacent to the source and their corresponding edge index
        std::unordered_map<uint32_t, uint32_t> idmap;
        uint32_t edge_idx = 0;
        for(warthog::graph::edge_iter it = source->outgoing_begin(); 
                it != source->outgoing_end(); it++)
        {
            idmap.insert(
                    std::pair<uint32_t, uint32_t>((*it).node_id_, edge_idx));
            edge_idx++;
        }

        // compute the extent of each rectangle bounding box
        std::function<void(warthog::search_node*)> bbox_fn = 
                [g_, &labels_, &source_id, &idmap](warthog::search_node* n) -> void
                {
                    // skip the source (it doesn't belong to any rectangle)
                    assert(n);
                    if(n->get_id() == source_id) { return; } 
                    assert(n->get_parent());

                    // grow the rectangle associated with each outgoing
                    // edge of the source node
                    std::unordered_map<uint32_t, uint32_t>::iterator it;
                    if(n->get_parent()->get_parent() == 0)
                    // special case for the successors of the source
                    { it = idmap.find(n->get_id()); }
                    else
                    // all the other nodes
                    { it = idmap.find(n->get_parent()->get_id()); }

                    if(it == idmap.end())
                    {
                        std::cerr << "err; invalid first edge; "
                             << " n id: "<< n->get_id() << " parent id: " 
                             << n->get_parent()->get_id() << std::endl;
                        exit(0);
                    }
                    int32_t x, y;
                    g_->get_xy(n->get_id(), x, y);
                    assert(x != warthog::INF && y != warthog::INF);
                    labels_.back().at((*it).second).grow(x, y);
                    assert(labels_.back().at((*it).second).is_valid());
                };
        dijkstra.apply_to_closed(bbox_fn);
    }
    std::cerr << "\nall done\n"<< std::endl;
    arclabels_bb_print(labels_, g_, out);
}

void
warthog::arclabels_ch_bb_compute(
        warthog::graph::planar_graph* g_, 
        std::vector<uint32_t>* rank, std::ostream& out)
{
    if(!g_ || !rank) { return; } 
    uint32_t start_id_ = 0;
    uint32_t last_id_ = g_->get_num_nodes()-1;
    t_arclabels_bb labels_;

    std::cerr << "computing bb labels\n";
    warthog::zero_heuristic heuristic;
    warthog::ch_expansion_policy expander(g_, rank, false, warthog::ch::DOWN);

    warthog::flexible_astar<
        warthog::zero_heuristic,
        warthog::ch_expansion_policy> dijkstra(&heuristic, &expander);

    // need to keep track of the first edge on the way to the current node
    // (the solution is a bit hacky as we break the chain of backpointers 
    // to achieve this; it doesn't matter, we don't care about the path)
    std::function<void(warthog::search_node*)> relax_fn = 
            [] (warthog::search_node* n) -> void
            {
                // the start node and its children don't need their 
                // parent pointers updated. for all other nodes the
                // grandparent becomes the parent
                if(n->get_parent()->get_parent() != 0)
                {
                    if(n->get_parent()->get_parent()->get_parent() != 0)
                    {
                        n->set_parent(n->get_parent()->get_parent());
                    }
                }
            };
    dijkstra.apply_on_relax(relax_fn);

    for(uint32_t i = start_id_; i <= last_id_; i++)
    {
        std::cerr 
            << "\rprocessing node: " << i << ". continue to " << last_id_;
        uint32_t source_id = i;
        uint32_t ext_source_id = g_->to_external_id(source_id);
        dijkstra.get_length( 
                warthog::problem_instance(ext_source_id, warthog::INF));
        
        // create a rectangle label for every outgoing edge
        warthog::graph::node* source = g_->get_node(source_id);
        labels_.push_back(
                std::vector<warthog::geom::rectangle>(
                    source->out_degree()));

        // we need an easy way to convert between the ids of nodes
        // adjacent to the source and their corresponding edge index
        std::unordered_map<uint32_t, uint32_t> idmap;
        for(warthog::graph::edge_iter it = source->outgoing_begin(); 
                it != source->outgoing_end(); it++)
        {
            uint32_t edge_idx = it - source->outgoing_begin();
            idmap.insert(
                    std::pair<uint32_t, uint32_t>((*it).node_id_, edge_idx));
        }
        assert(idmap.size() == source->out_degree());

        // compute the extent of each rectangle bounding box
        std::function<void(warthog::search_node*)> bbox_fn = 
            [g_, &labels_, &source_id, &idmap](warthog::search_node* n) 
            -> void
            {
                // skip the source (it doesn't belong to any rectangle)
                assert(n);
                if(n->get_id() == source_id) { return; } 
                assert(n->get_parent());

                // grow the rectangle associated with each outgoing
                // edge of the source node
                std::unordered_map<uint32_t, uint32_t>::iterator it;
                if(n->get_parent()->get_parent() == 0)
                // special case for the successors of the source
                { it = idmap.find(n->get_id()); }
                else
                // all the other nodes
                { it = idmap.find(n->get_parent()->get_id()); }



                if(it == idmap.end())
                {
                    std::cerr << "err; invalid first edge; "
                         << " n id: "<< n->get_id() << " parent id: " 
                         << n->get_parent()->get_id() << std::endl;
                    exit(0);
                }
                int32_t x, y;
                g_->get_xy(n->get_id(), x, y);
                assert(x != warthog::INF && y != warthog::INF);
                labels_.back().at((*it).second).grow(x, y);
                assert(labels_.back().at((*it).second).is_valid());
            };

        dijkstra.apply_to_closed(bbox_fn);
    }
    std::cerr << "\nall done\n"<< std::endl;
    warthog::arclabels_bb_print(labels_, g_, out);
}

void
warthog::arclabels_ch_bb_jpg_compute( 
        warthog::graph::corner_point_graph* g, 
        std::vector<uint32_t>* rank, std::ostream& out)
{
}
