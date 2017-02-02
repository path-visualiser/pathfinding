#include "arclabels.h"
#include "ch_expansion_policy.h"
#include "corner_point_graph.h"
#include "fch_jpg_expansion_policy.h"
#include "fch_expansion_policy.h"
#include "flexible_astar.h"
#include "graph_expansion_policy.h"
#include "planar_graph.h"
#include "problem_instance.h"
#include "zero_heuristic.h"

#include <ostream>

void
warthog::arclabels::bb_print(t_arclabels_bb& labels_, 
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
warthog::arclabels::bb_compute(
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
    arclabels::bb_print(labels_, g_, out);
}

void
warthog::arclabels::ch_bb_compute(
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
    warthog::arclabels::bb_print(labels_, g_, out);
}

void
warthog::arclabels::ch_bb_jpg_compute( 
        warthog::graph::corner_point_graph* g_, 
        std::vector<uint32_t>* rank, std::ostream& out)
{
    if(!g_ || !rank) { return; } 
    uint32_t start_id_ = 0;
    uint32_t last_id_ = g_->get_num_nodes()-1;
    t_arclabels_bb labels_;

    std::cerr << "computing bb labels\n";
    warthog::zero_heuristic heuristic;
    warthog::fch_jpg_expansion_policy expander(g_, rank);
    expander.down_successors_only(true);

    warthog::flexible_astar<
        warthog::zero_heuristic,
        warthog::fch_jpg_expansion_policy> dijkstra(&heuristic, &expander);

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
        uint32_t ext_source_id = 
            g_->get_planar_graph()->to_external_id(source_id);

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
    warthog::arclabels::bb_print(labels_, g_->get_planar_graph(), out);
}

warthog::arclabels::af_params
warthog::arclabels::get_af_params(std::vector<uint32_t>* part)
{
    warthog::arclabels::af_params retval;
    // scan all partition labels and collect them
    // into a single set (so we know how many there are)
    std::set<uint32_t> id_set;
    for(uint32_t i = 0; i < part->size(); i++)
    {
        id_set.insert(part->at(i));
    }
    retval.nparts_ = id_set.size();

    // partition ids should form a contiguous and zero-indexed set
    assert(id_set.size() > 0);
    assert(*(id_set.begin()) == 0);
    assert(*id_set.rbegin() == (retval.nparts_-1));

    retval.bytes_per_label_ = retval.nparts_ / 8;
    if((retval.nparts_ % 8) > 0)
    {
        retval.bytes_per_label_++;
    }
    return retval;
}

void
warthog::arclabels::af_print( t_arclabels_af& flags_, 
        warthog::arclabels::af_params par, std::ostream& out)
{
    if(flags_.size() == 0) { return; }

    uint32_t firstid_ = 0;
    uint32_t lastid_ = flags_.size()-1;
    out << "# to save some space, labels are written out using 64-bit words\n";
    out << "# with spaces between words for labels longer than 64 bits\n";

    out << "firstnode " << firstid_ << " lastnode " << lastid_
        << " partitions " << par.nparts_ << std::endl;

    // we split labels into 64bit words for printing
    const uint32_t word_sz = sizeof(uint64_t);
    uint32_t words_per_label = par.bytes_per_label_ / word_sz;
    if((par.bytes_per_label_ % word_sz) > 0) words_per_label++;

    // iterate over the labels for each outgoing arc of each node
    for(uint32_t i = firstid_; i <= lastid_; i++)
    {
        std::vector<uint8_t*>& node_flags = flags_.at(i);
        for(uint32_t j = 0; j < node_flags.size(); j++)
        {
            uint8_t* label = node_flags.at(j);
            for(uint32_t word = 0; word < words_per_label; word++)
            {
                uint8_t printed_word[word_sz];
                for(uint32_t k = 0; k < word_sz; k++)
                {
                    // read the label, byte by byte, one word at a time
                    if((word*word_sz+k) < par.bytes_per_label_)
                    {
                        printed_word[k] = label[word*word_sz+k];
                    }
                    // pad the last word with leading zeroes if necessary
                    else
                    {
                        printed_word[k] = 0;
                    }
                }
                out << *((uint64_t*)&(printed_word))
                    << (word < words_per_label ? " " : "");
            }
            out << std::endl;
        }
    }
}

warthog::arclabels::t_arclabels_af*
warthog::arclabels::af_compute(warthog::graph::planar_graph* g_, 
        std::vector<uint32_t>* part_, warthog::arclabels::af_params par)
{
    if(!g_) { return 0; } 

    uint32_t firstid_ = 0;
    uint32_t lastid_ = g_->get_num_nodes() - 1;

    // allocate memory for arcflag labels. 
    t_arclabels_af* flags_ = new t_arclabels_af();
    flags_->resize(g_->get_num_nodes());

    for(uint32_t i = 0; i < g_->get_num_nodes(); i++)
    {
        warthog::graph::node* n = g_->get_node(i);
        flags_->at(i).resize(n->out_degree());

        for(uint32_t j = 0; j < n->out_degree(); j++)
        {
            // allocate memory for each label and initialise 
            // each bit to zero
            uint8_t* label = new uint8_t[par.bytes_per_label_];
            flags_->at(i).at(j) = label;
            for(uint32_t k = 0; k < par.bytes_per_label_; k++)
            {
                label[k] = 0;
            }
        }
    }

    std::cerr << "computing arcflag labels\n";
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

    for(uint32_t i = firstid_; i <= lastid_; i++)
    {
        // run a dijkstra search from each node
        std::cerr << "\rprocessing node " << i << "; continues until node " 
            << lastid_ << "\r";
        uint32_t source_id = i;
        uint32_t ext_source_id = g_->to_external_id(source_id);
        warthog::problem_instance pi(ext_source_id, warthog::INF);
        dijkstra.get_length(pi);

        // now we analyse the closed list to compute arc flags
        warthog::graph::node* source = g_->get_node(source_id);

        // first, we need an easy way to convert between the ids of nodes
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

        // analyse the nodes on the closed list and label the edges of the 
        // source node accordingly
        std::function<void(warthog::search_node*)> fn_arcflags =
                [&flags_, part_, &source_id, &idmap](warthog::search_node* n)
                -> void
                {
                    // skip the source
                    assert(n);
                    if(n->get_id() == source_id) { return; } 
                    assert(n->get_parent());

                    // label the edges of the source
                    // (TODO: make this stuff faster)
                    uint32_t part_id = part_->at(n->get_id());
                    uint32_t edge_index  = (*idmap.find(
                            n->get_parent()->get_id() == source_id ? 
                            n->get_id() : n->get_parent()->get_id())).second;
                    flags_->at(source_id).at(edge_index)[part_id >> 3]
                        |= (1 << (part_id & 7));
                };
        dijkstra.apply_to_closed(fn_arcflags);
    }
    std::cerr << "\nall done\n"<< std::endl;
    return flags_;
}

warthog::arclabels::t_arclabels_af*
warthog::arclabels::ch_af_compute(warthog::graph::planar_graph* g_, 
        std::vector<uint32_t>* part_, std::vector<uint32_t>* rank, 
        warthog::arclabels::af_params par)
{
    if(!g_ || !rank) { return 0; } 

    uint32_t firstid = 0;
    uint32_t lastid = g_->get_num_nodes()-1; 

    // allocate memory for each node
    t_arclabels_af* flags_ = new t_arclabels_af();
    flags_->resize(g_->get_num_nodes());

    for(uint32_t i = 0; i < g_->get_num_nodes(); i++)
    {
        warthog::graph::node* n = g_->get_node(i);
        flags_->at(i).resize(n->out_degree());

        for(uint32_t j = 0; j < n->out_degree(); j++)
        {
            // allocate memory for each label and initialise 
            // each bit to zero
            uint8_t* label = new uint8_t[par.bytes_per_label_];
            flags_->at(i).at(j) = label;
            for(uint32_t k = 0; k < par.bytes_per_label_; k++)
            {
                label[k] = 0;
            }
        }
    }

    std::cerr << "computing ch arcflag labels\n";
    warthog::zero_heuristic heuristic;
    warthog::fch_expansion_policy expander(g_, rank);

    warthog::flexible_astar<
        warthog::zero_heuristic,
        warthog::fch_expansion_policy> dijkstra(&heuristic, &expander);

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

    for(uint32_t i = firstid; i <= lastid; i++)
    {
        // run a dijkstra search from each node
        std::cerr << "\rprocessing node " << i << "; continues until node " 
            << (lastid) << "\r";
        uint32_t source_id = i;
        uint32_t ext_source_id = g_->to_external_id(source_id);
        warthog::problem_instance pi(ext_source_id, warthog::INF);
        dijkstra.get_length(pi);

        // now we analyse the closed list to compute arc flags
        warthog::graph::node* source = g_->get_node(source_id);

        // first, we need an easy way to convert between the ids of nodes
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

        // analyse the nodes on the closed list and label the edges of the 
        // source node accordingly
        std::function<void(warthog::search_node*)> fn_arcflags =
                [&flags_, part_, &source_id, &idmap](warthog::search_node* n) -> void
                {
                    // skip the source
                    assert(n);
                    if(n->get_id() == source_id) { return; } 
                    assert(n->get_parent());

                    // label the edges of the source
                    // (TODO: make this stuff faster)
                    uint32_t part_id = part_->at(n->get_id());
                    uint32_t edge_index  = (*idmap.find(
                            n->get_parent()->get_id() == source_id ? 
                            n->get_id() : n->get_parent()->get_id())).second;
                    flags_->at(source_id).at(edge_index)[part_id >> 3]
                        |= (1 << (part_id & 7));
                };
        dijkstra.apply_to_closed(fn_arcflags);
    }
    std::cerr << "\nall done\n"<< std::endl;
    return flags_;
}

warthog::arclabels::t_arclabels_af*
warthog::arclabels::ch_af_jpg_compute(warthog::graph::corner_point_graph* g, 
        std::vector<uint32_t>* part, std::vector<uint32_t>* rank,
        warthog::arclabels::af_params par)
{
    return 0;
}
