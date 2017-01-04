#include "af_filter.h"
#include "constants.h"
#include "flexible_astar.h"
#include "fch_expansion_policy.h"
#include "graph_expansion_policy.h"
#include "planar_graph.h"
#include "problem_instance.h"
#include "zero_heuristic.h"

#include <fstream>
#include <iostream>
#include <set>

warthog::af_filter::af_filter(
        warthog::graph::planar_graph* g,
        std::vector<uint32_t>* part)
{
    init(g, part);
}

warthog::af_filter::af_filter(
        warthog::graph::planar_graph* g,
        std::vector<uint32_t>* part,
        const char* arcflags_file)
{
    init(g, part);
    load_labels(arcflags_file);
}

warthog::af_filter::~af_filter()
{
    for(uint32_t i = 0; i < flags_.size(); i++)
    {
        std::vector<uint8_t*>& node_flags = flags_.at(i);
        for(uint32_t j = 0; j < node_flags.size(); j++)
        {
            uint8_t* label = node_flags.back();
            delete [] label;
            node_flags.pop_back();
        }
    }
}

void 
warthog::af_filter::init(
        warthog::graph::planar_graph* g, 
        std::vector<uint32_t>* part)
{
    g_ = g;
    part_ = part;

    firstid_ = lastid_ = 0;

    // scan all partition labels and collect them
    // into a single set (so we know how many there are)
    std::set<uint32_t> id_set;
    for(uint32_t i = 0; i < part->size(); i++)
    {
        id_set.insert(part->at(i));
    }

    // partition ids should form a contiguous and zero-indexed set
    assert(id_set.size() > 0);
    nparts_ = id_set.size();
    assert(*(id_set.begin()) == 0);
    assert(*id_set.rbegin() == (nparts_-1));

    bytes_per_label_ = nparts_ / 8;
    if((nparts_ % 8) > 0)
    {
        bytes_per_label_++;
    }

    // allocate memory for arcflag labels. 
    for(uint32_t i = 0; i < g_->get_num_nodes(); i++)
    {
        flags_.push_back(std::vector<uint8_t*>());
        warthog::graph::node* n = g_->get_node(i);
        for(uint32_t j = 0; j < n->out_degree(); j++)
        {
            uint8_t* label = new uint8_t[bytes_per_label_];
            flags_.back().push_back(label);

            // each bit of each label is initialised to zero
            for(uint32_t k = 0; k < bytes_per_label_; k++)
            {
                label[k] = 0;
            }
        }
    }
    t_byte_ = 0;
    t_bitmask_ = 0;
}

void
warthog::af_filter::print(std::ostream& out)
{
    out << "# to save some space, labels are written out using 64-bit words\n";
    out << "# with spaces between words for labels longer than 64 bits\n";

    out << "firstnode " << firstid_ << " lastnode " << lastid_
        << " partitions " << nparts_ << std::endl;

    // we split labels into 64bit words for printing
    const uint32_t word_sz = sizeof(uint64_t);
    uint32_t words_per_label = bytes_per_label_ / word_sz;
    if((bytes_per_label_ % word_sz) > 0) words_per_label++;

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
                    if((word*word_sz+k) < bytes_per_label_)
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

bool
warthog::af_filter::load_labels(const char* filename)
{
    std::cerr << "loading arcflags file\n";
    std::ifstream ifs(filename);

    uint32_t num_parts;

    // skip comment lines
    while(ifs.peek() == '#')
    {
        while(ifs.get() != '\n');
    }

    // read header
    bool bad_header = false;
    std::string token;

    ifs >> token;
    if(token == "firstnode")
    {
        ifs >> firstid_;

        ifs >> token;
        if(token == "lastnode")
        {
            ifs >> lastid_;
        }
        else { bad_header = true; }

        ifs >> token;
        if(token == "partitions")
        {
            ifs >> num_parts;
            if(num_parts != nparts_)
            { 
                std::cerr 
                    << "err; number of partitions in arcflags file\n"
                    << " does not match input graph (read: " << num_parts 
                    << " expected: " << nparts_ << ")\n";
                return false;
            }
        }
        else { bad_header = true; }
    }
    else { bad_header = true; }

    if(bad_header)
    {
        std::cerr << "err; arcflags file has invalid header. required format: \n"
            << "firstnode [integer] lastnode [integer] partitions [integer]\n";
        return false;
    }

    // read labels for each outgoing arc
    const uint32_t word_sz = sizeof(uint64_t);
    uint32_t words_per_label = ceil(bytes_per_label_ / (double)word_sz);
    for(uint32_t i = firstid_; i < (lastid_+1); i++)
    {
        warthog::graph::node* n = g_->get_node(i);
        for(uint32_t j = 0; j < n->out_degree(); j++)
        {
            for(uint32_t word = 0; word < words_per_label; word++)
            {
                if(!ifs.good())
                {
                    std::cerr << "errr; unexpected read error loading arcflags; "
                        << "aborting\n";
                    return false;
                }

                uint64_t label = 0;
                ifs >> label;
                //std::cout << "node " << i << " arc " << j << 
                //" label " << label << std::endl;
                *((uint64_t*)&flags_.at(i).at(j)[word*word_sz]) = label;
                assert(*((uint64_t*)&flags_.at(i).at(j)[word*word_sz])
                            == label);
            }
        }
    }
    return true;
}

void
warthog::af_filter::compute_ch(std::vector<uint32_t>* rank)
{

    compute_ch(0, g_->get_num_nodes()-1, rank);

}

// compute arclabels for all nodes whose ids are in
// the range [firstid, lastid]
void
warthog::af_filter::compute_ch(
        uint32_t firstid, uint32_t lastid, std::vector<uint32_t>* rank)
{
    if(!g_ || !rank) { return; } 

    firstid_ = firstid;
    lastid_ = lastid;
    if(lastid_ >= g_->get_num_nodes()) { lastid_ = g_->get_num_nodes() - 1; }

    std::cerr << "computing ch arcflag labels " 
        << "for all nodes in the id-range [" 
        << firstid_ << ", " << lastid_ << "]\n";
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
        warthog::problem_instance pi(source_id, warthog::INF);
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
                [this, &source_id, &idmap](warthog::search_node* n) -> void
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
                    flags_.at(source_id).at(edge_index)[part_id >> 3]
                        |= (1 << (part_id & 7));
                };
        dijkstra.apply_to_closed(fn_arcflags);
    }
    std::cerr << "\nall done\n"<< std::endl;
}

void
warthog::af_filter::compute()
{

    compute(0, g_->get_num_nodes()-1);

}

// compute arclabels for all nodes whose ids are in
// the range [firstid, lastid]
void
warthog::af_filter::compute(uint32_t firstid, uint32_t lastid)
{
    if(!g_) { return; } 

    firstid_ = firstid;
    lastid_ = lastid;
    if(lastid_ >= g_->get_num_nodes()) { lastid_ = g_->get_num_nodes() - 1; }

    std::cerr << "computing ch arcflag labels " 
        << "for all nodes in the id-range [" 
        << firstid_ << ", " << lastid_ << "]\n";
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

    for(uint32_t i = 0; i < (lastid_-firstid_)+1; i++)
    {
        // run a dijkstra search from each node
        std::cerr << "\rprocessing node " << i << "; continues until node " 
            << (lastid_-firstid_) << "\r";
        uint32_t source_id = i + firstid_;
        warthog::problem_instance pi(source_id, warthog::INF);
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
                [this, &source_id, &idmap](warthog::search_node* n) -> void
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
                    flags_.at(source_id).at(edge_index)[part_id >> 3]
                        |= (1 << (part_id & 7));
                };
        dijkstra.apply_to_closed(fn_arcflags);
    }
    std::cerr << "\nall done\n"<< std::endl;
}
