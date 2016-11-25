#include "bb_af_filter.h"
#include "constants.h"
#include "flexible_astar.h"
#include "fwd_ch_expansion_policy.h"
#include "planar_graph.h"
#include "problem_instance.h"
#include "zero_heuristic.h"

#include <fstream>
#include <iostream>
#include <set>

warthog::bb_af_filter::bb_af_filter(
        warthog::graph::planar_graph* g,
        std::vector<uint32_t>* rank,
        std::vector<uint32_t>* part)
{
    init(g, rank, part);
}

warthog::bb_af_filter::bb_af_filter(
        warthog::graph::planar_graph* g,
        std::vector<uint32_t>* rank,
        std::vector<uint32_t>* part,
        const char* bb_af_file)
{
    init(g, rank, part);
    load_labels(bb_af_file);
}

warthog::bb_af_filter::~bb_af_filter()
{
    for(uint32_t i = 0; i < labels_.size(); i++)
    {
        std::vector<bbaf_label>& arclabs = labels_.at(i);
        for(uint32_t j = 0; j < arclabs.size(); j++)
        {
            bbaf_label& label = arclabs.back();
            delete [] label.flags_;
            arclabs.pop_back();
        }
    }
}

void 
warthog::bb_af_filter::init(
        warthog::graph::planar_graph* g, 
        std::vector<uint32_t>* rank,
        std::vector<uint32_t>* part)
{
    g_ = g;
    part_ = part;
    rank_ = rank;

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
    if((nparts_ & 7) > 0)
    {
        bytes_per_label_++;
    }

    // allocate memory for arcflag labels. 
    for(uint32_t i = 0; i < g_->get_num_nodes(); i++)
    {
        labels_.push_back(std::vector<bbaf_label>());
        warthog::graph::node* n = g_->get_node(i);
        for(uint32_t j = 0; j < n->out_degree(); j++)
        {
            bbaf_label label;
            label.flags_ = new uint8_t[bytes_per_label_];
            labels_.back().push_back(label);

            // each bit of each label is initialised to zero
            for(uint32_t k = 0; k < bytes_per_label_; k++)
            {
                label.flags_[k] = 0;
            }
        }
    }
    t_byte_ = 0;
    t_bitmask_ = 0;
    tx_ = ty_ = 0;
}

void
warthog::bb_af_filter::print(std::ostream& out)
{
    out 
    << "# Each line (aside from comments & header) is an edge label.\n"
    << "# A label comprises a set of arcflags (written out as 64bit words) \n"
    << "# and a rectangular bounding box which contains all nodes in the \n"
    << "# down-closure of the associated edge. The last 4 digits of each\n"
    << "# line describe the bounding box: minx, miny, maxx, maxy.\n"
    << "#\n"
    << "# NB: Before printing, labels are sorted by edge-tail index. \n"
    << "# Ties are broken using the order edges appear in the file \n"
    << "# " << g_->get_filename() << std::endl
    << "# (lower/earlier is better)\n";

    out << "firstnode " << firstid_ << " lastnode " << lastid_
        << " partitions " << nparts_ << std::endl;

    // we split labels into 64bit words for printing
    uint32_t words_per_label = bytes_per_label_ / 8;
    if((bytes_per_label_ % 8) > 0) words_per_label++;

    // iterate over the labels for each outgoing arc of each node
    for(uint32_t i = firstid_; i <= lastid_; i++)
    {
        std::vector<bbaf_label>& arclabs = labels_.at(i);
        for(uint32_t j = 0; j < arclabs.size(); j++)
        {
            // print the arcflags
            bbaf_label& label = arclabs.at(j);
            for(uint32_t word = 0; word < words_per_label; word++)
            {
                uint8_t printed_word[8];
                for(uint32_t k = 0; k < 8; k++)
                {
                    // read the flag label, byte by byte, one word at a time
                    if((word*8+k) < bytes_per_label_)
                    {
                        printed_word[k] = label.flags_[word*8+k];
                    }
                    // pad the last word with leading zeroes if necessary
                    else
                    {
                        printed_word[k] = 0;
                    }
                }
                out << *((uint64_t*)&(printed_word)) << " ";
            }

            // print the bounding box
            out << label.bbox_.x1 << " " << label.bbox_.y1 << " " 
                << label.bbox_.y1 << " " << label.bbox_.y2;

            out << std::endl;
        }
    }
}

bool
warthog::bb_af_filter::load_labels(const char* filename)
{
    std::cerr << "loading bb_af file\n";
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
                    << "err reading arcflags file; number of partitions \n"
                    << " does not match input graph (read: " << num_parts 
                    << " expected: " << nparts_ << ")\n";
                exit(0);
            }
        }
        else { bad_header = true; }
    }
    else { bad_header = true; }

    if(bad_header)
    {
        std::cerr << "invalid header; required format: \n"
            << "firstnode [integer] lastnode [integer] partitions [integer]\n";
        return false;
    }

    // read labels for each outgoing arc
    uint32_t words_per_label = (bytes_per_label_ / 8);
    if((bytes_per_label_ % 8) != 0) { bytes_per_label_++; }
    uint32_t lines = 0;
    for(uint32_t i = firstid_; i < (lastid_+1); i++)
    {
        warthog::graph::node* n = g_->get_node(i);
        for(uint32_t j = 0; j < n->out_degree(); j++)
        {
            for(uint32_t word = 0; word < words_per_label; word++)
            {
                if(!ifs.good())
                {
                    std::cerr 
                        << "unexpected error reading arcflags data on line "
                        << lines << "; aborting\n";
                    return false;
                }

                uint64_t label;
                ifs >> label;
                *((uint64_t*)&labels_.at(i).at(j).flags_[word*8]) = label;
                assert(*((uint64_t*)&labels_.at(i).at(j).flags_[word*8])
                            == label);
            }

            if(!ifs.good())
            {
                std::cerr 
                    << "unexpected error reading arcflags data on line "
                    << lines << "; aborting\n";
                return false;
            }

            uint32_t x1, y1, x2, y2;
            ifs >> x1 >> y1 >> x2 >> y2;
            labels_.at(i).at(j).bbox_.x1 = x1;
            labels_.at(i).at(j).bbox_.y1 = y1;
            labels_.at(i).at(j).bbox_.x2 = x2;
            labels_.at(i).at(j).bbox_.y2 = y2;
            lines++;
       }
    }
    return true;
}

void
warthog::bb_af_filter::compute()
{

    compute(0, g_->get_num_nodes()-1);

}

// compute arclabels for all nodes whose ids are in
// the range [firstid, lastid]
void
warthog::bb_af_filter::compute(uint32_t firstid, uint32_t lastid)
{
    if(!g_ || !rank_) { return; } 

    firstid_ = firstid;
    lastid_ = lastid;
    if(lastid_ >= g_->get_num_nodes()) { lastid_ = g_->get_num_nodes() - 1; }

    std::cerr << "computing ch arcflag labels " 
        << "for all nodes in the id-range [" 
        << firstid_ << ", " << lastid_ << "]\n";
    warthog::zero_heuristic heuristic;
    warthog::fwd_ch_expansion_policy expander(g_, rank_);

    warthog::flexible_astar<
        warthog::zero_heuristic,
        warthog::fwd_ch_expansion_policy> dijkstra(&heuristic, &expander);

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
        dijkstra.get_length(source_id, warthog::INF);

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
                    uint32_t e_idx  = (*idmap.find(
                            n->get_parent()->get_id() == source_id ? 
                            n->get_id() : n->get_parent()->get_id())).second;
                    labels_.at(source_id).at(e_idx).flags_[part_id >> 3]
                        |= (1 << (part_id & 7));

                    // only nodes that are in the down-closure get added
                    // to the bounding box
                    if(this->rank_->at(n->get_id()) < rank_->at(source_id))
                    {
                        int32_t x, y;
                        this->g_->get_xy(n->get_id(), x, y);
                        assert(x != warthog::INF && y != warthog::INF);
                        labels_.at(source_id).at(e_idx).bbox_.grow(x, y);
                        assert(labels_.at(source_id).at(e_idx).bbox_.is_valid());
                    }

                };
        dijkstra.apply_to_closed(fn_arcflags);
    }
    std::cerr << "\nall done\n"<< std::endl;
}

void
warthog::bb_af_filter::set_goal(uint32_t goalid)
{
    uint32_t t_part = part_->at(goalid);
    t_byte_ = t_part >> 3;
    t_bitmask_ = 1 << (t_part & 7);
    g_->get_xy(goalid, tx_, ty_);
}
