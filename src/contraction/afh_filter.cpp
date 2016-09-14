#include "afh_filter.h"
#include "constants.h"
#include "contraction.h"
#include "flexible_astar.h"
#include "fwd_ch_expansion_policy.h"
#include "planar_graph.h"
#include "problem_instance.h"
#include "zero_heuristic.h"

#include <fstream>
#include <iostream>
#include <set>
#include <unordered_map>

warthog::afh_filter::afh_filter(
        warthog::graph::planar_graph* g,
        std::vector<uint32_t>* rank,
        std::vector<uint32_t>* part)
{
    init(g, rank, part);
}

warthog::afh_filter::afh_filter(
        warthog::graph::planar_graph* g,
        std::vector<uint32_t>* rank,
        std::vector<uint32_t>* part,
        const char* afh_file)
{
    init(g, rank, part);
    load_afh_file(afh_file);
}

warthog::afh_filter::~afh_filter()
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
warthog::afh_filter::init(
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

    // allocate memory for arcflag labels
    for(uint32_t i = 0; i < g_->get_num_nodes(); i++)
    {
        // allocate memory for edge labels and initialise
        // each label to zero
        flags_.push_back(std::vector<uint8_t*>());
        warthog::graph::node* n = g_->get_node(i);
        for(uint32_t j = 0; j < n->out_degree(); j++)
        {
            uint8_t* label = new uint8_t[bytes_per_label_];
            flags_.back().push_back(label);

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
warthog::afh_filter::print(std::ostream& out)
{
    out << "# to save some space, labels are written out using 64-bit words\n";
    out << "# with spaces between words for labels longer than 64 bits\n";

    out << "firstnode " << firstid_ << " lastnode " << lastid_
        << " partitions " << nparts_ << std::endl;

    // we split labels into 64bit words for printing
    uint32_t words_per_label = bytes_per_label_ / 8;
    if((bytes_per_label_ % 8) > 0) words_per_label++;

    // iterate over the labels for each outgoing arc of each node
    for(uint32_t i = firstid_; i <= lastid_; i++)
    {
        std::vector<uint8_t*>& node_flags = flags_.at(i);
        for(uint32_t j = 0; j < node_flags.size(); j++)
        {
            uint8_t* label = node_flags.at(j);
            for(uint32_t word = 0; word < words_per_label; word++)
            {
                uint8_t printed_word[8];
                for(uint32_t k = 0; k < 8; k++)
                {
                    // read the label, byte by byte, one word at a time
                    if((word*8+k) < bytes_per_label_)
                    {
                        printed_word[k] = label[word*8+k];
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

void
warthog::afh_filter::load_afh_file(const char* filename)
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
        exit(0);
    }

    // read labels for each outgoing arc
    uint32_t words_per_label = (bytes_per_label_ / 8);
    if((bytes_per_label_ % 8) != 0) { bytes_per_label_++; }
    for(uint32_t i = firstid_; i < (lastid_+1); i++)
    {
        warthog::graph::node* n = g_->get_node(i);
        for(uint32_t j = 0; j < n->out_degree(); j++)
        {
            for(uint32_t word = 0; word < words_per_label; word++)
            {
                if(!ifs.good())
                {
                    std::cerr << "unexpected read error loading arcflags; "
                        << "aborting\n";
                    exit(0);
                }

                uint64_t label;
                ifs >> label;
                *((uint64_t*)&flags_.at(i).at(j)[word*8]) = label;
                assert(*((uint64_t*)&flags_.at(i).at(j)[word*8])
                            == label);
            }
        }
    }
}

void
warthog::afh_filter::compute()
{

    compute(0, g_->get_num_nodes()-1);

}

// compute arclabels for all nodes whose ids are in
// the range [firstid, lastid]
void
warthog::afh_filter::compute(uint32_t firstid, uint32_t lastid)
{
    if(!g_ || !rank_) { return; } 

    firstid_ = firstid;
    lastid_ = lastid;
    if(lastid_ >= g_->get_num_nodes()) { lastid_ = g_->get_num_nodes() - 1; }

    // create an array with node ids ordered by rank
    std::vector<uint32_t> ids_by_rank(*rank_); 
    warthog::ch::value_index_swap_dimacs(ids_by_rank);

    std::cerr << "computing afh down-labels\n";
    compute_down_flags(ids_by_rank);

    std::cerr << "computing afh up-labels\n";
    compute_up_flags(ids_by_rank);

    std::cerr << "\nall done\n"<< std::endl;
}

void
warthog::afh_filter::compute_up_flags(std::vector<uint32_t>& ids_by_rank)
{
    // traverse the hierarchy top-to-bottom
    for(int32_t idx = ids_by_rank.size()-1; idx >= 0; idx--)
    {
        std::cerr << "\rprocessing node ranked " << idx 
            << "; continues until rank 0";
        uint32_t n_id = ids_by_rank.at(idx);

        // iterate over up edges leading into node n;
        // i.e. (m, n) \in E | n > m
        warthog::graph::node* n = g_->get_node(n_id);
        for(warthog::graph::edge_iter in_it = n->incoming_begin(); 
                in_it != n->incoming_end(); in_it++)
        {
            // ignore edges where m > n (we want to label up only)
            uint32_t m_id = (*in_it).node_id_;
            if(rank_->at(m_id) > rank_->at(n_id)) { continue;  }
            warthog::graph::node* m = g_->get_node(m_id);

            // identify the outgoing version of (m, n)
            warthog::graph::edge_iter it_e_mn = m->find_edge(n_id);
            uint32_t e_mn_idx = it_e_mn - m->outgoing_begin();
            assert(it_e_mn != m->outgoing_end());

            std::set<uint32_t> intermediate;
            warthog::ch::unpack(m_id, it_e_mn, g_, intermediate);
            for(warthog::graph::edge_iter out_it = n->outgoing_begin(); 
                    out_it != n->outgoing_end(); out_it++)
            {
                uint32_t n_nei_id = (*out_it).node_id_;

                // <m, n, m> is redundant
                if(n_nei_id == m_id) { continue; }

                // <m, n, n_nei> is redundant due to edge (m, n_nei)
//                if(m->find_edge(n_nei_id) != m->outgoing_end()) { continue; }
                warthog::graph::edge_iter it_m_nei = m->find_edge(n_nei_id);
                if( it_m_nei != m->outgoing_end() && 
                    (*it_m_nei).wt_ <= ((*in_it).wt_ + (*out_it).wt_))
                { continue; }

                // <m, n, n_nei> is redundant if n_nei is bypassed by the
                // edge (m, n); i.e. a down path exists <m, ..., n_nei>
                if(intermediate.find(n_nei_id) == intermediate.end()) 
                {
                    // add to the label of (m, n) the flags of (n, n_nei)
                    for(uint32_t i = 0; i < bytes_per_label_; i++)
                    {
                        flags_.at(m_id).at(e_mn_idx)[i] |=
                            flags_.at(n_id).at(out_it - n->outgoing_begin())[i];
                    }
                }
            }

            // add the colour of n to the label of (m, n)
            uint32_t n_col = part_->at(n_id);
            flags_.at(m_id).at(e_mn_idx)[n_col >> 3] |= (1 << (n_col & 7));
        }
    }
    std::cerr <<"\ndone\n"<< std::endl;
}

//void
//warthog::afh_filter::unpack(uint32_t from_id,
//        warthog::graph::edge_iter it_e,
//        std::set<uint32_t>& intermediate)
//{
//    warthog::graph::node* from = g_->get_node(from_id);
//    assert(it_e >= from->outgoing_begin() && it_e < from->outgoing_end());
//
//    warthog::graph::edge* e_ft = &*it_e;
//    uint32_t to_id = e_ft->node_id_;
//
//    for(warthog::graph::edge_iter it = from->outgoing_begin(); 
//            it < from->outgoing_end(); it++)
//    {
//        warthog::graph::node* succ = g_->get_node((*it).node_id_);
//        warthog::graph::edge_iter it_e_succ = succ->outgoing_begin();
//        while(true)
//        {
//            it_e_succ = succ->find_edge(to_id, it_e_succ);
//            if(it_e_succ == succ->outgoing_end()) { break; }
//            assert( it_e_succ >= succ->outgoing_begin() && 
//                    it_e_succ <= succ->outgoing_end());
//
//            warthog::graph::edge* e_st = &*it_e_succ;
//            if(((*it).wt_ + e_st->wt_) == e_ft->wt_) { break; }
//            it_e_succ++;
//        }
//        if(it_e_succ == succ->outgoing_end()) { continue; } 
//
//        // recursively unpack the two edges being represented by
//        // the single shortcut edge (from_id, to_id)
//        unpack(from_id, it, intermediate);
//
//        succ = g_->get_node((*it).node_id_);
//        unpack((*it).node_id_, it_e_succ, intermediate);
//        break;
//    }
//}

void
warthog::afh_filter::compute_down_flags(std::vector<uint32_t>& ids_by_rank)
{
    // compute down-flags. 
    for(uint32_t idx = 0; idx < ids_by_rank.size(); idx++)
    {
        std::cerr << "\rprocessing node ranked " << idx 
            << "; continues until rank " << (ids_by_rank.size()-1);
        uint32_t n_id = ids_by_rank.at(idx);

        warthog::graph::node* n = g_->get_node(n_id);
        uint32_t n_col = part_->at(n_id);

        // label each down-edge, (m, n), being careful to avoid
        // redundant paths such as <m, n, q> when <m, q> exists
        for(warthog::graph::edge_iter ei_mn = n->incoming_begin(); 
            ei_mn != n->incoming_end(); ei_mn++)
        {
            // we focus on down edges where m > n 
            uint32_t m_id = (*ei_mn).node_id_;
            if(rank_->at(m_id) < rank_->at(n_id)) { continue; }

            // to avoid redundant paths we build a set with down-successors
            // of m which might also be down-successors of n
            warthog::graph::node* m = g_->get_node(m_id);
            std::unordered_map<uint32_t, warthog::graph::edge> m_down;
            for(warthog::graph::edge_iter em = m->outgoing_begin();
                    em != m->outgoing_end(); em++)
            {
                if(rank_->at((*em).node_id_) > rank_->at(n_id)) { continue; }
                m_down.insert(
                        std::pair<uint32_t, warthog::graph::edge>(
                            (*em).node_id_, *em));
            }

            // next we construct a label for the edge (m, n) by taking
            // a union of arcflags from the non-redundant down successors of n 
            warthog::graph::edge_iter it_e_mn = m->find_edge(n_id); 
            assert( it_e_mn >= m->outgoing_begin() && 
                    it_e_mn < m->outgoing_end());
            int32_t e_mn_index = it_e_mn - m->outgoing_begin();
            uint8_t* e_mn_flags = flags_.at(m_id).at(e_mn_index);

            for(warthog::graph::edge_iter en = n->outgoing_begin(); 
                    en != n->outgoing_end(); en++)
            {
                // only down successors
                if(rank_->at((*en).node_id_) > rank_->at(n_id)) { continue; }

                // only if the down path <m, n, n'> is not redundant
                std::unordered_map<uint32_t, warthog::graph::edge>::iterator mit
                    = m_down.find((*en).node_id_);
                if(mit != m_down.end() &&
                        (*mit).second.wt_ < ((*ei_mn).wt_ + (*en).wt_)) { continue; }

                // label (m, n) with colours that are reachable via 
                // non-redundant edges to the down-successors of n
                uint32_t en_idx = en - n->outgoing_begin();
                for(uint8_t j = 0; j < bytes_per_label_; j++)
                { e_mn_flags[j] |= flags_.at(n_id).at(en_idx)[j]; }
            }
            //finally, label (m, n) with the colour of n
            e_mn_flags[n_col >> 3] |= (1 << (n_col & 7));
        }
    }
    std::cerr << "\ndone\n";
}
