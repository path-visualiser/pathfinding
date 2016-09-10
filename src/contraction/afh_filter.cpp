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

//void
//warthog::afh_filter::compute_up_flags()
//{
//    // traverse the hierarchy top-to-bottom
//    std::cerr << "computing afh up-labels\n";
//    for(uint32_t source_id = last_id_; source_id >= firstid_; source_id--)
//    {
//        std::cerr << "\rprocessing node " << i << "; continues until node " 
//            << firstid_ << "\r";
//
//        std::set<uint32_t> neis;
//        std::vector<std::set<uint32_t>> reachable(nparts_);
//        warthog::search_node* source = g_->get_node(source_id);
//
//        // scan the up and down flags of every successor node and
//        // record which node is useful for reaching which colour
//        for(warthog::graph::edge_iter it = n->outgoing_begin(); 
//                it != n->outgoing_end(); it++)
//        {
//            uint32_t e_idx = it - n->outgoing_begin();
//            uint32_t nei_idx = (*it).node_id_;
//            neis.insert(nei_idx);
//            for(uint32_t col = 0; col < nparts_; col++)
//            {
//                uint8_t* eflag = flags_->at(nei_idx).at(e_idx);
//                if(eflag[col >> 3] & (1 << (col & 7)))
//                {
//                    reachable.at(col).insert(nei_idx);
//                }
//            }
//        }
//
//        // label up edges leading into the current source node
//        std::set<uint32_t> mycolours;
//        for(warthog::graph::edge_iter it = n->incoming_begin(); 
//                it != n->incoming_end(); it++)
//        {
//            uint32_t nei_idx = (*it).node_id_;
//            warthog::graph::node* nei = g_->get_node(nei_idx);
//            uint32_t e_idx = nei->find_edge_index(source_id);
//            assert(e_idx != warthog::INF);
//            // process colours we might be able to reach optimally
//            // by taking an up edge (to reach source) and then a down
//            // edge (from source to one of its successors)
//            for(uint32_t col = 0; col < nparts_; col++)
//            {
//                std::set<uint32_t>& col_neis = dn_reach.at(col);
//
//                // skip colours we cannot reach optimally
//                if(col_neis.size() == 0) { continue; } 
//
//                // always label the edge with the colour of 
//                // the source
//                if(col == part_->at(source_id))
//                {
//                    mycolours.insert(col);
//                    continue;
//                }
//
//                if(col_neis.find(nei_idx) != col_neis.end())
//                {
//                    // redundant (i): <nei_idx, source, nei_idx>
//                    if(col_neis.size() == 1) { continue; }
//
//                    // redundant (ii): we try to prove for each up/down 
//                    // path such as <nei_idx, source, nei2_idx>
//                    // that there is shorter alternative, namely
//                    // <nei_idx, nei2_idx>
//                    //
//                    // we do two things here:
//                    // 1. we try to look for a direct edge to nei2_idx
//                    // 2. we unpack shortcut edges and look for nei2_idx 
//                    // among the intermediate nodes on that path 
//                    std::set<uint32_t> intermediate;
//                    unpack(source_id, nei_idx, intermediate);
//                    for(uint32_t x = 0; x < col_neis.size() x++) 
//                    { 
//                        uint32_t nei2_idx = *(col_neis.begin() + x);
//
//                        // paths that backtrack on themselves are redundant
//                        // i.e. <nei_idx, source, nei_idx>
//                        if(nei_idx == nei2_idx) { continue; }
//
//                        // check if the two nodes are directly connected
//                        if(find_outgoing_edge(nei_idx, nei2_idx) 
//                                != warthog::INF) { continue; }
//
//                        // shortcut unpacking
//                        if(intermediate.find(x) == intermediate.end())
//                        {
//                            // it seems there is at least one node that could
//                            // be optimally reachanble via the up,down path 
//                            // <nei_idx, source, nei2_idx>
//                            // we thus label (nei_idx, source) with col
//                            mycolours.insert(col);
//                            break;
//                        }
//                    }
//                }
//                else
//                {
//                    // col doesn't appear in dn_flags_[nei_idx]
//                    // we need to pass through node n to reach it
//                    mycolours.push_back(col);
//                }
//            }
//        }
//    }
//    std::cerr <<"\ndone\n"<< std::endl;
//}

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
        for(warthog::graph::edge_iter it = n->incoming_begin(); 
                it != n->incoming_end(); it++)
        {
            // ignore edges where m > n (we want to label up only)
            uint32_t m_id = (*it).node_id_;
            if(rank_->at(m_id) > rank_->at(n_id)) { continue;  }
            warthog::graph::node* m = g_->get_node(m_id);

            // identify the outgoing version of (m, n)
            uint32_t e_mn_idx = m->find_edge_index(n_id);
            assert(e_mn_idx != warthog::INF);

            std::set<uint32_t> intermediate;
            unpack(m_id, e_mn_idx, intermediate);
            for(warthog::graph::edge_iter it = n->outgoing_begin(); 
                    it != n->outgoing_end(); it++)
            {
                uint32_t n_nei_id = (*it).node_id_;

                // <m, n, m> is redundant
                if(n_nei_id == m_id) { continue; }

                // <m, n, n_nei> is redundant due to edge (m, n_nei)
                if(m->find_edge_index(n_nei_id) != warthog::INF) { continue; }

                // <m, n, n_nei> is redundant if n_nei is bypassed by the
                // edge (m, n); i.e. a down path exists <m, ..., n_nei>
                if(intermediate.find(n_nei_id) == intermediate.end()) 
                {
                    // add to the label of (m, n) the flags of (n, n_nei)
                    for(uint32_t i = 0; i < bytes_per_label_; i++)
                    {
                        flags_.at(m_id).at(e_mn_idx)[i] |=
                            flags_.at(n_id).at(it - n->outgoing_begin())[i];
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

void
warthog::afh_filter::unpack(uint32_t from_id, uint32_t edge_idx,
        std::set<uint32_t>& intermediate)
{
    warthog::graph::node* from = g_->get_node(from_id);
    assert(edge_idx < from->out_degree());

    warthog::graph::edge* e_ft = &*(from->outgoing_begin() + edge_idx);
    uint32_t to_id = e_ft->node_id_;

    for(warthog::graph::edge_iter it = from->outgoing_begin(); 
            it < from->outgoing_end(); it++)
    {
        warthog::graph::node* succ = g_->get_node((*it).node_id_);
        uint32_t e_succ_idx = 0;
        while(true)
        {
            e_succ_idx = succ->find_edge_index(to_id, e_succ_idx);
            if(e_succ_idx == warthog::INF) { break; }
            assert(e_succ_idx < succ->out_degree());

            warthog::graph::edge* e_st = 
                &*(succ->outgoing_begin() + e_succ_idx);
            if(((*it).wt_ + e_st->wt_) == e_ft->wt_) { break; }
            e_succ_idx++;
        }
        if(e_succ_idx == warthog::INF) { continue; } 

        // recursively unpack the two edges being represented by
        // the single shortcut edge (from_id, to_id)
        unpack(from_id, it - from->outgoing_begin(), intermediate);

        succ = g_->get_node((*it).node_id_);
        assert(e_succ_idx < succ->out_degree());
        unpack((*it).node_id_, e_succ_idx, intermediate);
        break;
    }
}

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
            uint32_t e_mn_index = m->find_edge_index(n_id);
            assert(e_mn_index != warthog::INF);
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
