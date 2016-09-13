#include "afhd_filter.h"
#include "ch_expansion_policy.h"
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

warthog::afhd_filter::afhd_filter(
        warthog::graph::planar_graph* g,
        std::vector<uint32_t>* rank,
        std::vector<uint32_t>* part)
{
    init(g, rank, part);
}

warthog::afhd_filter::afhd_filter(
        warthog::graph::planar_graph* g,
        std::vector<uint32_t>* rank,
        std::vector<uint32_t>* part,
        const char* afhd_file)
{
    init(g, rank, part);
    load_afhd_file(afhd_file);
}

warthog::afhd_filter::~afhd_filter()
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
warthog::afhd_filter::init(
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
warthog::afhd_filter::print(std::ostream& out)
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
warthog::afhd_filter::load_afhd_file(const char* filename)
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
warthog::afhd_filter::compute()
{

    compute(0, g_->get_num_nodes()-1);

}

// compute arclabels for all nodes whose ids are in
// the range [firstid, lastid]
void
warthog::afhd_filter::compute(uint32_t firstid, uint32_t lastid)
{
    if(!g_ || !rank_) { return; } 

    firstid_ = firstid;
    lastid_ = lastid;
    if(lastid_ >= g_->get_num_nodes()) { lastid_ = g_->get_num_nodes() - 1; }

    // create an array with node ids ordered by rank
    std::vector<uint32_t> ids_by_rank(*rank_); 
    warthog::ch::value_index_swap_dimacs(ids_by_rank);

    std::cerr << "computing afhd down-labels\n";
    compute_down_flags(ids_by_rank);

    std::cerr << "computing afhd up-labels\n";
    compute_up_flags(ids_by_rank);

    std::cerr << "\nall done\n"<< std::endl;
}

//void
//warthog::afhd_filter::compute_up_flags()
//{
//    // traverse the hierarchy top-to-bottom
//    std::cerr << "computing afhd up-labels\n";
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
warthog::afhd_filter::compute_up_flags(std::vector<uint32_t>& ids_by_rank)
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
            unpack(m_id, it_e_mn, intermediate);
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

void
warthog::afhd_filter::unpack(uint32_t from_id,
        warthog::graph::edge_iter it_e,
        std::set<uint32_t>& intermediate)
{
    warthog::graph::node* from = g_->get_node(from_id);
    assert(it_e >= from->outgoing_begin() && it_e < from->outgoing_end());

    warthog::graph::edge* e_ft = &*it_e;
    uint32_t to_id = e_ft->node_id_;

    for(warthog::graph::edge_iter it = from->outgoing_begin(); 
            it < from->outgoing_end(); it++)
    {
        warthog::graph::node* succ = g_->get_node((*it).node_id_);
        warthog::graph::edge_iter it_e_succ = succ->outgoing_begin();
        while(true)
        {
            it_e_succ = succ->find_edge(to_id, it_e_succ);
            if(it_e_succ == succ->outgoing_end()) { break; }
            assert( it_e_succ >= succ->outgoing_begin() && 
                    it_e_succ <= succ->outgoing_end());

            warthog::graph::edge* e_st = &*it_e_succ;
            if(((*it).wt_ + e_st->wt_) == e_ft->wt_) { break; }
            it_e_succ++;
        }
        if(it_e_succ == succ->outgoing_end()) { continue; } 

        // recursively unpack the two edges being represented by
        // the single shortcut edge (from_id, to_id)
        unpack(from_id, it, intermediate);

        succ = g_->get_node((*it).node_id_);
        unpack((*it).node_id_, it_e_succ, intermediate);
        break;
    }
}

void
warthog::afhd_filter::compute_down_flags(std::vector<uint32_t>& ids_by_rank)
{
    warthog::zero_heuristic heuristic;
    warthog::ch_expansion_policy expander(g_, rank_, false, warthog::ch::DOWN);

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

    // label each down-edge, (m, n), being careful to avoid
    // redundant paths such as <m, n, q> when <m, q> exists
    for(uint32_t i = 0; i <= lastid_; i++)
    {
        std::cerr << "\rprocessing node " << i << "; continues"
            << " until node " << lastid_;
        uint32_t source_id = i;
        dijkstra.get_length(source_id, warthog::INF);
        
        // create a rectangle label for every outgoing edge
        warthog::graph::node* source = g_->get_node(source_id);

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

        // label each arc, e, of the source with the colour of every node 
        // that is reached by an optimal down path whose first edge is e
        std::function<void(warthog::search_node*)> label_fn = 
                [this, &source_id, &idmap](warthog::search_node* n) -> void
                {
                    // skip the source (it doesn't appear on any path)
                    assert(n);
                    if(n->get_id() == source_id) { return; } 
                    assert(n->get_parent());

                    // identify the first arc on the optimal down path 
                    // to node n
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
                    uint32_t n_col = part_->at(n->get_id());
                    flags_.at(source_id).at((*it).second)[n_col >> 3] 
                        |= (1 << (n_col & 7));
                };
        dijkstra.apply_to_closed(label_fn);
    }
    std::cerr << "\ndone\n";
}
