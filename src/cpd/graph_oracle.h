// cpd/graph_oracle.h
//
// A Compressed Path Database is an oracle that tells
// the first move on the optimal path: from any source
// node in graph V, to any other node, in graph V. 
// The implementation is based on descriptions from 
// the following paper:
//
// [A. Botea, B. Strasser and D. Harabor. 2015. 
//  Compressing Optimal Paths with Run Length Encoding.
//  Journal of Artificial Intelligence Research (JAIR)]
//
//
// @author: dharabor
// @created: 2020-02-26
//

#ifndef WARTHOG_CPD_GRAPH_ORACLE_H
#define WARTHOG_CPD_GRAPH_ORACLE_H

#include "constants.h"
#include "cpd.h"
#include "graph_expansion_policy.h"
#include "xy_graph.h"

namespace warthog
{

namespace cpd
{

class graph_oracle
{
    public: 
        graph_oracle(warthog::graph::xy_graph* g) 
             : g_(g)
        {
            order_.resize(g_->get_num_nodes());
            fm_.resize(g_->get_num_nodes());
        }

        graph_oracle() : g_(nullptr) { }

        virtual ~graph_oracle() { } 

        graph_oracle(const graph_oracle&) = default;

        bool
        operator==(const graph_oracle& other)
        {
            if (order_ != other.order_)
            {
                return false;
            }

            if (fm_.size() != other.fm_.size())
            {
                return false;
            }

            for (size_t i = 0; i < fm_.size(); i++)
            {
                std::vector<warthog::cpd::rle_run32> row1 = fm_.at(i);
                std::vector<warthog::cpd::rle_run32> row2 = other.fm_.at(i);

                if (row1.size() != row2.size())
                {
                    return false;
                }

                for (size_t j = 0; j < row1.size(); j++)
                {
                    if (row1.at(j).data_ != row2.at(j).data_)
                    {
                        return false;
                    }
                }
            }

            return true;
        }

        inline uint32_t 
        get_move(warthog::sn_id_t source_id, 
                 warthog::sn_id_t target_id)
        {
            if(fm_.at(source_id).size() == 0) { return warthog::cpd::CPD_FM_NONE; }

            std::vector<warthog::cpd::rle_run32>& row = fm_.at(source_id);
            uint32_t target_index = order_.at(target_id);
            uint32_t end = (uint32_t)row.size();
            uint32_t begin = 0;
            while(begin<(end-1))
            {
                uint32_t mid = begin + ((end-begin)>>1);
                if(target_index < row.at(mid).get_index()) { end = mid ;  }
                else { begin = mid; }
            }
            return row.at(begin).get_move();
        }

        inline void
        clear()
        {
            order_.clear();
            fm_.clear();
        }

        inline void
        compute_dfs_preorder()
        {
            warthog::cpd::compute_dfs_preorder(g_, &order_);
        }

        // convert the column order into a map: from vertex id to its ordered
        // index
        inline void
        value_index_swap_array()
        {
            warthog::helpers::value_index_swap_array(order_);
        }

        // compress a given first-move table @param row and associate
        // the compressed result with source node @param source_id
        void
        add_row(uint32_t source_id, 
                 std::vector<warthog::cpd::fm_coll>& row);

        inline warthog::graph::xy_graph* 
        get_graph() { return g_; } 

        inline void
        set_graph(warthog::graph::xy_graph* g)
        { g_ = g; }

        inline size_t
        mem()
        {
            size_t retval = 
                g_->mem() + 
                sizeof(uint32_t) * order_.size() + 
                sizeof(std::vector<warthog::cpd::rle_run32>) * fm_.size();

            for(uint32_t i = 0; i < fm_.size(); i++)
            {
                retval += sizeof(warthog::cpd::rle_run32) * fm_.at(i).size();
            }

            return retval; 
        }

        friend std::ostream&
        operator<<(std::ostream& out, warthog::cpd::graph_oracle& o);

        friend std::istream&
        operator>>(std::istream& in, warthog::cpd::graph_oracle& o);

        warthog::cpd::graph_oracle&
        operator+=(const warthog::cpd::graph_oracle &cpd);

    private:
        std::vector<std::vector<warthog::cpd::rle_run32>> fm_;
        std::vector<uint32_t> order_;
        warthog::graph::xy_graph* g_;
};

std::ostream&
operator<<(std::ostream& out, warthog::cpd::graph_oracle& o);

std::istream&
operator>>(std::istream& in, warthog::cpd::graph_oracle& o);

void
compute_row(uint32_t source_id, warthog::cpd::graph_oracle* cpd,
            warthog::search* dijk, std::vector<warthog::cpd::fm_coll> &s_row);

// helps to precompute first-move data
struct graph_oracle_listener
{
    inline void
    generate_node(warthog::search_node *from, warthog::search_node *succ,
                  warthog::cost_t edge_cost, uint32_t edge_id)

    {
        if(from == 0) { return; } // start node

        if(from->get_id() == *source_id_) // start node successors
        {
            //assert(s_row_.at(succ->get_id()) == 0);
            assert(edge_id <
                        oracle_->get_graph()->get_node(
                        (uint32_t)*source_id_)->out_degree());
            s_row_->at(succ->get_id()) = (1 << edge_id);
            assert(s_row_->at(succ->get_id()));
        }
        else // all other nodes
        {
            warthog::sn_id_t succ_id = succ->get_id();
            warthog::sn_id_t from_id = from->get_id();
            double alt_g = from->get_g() + edge_cost;
            double g_val =
                succ->get_search_number() == from->get_search_number() ?
                succ->get_g() : DBL_MAX;

            //  update first move
            if(alt_g < g_val)
            {
                s_row_->at(succ_id) = s_row_->at(from_id);
                assert(s_row_->at(succ_id) == s_row_->at(from_id));
            }

            // add to the list of optimal first moves
            if(alt_g == g_val)
            {
                s_row_->at(succ_id) |= s_row_->at(from_id);
                assert(s_row_->at(succ_id) >= s_row_->at(from_id));
            }

        }
    }

    inline void
    expand_node(warthog::search_node* current) { }

    inline void
    relax_node(warthog::search_node* current) { }

    warthog::cpd::graph_oracle* oracle_;
    warthog::sn_id_t* source_id_;
    std::vector<warthog::cpd::fm_coll>* s_row_;

};

}

}

#endif
