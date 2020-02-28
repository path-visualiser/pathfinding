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

        virtual ~graph_oracle() { } 

        inline uint32_t 
        get_move(warthog::sn_id_t source_id, 
                 warthog::sn_id_t target_id)
        {
            assert(source_id < g_->get_num_nodes());
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
    
        void
        add_row(uint32_t source_id, 
                 std::vector<warthog::cpd::fm_coll>& row)
        {
            // source gets a wildcard move
            row.at(source_id) = warthog::cpd::CPD_FM_NONE;

            // greedily compress the row w.r.t. the current column order
            warthog::cpd::fm_coll current = row.at(order_.at(0));
            uint32_t head = 0;
            for(uint32_t index = 0; index < row.size(); index++)
            {
                assert(current > 0);
                uint32_t next_id = order_.at(index);
                if(!(current & row.at(next_id)))
                {
                    uint32_t firstmove = __builtin_ffsl(current) - 1;
                    assert(firstmove < warthog::cpd::CPD_FM_MAX);
                    fm_.at(source_id).push_back(
                            warthog::cpd::rle_run32{ (head << 4) | firstmove} );
                    current = row.at(next_id);
                    head = index;
                }
                current = current & row.at(next_id);
            } 
            
            // add the last run
            uint32_t firstmove = __builtin_ffsl(current) - 1;
            assert(firstmove < warthog::cpd::CPD_FM_MAX);
            fm_.at(source_id).push_back(
                    warthog::cpd::rle_run32{ (head << 4) | firstmove} );
            
            std::cerr << "compressed source row " << source_id << " with " << fm_.at(source_id).size() << std::endl;
        }

        void
        precompute() 
        {
            warthog::cpd::compute_dfs_preorder(g_, &order_);
            warthog::simple_graph_expansion_policy expander(g_);

            std::vector<uint32_t> source_nodes;
            source_nodes.reserve(g_->get_num_nodes());
            for(uint32_t i = 0; i < g_->get_num_nodes(); i++)
            {
                source_nodes.push_back(i);
            }

            warthog::cpd::compute_firstmoves_and_compress<
                warthog::simple_graph_expansion_policy, 
                warthog::cpd::graph_oracle> (this, &expander, &source_nodes);
        }

        warthog::graph::xy_graph*
        get_graph() { return g_; } 


    private:
        std::vector<std::vector<warthog::cpd::rle_run32>> fm_;
        std::vector<uint32_t> order_;
        warthog::graph::xy_graph* g_;
};

//friend ostream&
//operator<<(std::ostream& out, warthog::cpd::graph_oracle& o);
//
//friend ostream&
//operator<<(std::ostream& in, warthog::cpd::graph_oracle& o);

}

}

#endif

