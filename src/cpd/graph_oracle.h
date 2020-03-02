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

        // compute first-move labels and compress the result
        void
        precompute();
    
        // compress a given first-move table @param row and associate
        // the compressed result with source node @param source_id
        void
        add_row(uint32_t source_id, 
                 std::vector<warthog::cpd::fm_coll>& row);

        inline warthog::graph::xy_graph* 
        get_graph() { return g_; } 

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



    private:
        std::vector<std::vector<warthog::cpd::rle_run32>> fm_;
        std::vector<uint32_t> order_;
        warthog::graph::xy_graph* g_;
};

std::ostream&
operator<<(std::ostream& out, warthog::cpd::graph_oracle& o);

std::istream&
operator>>(std::istream& in, warthog::cpd::graph_oracle& o);

}

}

#endif

