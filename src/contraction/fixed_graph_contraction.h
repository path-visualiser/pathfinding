#ifndef WARTHOG_FIXED_GRAPH_CONTRACTION_H
#define WARTHOG_FIXED_GRAPH_CONTRACTION_H

// fixed_graph_contraction.h
//
// a graph contraction scheme using a fixed 
// node ordering.
//
// @author: dharabor
// @created: 2016-06-15
//

#include "contraction.h"
#include "flexible_astar.h"

#include <vector>

using namespace warthog::ch;

namespace warthog
{

class apriori_filter;
class euclidean_heuristic;

template<class FILTER>
class graph_expansion_policy;

namespace graph
{

class planar_graph;

}

namespace ch
{

class fixed_graph_contraction 
{
    public:
        // default order is input order (i.e. crap)
        fixed_graph_contraction(
                warthog::graph::planar_graph* g);

        fixed_graph_contraction(
                warthog::graph::planar_graph* g,
                std::vector<uint32_t>* order);

        ~fixed_graph_contraction();

        void
        contract();

        void
        get_order(std::vector<uint32_t>& order)
        {
            order = *order_;
        }

        inline warthog::graph::planar_graph*
        get_graph() { return g_; }

        void
        set_partial_contraction_percentage(uint32_t pct)
        {
           c_pct_ =  std::min<uint32_t>(pct, 100);
        }
       
        // helper functions
        void
        set_verbose(bool verbose) { this->verbose_ = verbose; }


        size_t
        mem();

    protected:
        virtual void
        preliminaries();

        virtual void
        postliminaries() { } 

        virtual uint32_t
        next();

        virtual double
        witness_search(uint32_t from_id, uint32_t to_id, double via_len);

    private:
        void 
        init();

        uint32_t order_index_;
        std::vector<uint32_t>* order_;
        warthog::graph::planar_graph* g_;
        
        bool verbose_;
        bool done_;
        uint32_t c_pct_;

        // witness search stuff
        warthog::euclidean_heuristic* heuristic_;
        warthog::apriori_filter* filter_;
        warthog::graph_expansion_policy<warthog::apriori_filter>* expander_;
        warthog::flexible_astar<
           warthog::euclidean_heuristic,
           warthog::graph_expansion_policy<warthog::apriori_filter>,
           warthog::pqueue_min>* alg_;
        uint32_t total_expansions_;
        uint32_t total_searches_;
};

}
}

#endif
