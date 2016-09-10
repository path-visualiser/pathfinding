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
#include "graph_contraction.h"

#include <vector>

using namespace warthog::ch;

namespace warthog
{

class apriori_filter;
class euclidean_heuristic;
class graph_expansion_policy;

namespace graph
{

class planar_graph;

}

namespace ch
{

class fixed_graph_contraction : public warthog::ch::graph_contraction
{
    public:
        // default order is input order (i.e. crap)
        fixed_graph_contraction(
                warthog::graph::planar_graph* g);

        fixed_graph_contraction(
                warthog::graph::planar_graph* g,
                std::vector<uint32_t>* order);

        virtual ~fixed_graph_contraction();

        void
        get_order(std::vector<uint32_t>& order)
        {
            order = *order_;
        }

        virtual size_t
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

        // witness search stuff
        warthog::euclidean_heuristic* heuristic_;
        warthog::graph_expansion_policy* expander_;
        warthog::apriori_filter* filter_;
        warthog::flexible_astar<
           warthog::euclidean_heuristic,
           warthog::graph_expansion_policy,
           warthog::apriori_filter>* alg_;
        uint32_t total_expansions_;
        uint32_t total_searches_;
};

}
}

#endif
