#ifndef WARTHOG_LAZY_GRAPH_CONTRACTION_H
#define WARTHOG_LAZY_GRAPH_CONTRACTION_H

// lazy_graph_contraction.h
//
// Applies a contraction operation to each node in a 
// graph. 
//
// Contraction in this case refers to a localised graph
// operation where additional ``shortcut'' edges are added to 
// the graph in such a way that the ``contracted'' node can
// be bypassed entirely during search (i.e. its neighbours are
// connected directly).
//
// In a contraction hierarchy every node is assigned a specific
// ``level'' which is equal to the contraction priority of the
// node. The highest priority nodes (== lowest level) are 
// contracted first and the process continues until every node
// has been contracted.
//
// The node ordering is done in a lazy manner
// using a variety of heuristics.
//
// For more details see:
// [Geisbergerger, Sanders, Schultes and Delling. 
// Contraction Hierarchies: Faster and Simpler Hierarchical 
// Routing in Road Networks. In Proceedings of the 2008
// Workshop on Experimental Algorithms (WEA)]
//
// @author: dharabor
// @created: 2016-01-25
//

#include "constants.h"
#include "contraction.h"
#include "graph_contraction.h"
#include "heap.h"

#include <cstdint>
#include <unordered_map>

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

class ch_pair
{
    public:
        ch_pair() 
            : node_id_(0), cval_(warthog::INF) { } 

        ch_pair(uint32_t node_id, int32_t cval)
            : node_id_(node_id), cval_(cval) { }

        ch_pair(const ch_pair& other)
        { node_id_ = other.node_id_; cval_ = other.cval_; }

        ~ch_pair() { } 

        uint32_t node_id_; 
        int32_t cval_; // 'value' of contracting this node (lower is better)
};

bool
operator<(const ch_pair& first, const ch_pair& second);

class lazy_graph_contraction : public graph_contraction
{
    public:
        lazy_graph_contraction(warthog::graph::planar_graph* g);
        virtual ~lazy_graph_contraction();

        virtual size_t
        mem();

        // the order of contraction according to ::next
        void  
        get_order(std::vector<uint32_t>& order);

    protected:
        virtual void
        preliminaries();

        virtual void
        postliminaries();

        virtual uint32_t
        next();

        virtual double
        witness_search(uint32_t from_id, uint32_t to_id, double via_len);

    private:
        // node order stuff
        warthog::heap<ch_pair>* heap_;
        warthog::heap_node<ch_pair>* hn_pool_;
        std::vector<uint32_t> order_; // node ids, as retured by ::next

        // a variety of heuristics are used to order nodes
        struct node_order_terms
        {
            node_order_terms()
            {
                depth_ = 0; 
                nc_ = 0;
            }

            int32_t ed_; // edge difference 
            uint32_t depth_; // max search space depth
            uint32_t nc_; // neis contracted
        };
        node_order_terms* terms_;

        // witness search stuff
        uint32_t ws_max_expansions_; 
        warthog::euclidean_heuristic* heuristic_;
        warthog::graph_expansion_policy<warthog::apriori_filter>* expander_;
        warthog::apriori_filter* filter_;
        warthog::flexible_astar<
            warthog::euclidean_heuristic,
            warthog::graph_expansion_policy<warthog::apriori_filter>>* alg_;

        // metrics
        uint32_t total_expansions_;
        uint32_t total_searches_;
        uint32_t total_lazy_updates_;
        
        
        int32_t 
        edge_difference(uint32_t node_id, bool estimate);

        int32_t
        compute_contraction_value(uint32_t node_id);

};

}

}

#endif

