#ifndef WARTHOG_FCH_BBAF_EXPANSION_POLICY_H
#define WARTHOG_FCH_BBAF_EXPANSION_POLICY_H

// contraction/fch_bbaf_expansion_policy.h
//
// Forward-driven search in contraction hiearchies using 
// arc-flags + rectangular bounding boxes to prune redundant 
// up and down edges
//
// @author: dharabor
// @created: 2016-08-23
//

#include "expansion_policy.h"
#include <vector>

namespace warthog
{

namespace graph
{
class planar_graph;
}

class bbaf_filter;
class problem_instance;
class search_node;
class fch_bbaf_expansion_policy : public expansion_policy
{
    public:
        fch_bbaf_expansion_policy(
                warthog::graph::planar_graph* graph,
                std::vector<uint32_t>* rank, 
                warthog::bbaf_filter*);

        ~fch_bbaf_expansion_policy();

		virtual void 
		expand(warthog::search_node*, warthog::problem_instance*);

        virtual void
        get_xy(uint32_t node_id, int32_t& x, int32_t& y);

        virtual size_t
        mem()
        {
            return expansion_policy::mem() +
                sizeof(this);
        }

    private:
        std::vector<uint32_t>* rank_;
        warthog::graph::planar_graph* g_;
        warthog::bbaf_filter* filter_;
        uint32_t search_id_;

        inline uint32_t
        get_rank(uint32_t id)
        {
            return rank_->at(id);
        }
};
}

#endif
