#ifndef WARTHOG_FWD_CH_AF_EXPANSION_POLICY_H
#define WARTHOG_FWD_CH_AF_EXPANSION_POLICY_H

// contraction/fwd_ch_af_expansion_policy.h
//
// Forward-driven search in contraction hiearchies using 
// arc-flags to prune redundant up and down edges
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

class arcflags_filter;
class problem_instance;
class search_node;
class fwd_ch_af_expansion_policy : public expansion_policy
{
    public:
        fwd_ch_af_expansion_policy(
                warthog::graph::planar_graph* graph,
                std::vector<uint32_t>* rank, 
                warthog::arcflags_filter*);

        ~fwd_ch_af_expansion_policy();

		virtual void 
		expand(warthog::search_node*, warthog::problem_instance*);

        virtual void
        get_xy(warthog::search_node*, int32_t& x, int32_t& y);

        virtual uint32_t
        mem()
        {
            return expansion_policy::mem() +
                sizeof(this);
        }

    private:
        std::vector<uint32_t>* rank_;
        warthog::graph::planar_graph* g_;
        warthog::arcflags_filter* filter_;

        inline uint32_t
        get_rank(uint32_t id)
        {
            return rank_->at(id);
        }
};
}

#endif
