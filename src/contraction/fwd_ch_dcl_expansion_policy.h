#ifndef WARTHOG_FWD_CH_DCL_EXPANSION_POLICY
#define WARTHOG_FWD_CH_DCL_EXPANSION_POLICY

// contraction/fwd_ch_dcl_expansion_policy.h
//
// An expansion policy for forward-driven
// search in contraction hiearchies combined 
// with a bounding-box filtering scheme
//
// @author: dharabor
// @created: 2016-08-02
//

#include "expansion_policy.h"
#include <vector>

namespace warthog
{

namespace graph
{
class planar_graph;
}

class problem_instance;
class search_node;
class euclidean_heuristic;
class dcl_filter;

class fwd_ch_dcl_expansion_policy : public expansion_policy
{
    public:
        fwd_ch_dcl_expansion_policy(
                warthog::graph::planar_graph* graph,
                std::vector<uint32_t>* rank, 
                warthog::dcl_filter* nf);

        ~fwd_ch_dcl_expansion_policy();

		virtual void 
		expand(warthog::search_node*, warthog::problem_instance*);

        virtual void
        get_xy(warthog::search_node*, int32_t& x, int32_t& y);

        void
        set_apex(uint32_t apex) 
        { 
            apex_ = apex; 
            apex_reached_ = (apex == warthog::INF) ? true : false; 
        }

        virtual size_t
        mem()
        {
            return expansion_policy::mem() +
                sizeof(this);
        }

    private:
        std::vector<uint32_t>* rank_;
        warthog::graph::planar_graph* g_;
        warthog::dcl_filter* nf_;
        uint32_t apex_;
        bool apex_reached_;

        inline uint32_t
        get_rank(uint32_t id)
        {
            return rank_->at(id);
        }
};
}

#endif
