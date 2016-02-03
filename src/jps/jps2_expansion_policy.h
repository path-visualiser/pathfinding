#ifndef WARTHOG_JPS2_EXPANSION_POLICY_H
#define WARTHOG_JPS2_EXPANSION_POLICY_H

// jps2_expansion_policy.h
//
// An experimental variation of warthog::jps2_expansion_policy,
// this version works with a modified version of online JPS
// which skips intermediate jump points (i.e. those jps 
// that do not have any forced neighbours)
//
// @author: dharabor
// @created: 06/01/2010

#include "expansion_policy.h"
#include "gridmap.h"
#include "helpers.h"
#include "jps.h"
#include "online_jump_point_locator2.h"
#include "problem_instance.h"
#include "search_node.h"

#include "stdint.h"

namespace warthog
{

class jps2_expansion_policy : public expansion_policy
{
	public:
		jps2_expansion_policy(warthog::gridmap* map);
		virtual ~jps2_expansion_policy();

		virtual void 
		expand(warthog::search_node*, warthog::problem_instance*);

		virtual inline uint32_t
		mem()
		{
			return expansion_policy::mem() + 
                sizeof(*this) + map_->mem() + jpl_->mem();
		}

        virtual void
        get_xy(warthog::search_node* n, int32_t& x, int32_t& y); 

	private:
		warthog::gridmap* map_;
		online_jump_point_locator2* jpl_;
		std::vector<uint32_t> jp_ids_;
        std::vector<warthog::cost_t> jp_costs_;
};

}

#endif

