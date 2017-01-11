#ifndef WARTHOG_WGRIDMAP_EXPANSION_POLICY_H
#define WARTHOG_WGRIDMAP_EXPANSION_POLICY_H

// wgridmap_expansion_policy.h
//
// An expansion policy for square gridmaps with
// non-uniform costs.
//
// @author: dharabor
// @created: 2014-09-17
//

#include "expansion_policy.h"
#include "search_node.h"
#include "weighted_gridmap.h"

#include <memory>

namespace warthog
{

class wgridmap_expansion_policy : public expansion_policy
{
	public:
		wgridmap_expansion_policy(warthog::weighted_gridmap* map);
		virtual ~wgridmap_expansion_policy();

		virtual void 
		expand(warthog::search_node*, warthog::problem_instance*);

		virtual size_t
		mem()
		{
			return 
                expansion_policy::mem() + 
                sizeof(*this) + map_->mem();
		}

        virtual void
        get_xy(uint32_t node_id, int32_t& x, int32_t& y);

        virtual warthog::search_node* 
        generate_start_node(warthog::problem_instance* pi);

        virtual warthog::search_node*
        generate_target_node(warthog::problem_instance* pi);


	private:
		warthog::weighted_gridmap* map_;
};

}

#endif

