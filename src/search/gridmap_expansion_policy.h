#ifndef WARTHOG_GRIDMAP_EXPANSION_POLICY_H
#define WARTHOG_GRIDMAP_EXPANSION_POLICY_H

// gridmap_expansion_policy.h
//
// An ExpansionPolicy for square uniform-cost grids 
//
// @author: dharabor
// @created: 28/10/2010
//

#include "blocklist.h"
#include "expansion_policy.h"
#include "gridmap.h"
#include "search_node.h"

#include <memory>

namespace warthog
{

class problem_instance;
class gridmap_expansion_policy : public expansion_policy
{
	public:
		gridmap_expansion_policy(warthog::gridmap* map);
		virtual ~gridmap_expansion_policy() { }

		virtual void 
		expand(warthog::search_node*, warthog::problem_instance*);

        virtual void
        get_xy(warthog::search_node*, int32_t& x, int32_t& y);

		virtual uint32_t
		mem()
		{
            return 
                expansion_policy::mem() + 
                sizeof(*this) + 
                map_->mem();
		}
	
	private:
		warthog::gridmap* map_;
};

}

#endif

