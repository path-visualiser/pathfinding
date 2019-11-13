#ifndef WARTHOG_JPST_LOCATOR_H
#define WARTHOG_JPST_LOCATOR_H

// sipp/jpst_locator.h
//
// A wrapper for Temporal JPS grid scanning operations.
// During each such scan the objective is to identify 
// grid cells which contain temporal obstacles.
// 
// This implementation focuses on 4 connected grids.
//
// @author: dharabor
// @created: 2019-11-09
//

#include "jps.h"
#include "sipp/jpst_gridmap.h"
#include "gridmap.h"

namespace warthog
{

class jpst_locator 
{
	public: 
		jpst_locator(warthog::jpst_gridmap* jpst_map);
		~jpst_locator();

		void
		jump(warthog::jps::direction d, uint32_t node_id, uint32_t goalid, 
				uint32_t& jumpnode_id, double& jumpcost);

		inline size_t 
		mem()
		{
			return sizeof(this);
		}

	private:
		void
		jump_north(uint32_t node_id, uint32_t goal_id, 
				uint32_t& jumpnode_id, double& jumpcost);
		void
		jump_south(uint32_t node_id, uint32_t goal_id, 
				uint32_t& jumpnode_id, double& jumpcost);
		void
		jump_east(uint32_t node_id, uint32_t goal_id, 
				uint32_t& jumpnode_id, double& jumpcost);
		void
		jump_west(uint32_t node_id, uint32_t goal_id, 
				uint32_t& jumpnode_id, double& jumpcost);

		// these versions can be passed a map parameter to
		// use when jumping. they allow switching between
		// map_ and rmap_ (a rotated counterpart).
		void
		__jump_east(uint32_t node_id, uint32_t goal_id, 
				uint32_t& jumpnode_id, double& jumpcost, 
				warthog::gridmap* mymap);
		void
		__jump_west(uint32_t node_id, uint32_t goal_id, 
				uint32_t& jumpnode_id, double& jumpcost, 
				warthog::gridmap* mymap);
		void
		__jump_north(uint32_t node_id, uint32_t goal_id, 
				uint32_t& jumpnode_id, double& jumpcost,
				warthog::gridmap* mymap);
		void
		__jump_south(uint32_t node_id, uint32_t goal_id, 
				uint32_t& jumpnode_id, double& jumpcost, 
				warthog::gridmap* mymap);

        warthog::jpst_gridmap* jpst_gm_;
};

}

#endif

