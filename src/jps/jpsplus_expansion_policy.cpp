#include "jpsplus_expansion_policy.h"

warthog::jpsplus_expansion_policy::jpsplus_expansion_policy(warthog::gridmap* map)
    : expansion_policy(map->width() * map->height())
{
	map_ = map;
	jpl_ = new warthog::offline_jump_point_locator(map);
}

warthog::jpsplus_expansion_policy::~jpsplus_expansion_policy()
{
	delete jpl_;
}

void 
warthog::jpsplus_expansion_policy::expand(
		warthog::search_node* current, warthog::problem_instance* problem)
{
    reset();

	// compute the direction of travel used to reach the current node.
	warthog::jps::direction dir_c =
	   	this->compute_direction(current->get_parent(), current);

	// get the tiles around the current node c
	uint32_t c_tiles;
	uint32_t current_id = current->get_id();
	map_->get_neighbours(current_id, (uint8_t*)&c_tiles);

	// look for jump points in the direction of each natural 
	// and forced neighbour
	uint32_t succ_dirs = warthog::jps::compute_successors(dir_c, c_tiles);
	uint32_t goal_id = problem->get_goal();
	for(uint32_t i = 0; i < 8; i++)
	{
		warthog::jps::direction d = (warthog::jps::direction) (1 << i);
		if(succ_dirs & d)
		{
			warthog::cost_t jumpcost;
			uint32_t succ_id;
			jpl_->jump(d, current_id, goal_id, succ_id, jumpcost);

			if(succ_id != warthog::INF)
			{
				add_neighbour(this->generate(succ_id), jumpcost);
			}
		}
	}
}

