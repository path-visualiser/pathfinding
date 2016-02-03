#include "jps2plus_expansion_policy.h"

warthog::jps2plus_expansion_policy::jps2plus_expansion_policy(warthog::gridmap* map)
    : expansion_policy(map->height() * map->width())
{
	map_ = map;
	jpl_ = new warthog::offline_jump_point_locator2(map);

	costs_.reserve(100);
	jp_ids_.reserve(100);
    reset();
}

warthog::jps2plus_expansion_policy::~jps2plus_expansion_policy()
{
	delete jpl_;
}

void 
warthog::jps2plus_expansion_policy::expand(
		warthog::search_node* current, warthog::problem_instance* problem)
{
	reset();
    costs_.clear();
    jp_ids_.clear();

	// compute the direction of travel used to reach the current node.
	warthog::jps::direction dir_c = current->get_pdir();

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
			jpl_->jump(d, current_id, goal_id, jp_ids_, costs_);
		}
	}

	uint32_t searchid = problem->get_searchid();
	uint32_t id_mask = (1 << 24)-1;
	for(uint32_t i = 0; i < jp_ids_.size(); i++)
	{
		// bits 0-23 store the id of the jump point
		// bits 24-31 store the direction to the parent
		uint32_t jp_id = jp_ids_.at(i);
		warthog::jps::direction pdir = (warthog::jps::direction)*(((uint8_t*)(&jp_id))+3);

		warthog::search_node* mynode = generate(jp_id & id_mask);
		if(mynode->get_searchid() != searchid) { mynode->reset(searchid); }
        add_neighbour(mynode, costs_.at(i));

		// stupid hack
		if((current->get_g() + costs_.at(i)) < mynode->get_g())
		{
			mynode->set_pdir(pdir);
		}
	}
}

void
warthog::jps2plus_expansion_policy::get_xy(
        warthog::search_node* n, int32_t& x, int32_t& y)
{
    map_->to_unpadded_xy(n->get_id(), (uint32_t&)x, (uint32_t&)y);
}

