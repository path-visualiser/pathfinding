#include "gridmap.h"
#include "jps.h"
#include "sipp/jpst_locator.h"

#include <cassert>
#include <climits>

warthog::jpst_locator::jpst_locator(warthog::jpst_gridmap* jpst_map)
	: jpst_gm_(jpst_map)
{ }

warthog::jpst_locator::~jpst_locator()
{ }


// Finds a jump point successor of node (x, y) in Direction d.
// Also given is the location of the goal node (goalx, goaly) for a particular
// search instance. If encountered, the goal node is always returned as a 
// jump point successor.
//
// @return: the id of a jump point successor or warthog::INF32 if no jp exists.
void
warthog::jpst_locator::jump(warthog::jps::direction d,
	   	uint32_t node_id, uint32_t goal_id, uint32_t& jumpnode_id, 
		double& jumpcost)
{
	switch(d)
	{
		case warthog::jps::NORTH:
			jump_north(node_id, goal_id, jumpnode_id, jumpcost);
			break;
		case warthog::jps::SOUTH:
			jump_south(node_id, goal_id, jumpnode_id, jumpcost);
			break;
		case warthog::jps::EAST:
			jump_east(node_id, goal_id, jumpnode_id, jumpcost);
			break;
		case warthog::jps::WEST:
			jump_west(node_id, goal_id, jumpnode_id, jumpcost);
			break;
		default:
			break;
	}
}

void
warthog::jpst_locator::jump_north(uint32_t node_id, 
		uint32_t goal_id, uint32_t& jumpnode_id, double& jumpcost)
{
    uint32_t num_steps = 0;
    uint32_t mapw = jpst_gm_->t_gm_->width();

    uint32_t jp_w_id;
    uint32_t jp_e_id;
    double jp_w_cost;
    double jp_e_cost;

    uint32_t x, y; 
    jpst_gm_->t_gm_->to_unpadded_xy(node_id, x, y);
    uint32_t max_steps_west = x;
    uint32_t max_steps_east = (jpst_gm_->t_gm_->header_width()-1) - x;

    uint32_t next_id = node_id;
    while(true)
    {
        __jump_east(next_id, goal_id, jp_e_id, jp_e_cost, jpst_gm_->t_gm_);
        __jump_west(next_id, goal_id, jp_w_id, jp_w_cost, jpst_gm_->t_gm_);
        if( (jp_e_cost <= max_steps_east) || (jp_w_cost <= max_steps_west) ) { break; }

        next_id -= mapw;
        num_steps++;
        if(!jpst_gm_->gm_->get_label(next_id)) 
        { 
            num_steps--;
            break;
        }
    }
    jumpnode_id = node_id - (mapw*num_steps);
    jumpcost = num_steps;
}

void
warthog::jpst_locator::__jump_north(uint32_t node_id, 
		uint32_t goal_id, uint32_t& jumpnode_id, double& jumpcost,
		warthog::gridmap* mymap)
{
	// jumping north in the original map is the same as jumping
	// east when we use a version of the map rotated 90 degrees.
	__jump_east(node_id, goal_id, jumpnode_id, jumpcost, jpst_gm_->t_gm_r_);
}

void
warthog::jpst_locator::jump_south(uint32_t node_id, 
		uint32_t goal_id, uint32_t& jumpnode_id, double& jumpcost)
{
    uint32_t num_steps = 0;
    uint32_t mapw = jpst_gm_->t_gm_->width();

    uint32_t jp_w_id;
    uint32_t jp_e_id;
    double jp_w_cost;
    double jp_e_cost;

    uint32_t x, y; 
    jpst_gm_->t_gm_->to_unpadded_xy(node_id, x, y);
    uint32_t max_steps_west = x;
    uint32_t max_steps_east = (jpst_gm_->t_gm_->header_width()-1) - x;

    uint32_t next_id = node_id;
    while(true)
    {
        __jump_east(next_id, goal_id, jp_e_id, jp_e_cost, jpst_gm_->t_gm_);
        __jump_west(next_id, goal_id, jp_w_id, jp_w_cost, jpst_gm_->t_gm_);
        if( (jp_e_cost <= max_steps_east) || (jp_w_cost <= max_steps_west) ) { break; }

        next_id += mapw;
        num_steps++;
        if(!jpst_gm_->gm_->get_label(next_id)) 
        { 
            num_steps--;
            break;
        }
    }
    jumpnode_id = node_id + (mapw*num_steps);
    jumpcost = num_steps;
}

void
warthog::jpst_locator::__jump_south(uint32_t node_id, 
		uint32_t goal_id, uint32_t& jumpnode_id, double& jumpcost,
		warthog::gridmap* mymap)
{
	// jumping north in the original map is the same as jumping
	// west when we use a version of the map rotated 90 degrees.
	__jump_west(node_id, goal_id, jumpnode_id, jumpcost, jpst_gm_->t_gm_r_);
}

void
warthog::jpst_locator::jump_east(uint32_t node_id, 
		uint32_t goal_id, uint32_t& jumpnode_id, double& jumpcost)
{
	__jump_east(node_id, goal_id, jumpnode_id, jumpcost, jpst_gm_->t_gm_);
}


void
warthog::jpst_locator::__jump_east(uint32_t node_id, 
		uint32_t goal_id, uint32_t& jumpnode_id, double& jumpcost, 
		warthog::gridmap* mymap)
{
	jumpnode_id = node_id;

	uint32_t neis[3] = {0, 0, 0};

	jumpnode_id = node_id;
	while(true)
	{
		// read in tiles from 3 adjacent rows. the curent node 
		// is in the low byte of the middle row
		mymap->get_neighbours_32bit(jumpnode_id, neis);

		// stop if we found any forced or dead-end tiles
		uint32_t stop_bits = neis[0] | neis[1] | neis[2];
		if(stop_bits)
		{
			uint32_t stop_pos = (uint32_t)__builtin_ffs((int)stop_bits)-1; // returns idx+1
			jumpnode_id += (uint32_t)stop_pos; 
			break;
		}
		jumpnode_id += 32;
	}
	jumpcost = jumpnode_id - node_id;
}

// analogous to ::jump_east 
void
warthog::jpst_locator::jump_west(uint32_t node_id, 
		uint32_t goal_id, uint32_t& jumpnode_id, double& jumpcost)
{
	__jump_west(node_id, goal_id, jumpnode_id, jumpcost, jpst_gm_->t_gm_);
}

void
warthog::jpst_locator::__jump_west(uint32_t node_id, 
		uint32_t goal_id, uint32_t& jumpnode_id, double& jumpcost, 
		warthog::gridmap* mymap)
{
	uint32_t neis[3] = {0, 0, 0};

	jumpnode_id = node_id;
	while(true)
	{
		// cache 32 tiles from three adjacent rows.
		// current tile is in the high byte of the middle row
		mymap->get_neighbours_upper_32bit(jumpnode_id, neis);

		uint32_t stop_bits = neis[0] | neis[1] | neis[2];
		if(stop_bits)
		{
			uint32_t stop_pos = (uint32_t)__builtin_clz(stop_bits);
			jumpnode_id -= stop_pos;
			break;
		}
		jumpnode_id -= 32;
	}
	jumpcost = node_id - jumpnode_id;
}

