#include "helpers.h"
#include "problem_instance.h"
#include "wgridmap_expansion_policy.h"

warthog::wgridmap_expansion_policy::wgridmap_expansion_policy(
		warthog::weighted_gridmap* map) 
    : expansion_policy(map->height() * map->width()), map_(map)
{
}

warthog::wgridmap_expansion_policy::~wgridmap_expansion_policy()
{
}

void 
warthog::wgridmap_expansion_policy::expand(warthog::search_node* current,
		warthog::problem_instance* problem)
{
    reset();

	// get terrain type of each tile in the 3x3 square around (x, y)
    uint32_t tile_ids[9];
    warthog::dbword tiles[9];
	uint32_t nodeid = current->get_id();
	map_->get_neighbours(nodeid, tile_ids, tiles);

    // NB: 
    // 1. neighbours generated in clockwise order starting from direction N.
    // 2. transition costs to each neighbour are calculated as the average of 
    // terrain values for all tiles touched by the agent when moving.
    // CAVEAT EMPTOR: precision loss from integer maths.
    if(tiles[1]) // N neighbour
    {
        warthog::search_node* n = generate(tile_ids[1]);
        warthog::cost_t cost = ((tiles[1] + tiles[4]) * warthog::ONE) >> 1; 
        add_neighbour(n, cost);
    }
    if(tiles[1] & tiles[2] & tiles[5]) // NE neighbour
    {
        warthog::search_node* n =  generate(tile_ids[2]);
        warthog::cost_t cost = ((tiles[1] + tiles[2] + tiles[4] + tiles[5]) * 
                warthog::ROOT_TWO) >> 2; 
        add_neighbour(n, cost);
    }
    if(tiles[5]) // E
    {
        warthog::search_node* n = generate(tile_ids[5]);
		warthog::cost_t cost = ((tiles[5] + tiles[4]) * warthog::ONE) >> 1; 
        add_neighbour(n, cost);
    }
    if(tiles[5] & tiles[7] & tiles[8]) // SE
    {
        warthog::search_node* n = generate(tile_ids[8]);
        warthog::cost_t cost = ((tiles[4] + tiles[5] + tiles[7] + tiles[8]) * 
                warthog::ROOT_TWO) >> 2; 
        add_neighbour(n, cost);
    }
    if(tiles[7]) // S
    {
        warthog::search_node* n = generate(tile_ids[7]);
        warthog::cost_t cost = ((tiles[7] + tiles[4]) * warthog::ONE) >> 1; 
        add_neighbour(n, cost);
    }
    if(tiles[3] & tiles[6] & tiles[7]) // SW
    {
        warthog::search_node* n = generate(tile_ids[6]);
        warthog::cost_t cost =  ((tiles[3] + tiles[4] + tiles[6] + tiles[7]) * 
                warthog::ROOT_TWO) >> 2; 
        add_neighbour(n, cost);
    }
    if(tiles[3]) // W
    {
        warthog::search_node* n = generate(tile_ids[3]);
        warthog::cost_t cost = ((tiles[3] + tiles[4]) * warthog::ONE) >> 1; 
        add_neighbour(n, cost);
    }
    if(tiles[0] & tiles[1] & tiles[3]) // NW neighbour
    {
        warthog::search_node* n = generate(tile_ids[0]);
        warthog::cost_t cost = ((tiles[0] + tiles[1] + tiles[3] + tiles[4]) * 
                warthog::ROOT_TWO) >> 2; 
        add_neighbour(n, cost);
    }
}

void
warthog::wgridmap_expansion_policy::get_xy(warthog::search_node* n,
        int32_t& x, int32_t& y)
{
    
    map_->to_unpadded_xy(n->get_id(), (uint32_t&)x, (uint32_t&)y);
}
