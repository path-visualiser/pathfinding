#include "gridmap_expansion_policy.h"
#include "helpers.h"
#include "problem_instance.h"

warthog::gridmap_expansion_policy::gridmap_expansion_policy(
		warthog::gridmap* map)
: expansion_policy(map->height()*map->width()), map_(map)
{
}

void 
warthog::gridmap_expansion_policy::expand(warthog::search_node* current,
		warthog::problem_instance* problem)
{
	reset();

	// get terrain type of each tile in the 3x3 square around (x, y)
	uint32_t tiles = 0;
	uint32_t nodeid = current->get_id();
	map_->get_neighbours(nodeid, (uint8_t*)&tiles);

//	#ifndef NDEBUG
//	uint32_t cx_, cy_;
//	warthog::helpers::index_to_xy(nodeid, map_->width(), cx_, cy_);
//	assert(tiles[0] == map_->get_label(cx_-1, cy_-1));
//	assert(tiles[1] == map_->get_label(cx_, cy_-1));
//	assert(tiles[2] == map_->get_label(cx_+1, cy_-1));
//	assert(tiles[3] == map_->get_label(cx_-1, cy_));
//	assert(tiles[4] == map_->get_label(cx_, cy_));
//	assert(tiles[5] == map_->get_label(cx_+1, cy_));
//	assert(tiles[6] == map_->get_label(cx_-1, cy_+1));
//	assert(tiles[7] == map_->get_label(cx_, cy_+1));
//	assert(tiles[8] == map_->get_label(cx_+1, cy_+1));
//	#endif

	// NB: no corner cutting or squeezing between obstacles!
	uint32_t nid_m_w = nodeid - map_->width();
	uint32_t nid_p_w = nodeid + map_->width();

	if((tiles & 514) == 514) // N
	{  
		add_neighbour(this->generate(nid_m_w), warthog::ONE);
	} 

	if((tiles & 1542) == 1542) // NE
	{ 
        add_neighbour(this->generate(nid_m_w + 1), warthog::ROOT_TWO);
	}

	if((tiles & 1536) == 1536) // E
	{
		add_neighbour(this->generate(nodeid + 1), warthog::ONE);
	}
	
	if((tiles & 394752) == 394752) // SE
	{	
        add_neighbour(this->generate(nid_p_w + 1), warthog::ROOT_TWO);
	}

	if((tiles & 131584) == 131584) // S
	{ 
		add_neighbour(this->generate(nid_p_w), warthog::ONE);
	}

	if((tiles & 197376) == 197376) // SW
	{ 
        add_neighbour(this->generate(nid_p_w - 1), warthog::ROOT_TWO);
	}

	if((tiles & 768) == 768) // W
	{ 
		add_neighbour(this->generate(nodeid - 1), warthog::ONE);
	}

	if((tiles & 771) == 771) // NW
	{ 
		add_neighbour(this->generate(nid_m_w - 1), warthog::ROOT_TWO);
	}
}

void
warthog::gridmap_expansion_policy::get_xy(warthog::search_node* n,
        int32_t& x, int32_t& y)
{
    
    map_->to_unpadded_xy(n->get_id(), (uint32_t&)x, (uint32_t&)y);
}

