#include "sipp/jpst_gridmap.h"
#include <algorithm>

warthog::jpst_gridmap::jpst_gridmap(warthog::gridmap* gm)
{
    gm_ = gm;
    sipp_map_ = new warthog::sipp_gridmap(gm);

    // one copy of the map for jumping E<->W; one copy for jumping N<->S
    t_gm_ = new warthog::gridmap(gm_->header_height(), gm_->header_width());
    t_gm_r_ = new warthog::gridmap(gm_->header_width(), gm_->header_height());
    
    // we want to jump over locations which have no temporal obstacles and
    // stop at locations that have none or which are simply blocked
    // we use a 1 to indicate temporal obstacles and locations beyond
    // the edge of the map. we use a 0 to indicate other locations
    t_gm_->invert(); 
    t_gm_r_->invert(); 

    // in the beginning every grid tile (including obstacles) is associated
    // with zero temporal obstacles. so we put a 0 for these locations
    for(uint32_t i = 0; i < t_gm_->header_width() * t_gm_->header_height(); i++)
    { 
        uint32_t gm_id = t_gm_->to_padded_id(i);
        t_gm_->set_label(gm_id, 0); 
        t_gm_r_->set_label(this->map_id_to_rmap_id(gm_id), 0);
    }
}

warthog::jpst_gridmap::~jpst_gridmap()
{
    delete t_gm_r_;
    delete t_gm_;
    delete sipp_map_;
}

// create a copy of the grid map which is rotated by 90 degrees clockwise.
// this version will be used when jumping North or South. 
warthog::gridmap*
warthog::jpst_gridmap::create_rmap()
{
	uint32_t maph = t_gm_->header_height();
	uint32_t mapw = t_gm_->header_width();
	uint32_t rmaph = mapw;
	uint32_t rmapw = maph;
	warthog::gridmap* rmap = new warthog::gridmap(rmaph, rmapw);

	for(uint32_t x = 0; x < mapw; x++) 
	{
		for(uint32_t y = 0; y < maph; y++)
		{
			uint32_t label = t_gm_->get_label(t_gm_->to_padded_id(x, y));
			uint32_t rx = ((rmapw-1) - y);
			uint32_t ry = x;
			uint32_t rid = rmap->to_padded_id(rx, ry);
			rmap->set_label(rid, label);
		}
	}
	return rmap;
}
