#include "gridmap_time_expansion_policy.h"
#include "helpers.h"
#include "problem_instance.h"

warthog::gridmap_time_expansion_policy::gridmap_time_expansion_policy(
		warthog::gridmap* map, bool manhattan) : map_(map), manhattan_(manhattan)
{
    neis_ = new warthog::arraylist<neighbour_record>(32);

    uint32_t map_sz = map->height() * map->width();
    assert(map_sz > 0);

    time_map_ = new std::vector<warthog::node_pool*>();
    for(uint32_t i = 0; i < 256; i++)
    {
        time_map_->push_back(new warthog::node_pool(map_sz));
    }

    // setup some constants to quickly compute the current timestep 
    // and xy-index
    bitwidth_map_ = 32 - __builtin_clz(map->height()*map->width());
    id_mask_ = (1 << bitwidth_map_)-1;
}

warthog::gridmap_time_expansion_policy::~gridmap_time_expansion_policy()
{
    for(uint32_t i = 0; i < time_map_->size(); i++)
    {
        delete time_map_->at(i);
    }
    time_map_->clear();
    delete time_map_;
    delete neis_;
}


void 
warthog::gridmap_time_expansion_policy::expand(warthog::search_node* current,
		warthog::problem_instance* problem)
{
	reset();

	// get terrain type of each tile in the 3x3 square around (x, y)
	uint32_t tiles = 0;
	uint32_t nodeid = current->get_id() & id_mask_;
	map_->get_neighbours(nodeid, (uint8_t*)&tiles);

    // advance time
    uint32_t timestep = current->get_id() >> bitwidth_map_;
    timestep++;

	uint32_t nid_m_w = nodeid - map_->width();
	uint32_t nid_p_w = nodeid + map_->width();

    // cardinal successors
	if((tiles & 514) == 514) // N
	{  
		add_neighbour(__generate(nid_m_w, timestep), 1);
	} 
	if((tiles & 1536) == 1536) // E
	{
		add_neighbour(__generate(nodeid + 1, timestep), 1);
	}
	if((tiles & 131584) == 131584) // S
	{ 
		add_neighbour(__generate(nid_p_w, timestep), 1);
	}
	if((tiles & 768) == 768) // W
	{ 
		add_neighbour(__generate(nodeid - 1, timestep), 1);
	}

    if(manhattan_) { return; }

    // diagonal successors
	// NB: no corner cutting or squeezing between obstacles!
	if((tiles & 1542) == 1542) // NE
	{ 
        add_neighbour(__generate(nid_m_w + 1, timestep), warthog::DBL_ROOT_TWO);
	}

	
	if((tiles & 394752) == 394752) // SE
	{	
        add_neighbour(__generate(nid_p_w + 1, timestep), warthog::DBL_ROOT_TWO);
	}


	if((tiles & 197376) == 197376) // SW
	{ 
        add_neighbour(__generate(nid_p_w - 1, timestep), warthog::DBL_ROOT_TWO);
	}

	if((tiles & 771) == 771) // NW
	{ 
		add_neighbour(__generate(nid_m_w - 1, timestep), warthog::DBL_ROOT_TWO);
	}
}

void
warthog::gridmap_time_expansion_policy::get_xy(uint32_t nid, int32_t& x, int32_t& y)
{
    map_->to_unpadded_xy(nid & id_mask_, (uint32_t&)x, (uint32_t&)y);
}

warthog::search_node* 
warthog::gridmap_time_expansion_policy::generate_start_node(
        warthog::problem_instance* pi)
{ 
    uint32_t max_id = map_->header_width() * map_->header_height();
    if(pi->start_id_ >= max_id) { return 0; }
    uint32_t padded_id = map_->to_padded_id(pi->start_id_);
    return __generate(padded_id, 0);
}

warthog::search_node*
warthog::gridmap_time_expansion_policy::generate_target_node(
        warthog::problem_instance* pi)
{
    uint32_t max_id = map_->header_width() * map_->header_height();
    if(pi->target_id_ >= max_id) { return 0; }
    uint32_t padded_id = map_->to_padded_id(pi->target_id_);
    return __generate(padded_id, 0);
}

size_t
warthog::gridmap_time_expansion_policy::mem()
{
   size_t total = sizeof(*this) + map_->mem();
   for(uint32_t i = 0; i < time_map_->size(); i++)
   {
       total += time_map_->at(i)->mem();
   }
   total += sizeof(neighbour_record) * neis_->capacity();
   return total;
}
