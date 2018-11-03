#include "grid.h"
#include "helpers.h"
#include "manhattan_time_expansion_policy.h"
#include "problem_instance.h"

#include <algorithm>

warthog::manhattan_time_expansion_policy::manhattan_time_expansion_policy(
		warthog::gridmap* map, bool manhattan) : map_(map)
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

warthog::manhattan_time_expansion_policy::~manhattan_time_expansion_policy()
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
warthog::manhattan_time_expansion_policy::expand(warthog::search_node* current,
		warthog::problem_instance* problem)
{
	reset();

    // get the xy id of the current node and extract current timestep
	uint32_t nodeid = current->get_id() & id_mask_;
    uint32_t timestep = current->get_id() >> bitwidth_map_;

	// get adjacent grid tiles (bitpacked into one 32bit word)
	uint32_t tiles = 0;
	map_->get_neighbours(nodeid, (uint8_t*)&tiles);

    // neighbour ids are calculated using nodeid offsets
	uint32_t nid_m_w = nodeid - map_->width();
	uint32_t nid_p_w = nodeid + map_->width();

    auto get_constraints = [problem] (uint32_t padded_id, uint32_t timestep)
        -> cell_constraints
    {
        warthog::time_constraints* cons = 
            (warthog::time_constraints*)problem->extra_params_;
        if(cons)
        {
            return cons->get_constraints(padded_id, timestep);
        }
        return cell_constraints();
    };

    // edge constraints for the current node
    cell_constraints cur_cc = get_constraints(nodeid, timestep);

    // cardinal successors
    cell_constraints succ_cc = get_constraints(nid_m_w, timestep+1);
    if( ((tiles & 514) == 514) && // N
        !(cur_cc.e_ & warthog::grid::NORTH) &&  
        !succ_cc.v_ ) 
	{  
		add_neighbour(__generate(nid_m_w, timestep+1), 1);
	} 

    succ_cc = get_constraints(nodeid + 1, timestep+1);
	if( ((tiles & 1536) == 1536) && // E
        !(cur_cc.e_ & warthog::grid::EAST) && 
        !succ_cc.v_ )
	{
		add_neighbour(__generate(nodeid + 1, timestep+1), 1);
	}

    succ_cc = get_constraints(nid_p_w, timestep+1);
	if( ((tiles & 131584) == 131584) && // S
        !(cur_cc.e_ & warthog::grid::SOUTH) && 
        !succ_cc.v_ )
	{ 

		add_neighbour(__generate(nid_p_w, timestep+1), 1);
	}

    succ_cc = get_constraints(nodeid - 1, timestep+1);
	if( ((tiles & 768) == 768) && // W
        !(cur_cc.e_ & warthog::grid::WEST) && 
        !succ_cc.v_ )
	{ 
		add_neighbour(__generate(nodeid - 1, timestep+1), 1);
	}

    // wait successor
    succ_cc = get_constraints(nodeid, timestep+1);
    if(!succ_cc.v_)
    {
        add_neighbour(__generate(nodeid, timestep+1), 1);
    }
}

void
warthog::manhattan_time_expansion_policy::get_xy(uint32_t nid, int32_t& x, int32_t& y)
{
    map_->to_unpadded_xy(nid & id_mask_, (uint32_t&)x, (uint32_t&)y);
}

warthog::search_node* 
warthog::manhattan_time_expansion_policy::generate_start_node(
        warthog::problem_instance* pi)
{ 
    uint32_t max_id = map_->header_width() * map_->header_height();
    if(pi->start_id_ >= max_id) { return 0; }
    uint32_t padded_id = map_->to_padded_id(pi->start_id_);
    return __generate(padded_id, 0);
}

warthog::search_node*
warthog::manhattan_time_expansion_policy::generate_target_node(
        warthog::problem_instance* pi)
{
    uint32_t max_id = map_->header_width() * map_->header_height();
    if(pi->target_id_ >= max_id) { return 0; }
    uint32_t padded_id = map_->to_padded_id(pi->target_id_);
    return __generate(padded_id, 0);
}

size_t
warthog::manhattan_time_expansion_policy::mem()
{
   size_t total = sizeof(*this) + map_->mem();
   for(uint32_t i = 0; i < time_map_->size(); i++)
   {
       total += time_map_->at(i)->mem();
   }
   total += sizeof(neighbour_record) * neis_->capacity();
   return total;
}
