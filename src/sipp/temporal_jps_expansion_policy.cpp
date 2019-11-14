#include "sipp/temporal_jps_expansion_policy.h"

warthog::temporal_jps_expansion_policy::temporal_jps_expansion_policy(
    warthog::jpst_gridmap* jpst_gm) : jpst_gm_(jpst_gm)
{
    gm_map_width_ = jpst_gm_->gm_->width();
    map_width_ = jpst_gm_->gm_->header_width();
    sz_xy = jpst_gm_->gm_->header_height() * jpst_gm_->gm_->header_width();

    pool_.push_back(new warthog::mem::node_pool(sz_xy));
    neis_ = new warthog::arraylist<neighbour_record>(32);
    jpl_t_ = new warthog::jpst_locator(jpst_gm_);

    assert(__builtin_ffs(warthog::jps::NONE) == 0);
    assert(__builtin_ffs(warthog::jps::NORTH) == 1);
    assert(__builtin_ffs(warthog::jps::SOUTH) == 2);
    assert(__builtin_ffs(warthog::jps::EAST) == 3);
    assert(__builtin_ffs(warthog::jps::WEST) == 4);

    // incompatible moves that can lead to 
    // edge collisions
    ec_moves_[0] = warthog::cbs::NONE;
    ec_moves_[1] = warthog::cbs::SOUTH;
    ec_moves_[2] = warthog::cbs::NORTH;
    ec_moves_[3] = warthog::cbs::WEST;
    ec_moves_[4] = warthog::cbs::EAST;

    // for quickly computing successor xy_id values
    xy_id_offsets_[0] = 0;
    xy_id_offsets_[1] = -1 * (int32_t)map_width_;
    xy_id_offsets_[2] = (int32_t)map_width_;
    xy_id_offsets_[3] = 1;
    xy_id_offsets_[4] = -1;

    opposite_dir[0] = warthog::jps::NONE;
    opposite_dir[1] = warthog::jps::SOUTH;
    opposite_dir[2] = warthog::jps::NORTH;
    opposite_dir[3] = warthog::jps::WEST;
    opposite_dir[4] = warthog::jps::EAST;
}

warthog::temporal_jps_expansion_policy::~temporal_jps_expansion_policy()
{
    reset();

    delete jpl_t_;
    delete neis_;

    for(uint32_t i = 0; i < pool_.size(); i++)
    {
        delete pool_.at(i);
    }
}
