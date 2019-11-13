#include "sipp/temporal_jps_expansion_policy.h"

warthog::temporal_jps_expansion_policy::temporal_jps_expansion_policy(
    warthog::jpst_gridmap* jpst_gm) : jpst_gm_(jpst_gm)
{
    sz_xy = jpst_gm_->gm_->header_height() * jpst_gm_->gm_->header_width();
    pool_.push_back(new warthog::mem::node_pool(sz_xy));
    neis_ = new warthog::arraylist<neighbour_record>(32);
    

    gm_ = jpst_gm_->gm_;
    t_gm_ = jpst_gm_->t_gm_;

    map_width_ = gm_->header_width();
    gm_map_width_ = gm_->width();

    jpl_spatial_ = new warthog::online_jump_point_locator(gm_);
    jpl_temporal_ = new warthog::jpst_locator(jpst_gm_);
}

warthog::temporal_jps_expansion_policy::~temporal_jps_expansion_policy()
{
    reset();

    delete jpl_temporal_;
    delete jpl_spatial_;
    delete neis_;

    for(uint32_t i = 0; i < pool_.size(); i++)
    {
        delete pool_.at(i);
    }
}
