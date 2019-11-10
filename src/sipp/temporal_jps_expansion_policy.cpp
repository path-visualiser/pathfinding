#include "sipp/temporal_jps_expansion_policy.h"

warthog::temporal_jps_expansion_policy::temporal_jps_expansion_policy(
    warthog::sipp_gridmap* sipp_map) : sipp_map_(sipp_map)
{
    sz_xy = sipp_map_->gm_->header_height() * sipp_map_->gm_->header_width();
    pool_.push_back(new warthog::mem::node_pool(sz_xy));
    neis_ = new warthog::arraylist<neighbour_record>(32);

    gm_ = sipp_map_->gm_;
    t_gm_ = sipp_map_->t_gm_;
    map_width_ = gm_->header_width();
    gm_map_width_ = gm_->width();
}

warthog::temporal_jps_expansion_policy::~temporal_jps_expansion_policy()
{
    reset();

    delete jpl_temporal_;
    delete jpl_spatial_;
    delete t_gm_;
    delete neis_;

    for(uint32_t i = 0; i < pool_.size(); i++)
    {
        delete pool_.at(i);
    }
}
