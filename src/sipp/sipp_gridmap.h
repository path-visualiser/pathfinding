#ifndef WARTHOG_SIPP_GRIDMAP_H
#define WARTHOG_SIPP_GRIDMAP_H

// domains/sipp_gridmap.h
// 
// A time-expanded gridmap where the temporal dimension is represented
// using "safe intervals" as described in the folllwing paper:
//
// [Phillips, M. and Likhachev, M., 2011, 
// SIPP: Safe interval path planning for Dynamic Environments. 
// In 2011 IEEE International Conference on Robotics and Automation 
// (pp. 5628-5635)]
//
// This implementation adds to SIPP support for moving obstacles.
//
// @author: dharabor
// @created: 2019-10-25
//

#include "domains/gridmap.h"
#include "mapf/cbs.h"

#include <vector>

namespace warthog
{

namespace sipp
{

struct safe_interval
{
    safe_interval()
    {
        s_time_ = 0;
        e_time_ = warthog::COST_MAX;
        action_ = warthog::cbs::move::WAIT;
        reach_time_ = COST_MAX;
    }

    // start and end times for the safe interval
    cost_t s_time_;
    cost_t e_time_;

    // the earliest time this interval has been reached
    // during the most recent search episode
    cost_t reach_time_; 

    // a safe interval is created when an obstacle disappears.
    // when the obstacle is due to some moving agent, we need
    // to track in which direction the agent moved when it
    // stepped away at time s_time_.
    warthog::cbs::move action_;
};

}

// id is padded_id (4 bytes) and then an index for the interval (
class sipp_gridmap
{
    public:
        sipp_gridmap(warthog::gridmap* gm);

        // add a temporal obstacle at location (@param x, @param y). 
        // after the call, the location is blocked for the duration of the
        // **OPEN** interval (@param start_time, @param end_time).
        // 
        // NB: runs in time _linear_ to the number of intervals at (x, y)
        void
        add_obstacle(uint32_t x, uint32_t y, cost_t start_time, cost_t end_time, 
                     warthog::cbs::move action = warthog::cbs::move::WAIT);

        // remove all temporal obstacles at location (@param x, @param y)
        // after the call this location has a single safe interval. 
        // for traversable tile, this interval begins
        // at time 0 and finishes at time INF.
        // for obstacle tiles, this interval begins and ends at time
        // COST_MAX (i.e. infinity)
        void
        clear_obstacles(uint32_t x, uint32_t y);

        // return a reference to the (@param index)th safe interval 
        // associated with the grid cell index @param node_id. 
        warthog::sipp::safe_interval&
        get_safe_interval(uint32_t node_id, uint32_t index)
        {
            assert(index < intervals_.at(node_id).size());
            return intervals_.at(node_id).at(index);
        }

        // return a reference to the list of safe intervals for the
        // grid cell with index @param node_id
        std::vector<warthog::sipp::safe_interval>& 
        get_all_intervals(uint32_t node_id)
        { return intervals_.at(node_id); }

        size_t
        mem()
        {
            size_t retval = 0;
            for(uint32_t i = 0; i < intervals_.size(); i++)
            {
                retval += sizeof(warthog::sipp::safe_interval)*intervals_.at(i).size();
                retval += sizeof(intervals_.at(i));
            }
            retval += gm_->mem();
            retval += sizeof(this);
            return retval;
        }


        warthog::gridmap* gm_;

    private:
        std::vector<std::vector<warthog::sipp::safe_interval>> intervals_;
};

}

#endif
