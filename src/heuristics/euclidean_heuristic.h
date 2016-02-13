#ifndef WARTHOG_EUCLIDEAN_HEURISTIC_H
#define WARTHOG_EUCLIDEAN_HEURISTIC_H

// euclidean_heuristic.h
//
// Straight-line heuristic for measuring distances in the plane.
//
// @author: dharabor
// @created: 2016-02-11
//
//

#include "constants.h"
#include "planar_graph.h"

namespace warthog
{

class euclidean_heuristic
{
    public:
        euclidean_heuristic(warthog::planar_graph* g) : g_(g) { }
        ~euclidean_heuristic() { }

		inline warthog::cost_t
		h(int32_t x, int32_t y, 
				int32_t x2, int32_t y2)
		{
            // NB: precision loss when warthog::cost_t is an integer
			double dx = x-x2;
			double dy = y-y2;
            return sqrt(dx*dx + dy*dy);
		}

		inline warthog::cost_t
		h(uint32_t id, uint32_t id2)
		{
			int32_t x, x2;
			int32_t y, y2;
            g_->get_xy(id, x, y);
            g_->get_xy(id2, x2, y2);
			return this->h(x, y, x2, y2);
		}

	private:
        warthog::planar_graph* g_;

};

}

#endif

