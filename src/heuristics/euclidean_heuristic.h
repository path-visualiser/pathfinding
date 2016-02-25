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

typedef void (*xyFn)(uint32_t id, int32_t& x, int32_t& y);
class euclidean_heuristic
{
    public:
        euclidean_heuristic(warthog::planar_graph* g) : g_(g) { hscale_ = 1; }
        ~euclidean_heuristic() { }

		inline warthog::cost_t
		h(int32_t x, int32_t y, 
				int32_t x2, int32_t y2)
		{
            // NB: precision loss when warthog::cost_t is an integer
			double dx = x-x2;
			double dy = y-y2;
            return sqrt(dx*dx + dy*dy) * hscale_;
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

        inline void
        set_hscale(double hscale)
        {
            if(hscale > 0)
            {
                hscale_ = hscale;
            }
        }

        inline double
        get_hscale() { return hscale_; }



	private:
        warthog::planar_graph* g_;
        double hscale_;

};

}

#endif

