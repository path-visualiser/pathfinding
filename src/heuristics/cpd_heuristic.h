#ifndef WARTHOG_CPD_HEURISTIC_H
#define WARTHOG_CPD_HEURISTIC_H

// cpd_heuristic.h
//
// @author: amaheo
// @created: 27/02/2020
//

#include "constants.h"
#include "helpers.h"
#include "xy_graph.h"

namespace warthog
{

// TODO Replace with actual CPD
class cpd_heuristic : public euclidean_heuristic
{
  public:
    cpd_heuristic(warthog::graph::xy_graph* g)
            : euclidean_heuristic(g) {}

    void
    h(warthog::sn_id_t id,
      warthog::sn_id_t id2,
      warthog::cost_t &h,
      warthog::cost_t &ub);

    warthog::sn_id_t
    next(warthog::sn_id_t current, warthog::sn_id_t target);
};

}

#endif
