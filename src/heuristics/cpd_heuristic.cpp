#include "cpd_heuristic.h"

void
warthog::cpd_heuristic::h(warthog::sn_id_t current, warthog::sn_id_t target,
                          warthog::cost_t &h, warthog::cost_t &ub)
{
    h = warthog::euclidean_heuristic::h(current, target);
    ub = warthog::COST_MAX;
}

// Should not be used for now
warthog::sn_id_t
warthog::cpd_heuristic::next(warthog::sn_id_t current, warthog::sn_id_t target)
{
    return warthog::SN_ID_MAX;
}
