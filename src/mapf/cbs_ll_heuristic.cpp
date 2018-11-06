#include "cbs_ll_heuristic.h"
#include "flexible_astar.h"
#include "gridmap_expansion_policy.h"
#include "zero_heuristic.h"

#include <functional>
#include <stdint.h>


warthog::cbs_ll_heuristic::cbs_ll_heuristic(warthog::gridmap* gm)
    : gm_(gm)
{
    t_index_ = 0;
    uint32_t bitwidth_map = 
        32 - __builtin_clz(gm->width() * gm->height());
    id_mask_ = (1 << bitwidth_map)-1;
}

void
warthog::cbs_ll_heuristic::compute_h_values(
        std::vector<uint32_t>& targets)
{
    t_map_.clear();
    h_.clear();
    h_.resize(targets.size());

    warthog::pqueue_min open;
    warthog::zero_heuristic h;
    warthog::gridmap_expansion_policy expander(gm_, true);

    warthog::flexible_astar<
        warthog::zero_heuristic, 
        warthog::gridmap_expansion_policy,
        warthog::pqueue_min> 
            alg(&h, &expander, &open);

    std::function<void(warthog::search_node*)> on_expand_fn = 
        [&] (warthog::search_node* current) -> void
        {
            uint32_t t_index = t_map_.size() - 1;
            uint32_t h_sz = h_.size();
            uint32_t h_sz_2 = h_.at(t_index).size();
            h_sz = h_sz + 1 - 1;
            h_sz_2 = h_sz_2 + 1 - 1;
            h_[t_index][current->get_id()] = current->get_g();
        };
    alg.apply_on_expand(on_expand_fn);

    for(uint32_t target_id : targets)
    {
        uint32_t t_index = t_map_.size();
        t_map_[target_id] = t_index;
        h_[t_index_].resize(gm_->width() * gm_->height(), warthog::INF);
        warthog::problem_instance problem(target_id, warthog::INF);
        warthog::solution sol;
        alg.get_distance(problem, sol);
    }
}
