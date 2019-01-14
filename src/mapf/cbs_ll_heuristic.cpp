#include "cbs_ll_heuristic.h"
#include "flexible_astar.h"
#include "gridmap.h"
#include "gridmap_expansion_policy.h"
#include "labelled_gridmap.h"
#include "vl_gridmap_expansion_policy.h"
#include "zero_heuristic.h"

#include <functional>
#include <stdint.h>

void
warthog::cbs_ll_heuristic::compute_h_values(
        std::vector<uint32_t>& targets,
        warthog::gridmap* map)
{
    t_map_.clear();
    h_.clear();
    h_.resize(targets.size());
    t_index_ = 0;

    warthog::pqueue_min open;
    warthog::zero_heuristic h;

    bool manhattan = true;
    warthog::gridmap_expansion_policy expander(map, manhattan);

    warthog::flexible_astar<
        warthog::zero_heuristic, 
        warthog::gridmap_expansion_policy,
        warthog::pqueue_min> 
            alg(&h, &expander, &open);

    uint32_t tmp_idx;
    std::function<void(warthog::search_node*)> on_expand_fn = 
        [&] (warthog::search_node* current) -> void
        {
            h_[tmp_idx][(uint32_t)current->get_id()] = current->get_g();
        };
    alg.apply_on_expand(on_expand_fn);

    for(uint32_t target_id : targets)
    {
        // skip targets already processed
        if(t_map_.find(target_id) != t_map_.end()) { continue; }

        tmp_idx = (uint32_t)t_map_.size();
        t_map_[target_id] = tmp_idx;

        h_[tmp_idx].resize(map->width() * map->height(), warthog::INF32);
        warthog::problem_instance problem(target_id);
        warthog::solution sol;
        alg.get_distance(problem, sol);
    }
}

