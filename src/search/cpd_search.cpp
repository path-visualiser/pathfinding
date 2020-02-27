#include "cpd_search.h"

// TODO refactor node information inside Stats

template<class H, class E, class Q>
warthog::search_node*
warthog::cpd_search<H, E, Q>::search(warthog::solution& sol)
{
    warthog::timer mytimer;
    mytimer.start();
    open_->clear();

    warthog::search_node* start;
    warthog::search_node* incumbent = nullptr;

    // get the internal target id
    if(pi_.target_id_ != warthog::SN_ID_MAX)
    {
        warthog::search_node* target =
                expander_->generate_target_node(&pi_);
        if(!target) { return nullptr; } // invalid target location
        pi_.target_id_ = target->get_id();
    }

    // initialise and push the start node
    if(pi_.start_id_ == warthog::SN_ID_MAX) { return nullptr; }
    start = expander_->generate_start_node(&pi_);
    if(!start) { return nullptr; } // invalid start location
    pi_.start_id_ = start->get_id();

    // This heuristic also returns its upper bound
    warthog::cost_t start_h, start_ub;
    heuristic_->h(pi_.start_id_, pi_.target_id_, start_h, start_ub);
    upper_bound_ = start_ub;

    // `hscale` is contained in the heuristic
    start->init(pi_.instance_id_, warthog::SN_ID_MAX, 0, start_h, start_ub);

    open_->push(start);
    sol.nodes_inserted_++;

    if(on_generate_fn_)
    { (*on_generate_fn_)(start, 0, 0, UINT32_MAX); }

    info(pi_.verbose_, pi_);
    info(pi_.verbose_, "cut-off =", cost_cutoff_, "- tlim =", time_lim_,
         "- k-move =", max_k_move_);

    // begin expanding
    while(open_->size())
    {
        warthog::search_node* current = open_->pop();
        warthog::cost_t fval = current->get_f();

        current->set_expanded(true); // NB: set before generating
        assert(current->get_expanded());
        sol.nodes_expanded_++;

        if(on_expand_fn_) { (*on_expand_fn_)(current); }

        // goal test
        if(expander_->is_target(current, &pi_))
        {
            incumbent = current;
            break;
        }

        // early termination: in case we want bounded-cost
        // search or if we want to impose some memory limit
        if(fval > cost_cutoff_) { break; }
        if(sol.nodes_expanded_ >= exp_cutoff_) { break; }
        // Early stopping criterion: the heuristic has found a, possibly
        // pertubed, path to the destination which is acceptable.
        if(current->get_ub() <= fval) { break; }

        if (fval > upper_bound_)
        {
            debug(pi_.verbose_, "Pruning:", *current);

            continue;
        }
        else if (fval < upper_bound_)
        {
            debug(pi_.verbose_, "Updating UB to", fval);

            upper_bound_ = fval;
            incumbent = current;
        }

        trace(pi_.verbose_, sol.nodes_expanded_, "- Expanding:", *current);

        // generate successors
        expander_->expand(current, &pi_);
        warthog::cost_t cost_to_n = 0;
        uint32_t edge_id = 0;
        for(warthog::search_node* n = expander_->first(n, cost_to_n);
            n != nullptr;
            expander_->next(n, cost_to_n))
        {
            sol.nodes_touched_++;

            if(on_generate_fn_)
            { (*on_generate_fn_)(n, current, cost_to_n, edge_id++); }

            // add new nodes to the fringe
            if(n->get_search_number() != current->get_search_number())
            {
                warthog::cost_t gval = current->get_g() + cost_to_n;
                warthog::cost_t hval;
                warthog::cost_t ub;

                heuristic_->h(n->get_id(), pi_.target_id_, hval, ub);

                if (ub < warthog::COST_MAX)
                {
                    ub += gval;
                }

                // Should we check for overflow here?
                n->init(current->get_search_number(), current->get_id(),
                        gval, gval + hval, ub);

                open_->push(n);
                sol.nodes_inserted_++;

                debug(pi_.verbose_, "Generating:", *current);

                if(on_relax_fn_) { (*on_relax_fn_)(n); }
            }
            // The neighbour may have been expanded on the way to the target
            // following an inflated heuristic
            else if(n->get_expanded())
            {
                warthog::cost_t gval = current->get_g() + cost_to_n;
                warthog::cost_t ub = n->get_ub();

                // We need to check whether `n` actually improves the UB
                if (gval < n->get_g() && ub < warthog::COST_MAX &&
                    gval + ub < upper_bound_)
                {
                    open_->push(n);
                    debug(pi_.verbose_, "Reinsert:", *n);
                }
                else
                {
                    debug(pi_.verbose_, "Closed:", *n);
                }
            }
            // update a node from the fringe
            else if(open_->contains(n))
            {
                warthog::cost_t gval = current->get_g() + cost_to_n;
                if(gval < n->get_g())
                {
                    n->relax(gval, current->get_id());
                    open_->decrease_key(n);
                    sol.nodes_updated_++;

                    debug(pi_.verbose_, "Updating:", *n);

                    if(on_relax_fn_) { (*on_relax_fn_)(n); }
                }
                else
                {
                    debug(pi_.verbose_, "Not updating:", *n);
                }
            }
        }
    }

    mytimer.stop();
    sol.time_elapsed_nano_ = mytimer.elapsed_time_nano();

    DO_ON_DEBUG_IF(pi_.verbose_)
    {
        if(incumbent == nullptr)
        {
            warning(pi_.verbose_, "Search failed; no solution exists.");
        }
        else
        {
            user(pi_.verbose_, "Best incumbent", *incumbent);
        }
    }

    // Rebuild path from incumbent to solution
    while (!expander_->is_target(incumbent, &pi_))
    {
        warthog::sn_id_t parent = incumbent->get_id();
        warthog::sn_id_t n =
                heuristic_->next(parent, pi_.target_id_);

        if (n == warthog::SN_ID_MAX)
        {
            return nullptr;
        }

        incumbent = expander_->generate(n);
    }

    return incumbent;
}
