#ifndef WARTHOG_CHAF_EXPANSION_POLICY_H
#define WARTHOG_CHAF_EXPANSION_POLICY_H

// contraction/chaf_expansion_policy.h
//
// Bi-directional contraction hierarchies + arc flags
//
// @author: dharabor
// @created: 2016-05-10
//

#include "contraction.h"
#include "expansion_policy.h"
#include "planar_graph.h"

#include <vector>

namespace warthog{

class af_filter;
class problem_instance;
class search_node;

class chaf_expansion_policy : public  expansion_policy
{
    public:
        // @param backward: when true successors are generated by following 
        // incoming arcs rather than outgoing arcs (default is outgoing)
        //
        // @param filter: the arc-flags filter used to prune arcs
        //
        chaf_expansion_policy(warthog::graph::planar_graph* g, 
                warthog::af_filter* filter,
                bool backward=false);

        virtual 
        ~chaf_expansion_policy() { }

		virtual void 
		expand(warthog::search_node*, warthog::problem_instance*);

        virtual void
        get_xy(uint32_t node_id, int32_t& x, int32_t& y);

        virtual warthog::search_node* 
        generate_start_node(warthog::problem_instance* pi);

        virtual warthog::search_node*
        generate_target_node(warthog::problem_instance* pi);

        virtual size_t
        mem();

        inline uint32_t
        get_num_nodes() { return g_->get_num_nodes(); }

    private:
        bool backward_;
        warthog::graph::planar_graph* g_;
        warthog::af_filter* filter_;

        typedef warthog::graph::edge_iter
                (warthog::chaf_expansion_policy::*chep_get_iter_fn) 
                (warthog::graph::node* n);

        // pointers to the neighbours in the direction of the search
        chep_get_iter_fn fn_begin_iter_;
        chep_get_iter_fn fn_end_iter_;

        // pointers to neighbours in the reverse direction to the search
        chep_get_iter_fn fn_rev_begin_iter_;
        chep_get_iter_fn fn_rev_end_iter_;

        inline warthog::graph::edge_iter
        get_fwd_begin_iter(warthog::graph::node* n) 
        { return n->outgoing_begin(); }

        inline warthog::graph::edge_iter
        get_fwd_end_iter(warthog::graph::node* n) 
        { return n->outgoing_end(); }

        inline warthog::graph::edge_iter
        get_bwd_begin_iter(warthog::graph::node* n) 
        { return n->incoming_begin(); }

        inline warthog::graph::edge_iter
        get_bwd_end_iter(warthog::graph::node* n) 
        { return n->incoming_end(); }
};

}
#endif

