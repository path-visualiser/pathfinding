#ifndef WARTHOG_GRAPH_EXPANSION_POLICY_H
#define WARTHOG_GRAPH_EXPANSION_POLICY_H

#include "constants.h"
#include "expansion_policy.h"

// search/graph_expansion_policy.cpp
//
// an expansion policy for planar graphs. includes support for a node 
// filtering mechanism (i.e. it can be configured to prune successor
// nodes that do not match some specified criteria)
//
// @author: dharabor
// @created: 2016-02-13
// 


#include "dummy_filter.h"
#include "planar_graph.h"
#include "problem_instance.h"
#include "search_node.h"

namespace warthog
{

template <class FILTER = warthog::dummy_filter>
class graph_expansion_policy : public expansion_policy
{
    public:
        graph_expansion_policy(
                warthog::graph::planar_graph* g, FILTER* filter = 0)
            : expansion_policy(g->get_num_nodes()), g_(g), filter_(filter)
        {
            assert(g);
            if(filter == 0)
            {
                fn_generate_successor = &warthog::graph_expansion_policy
                                        <FILTER>::fn_generate_no_filter;
            }
            else
            {
                fn_generate_successor = &warthog::graph_expansion_policy
                                        <FILTER>::fn_generate_with_filter;
            }
        }

        virtual ~graph_expansion_policy() { }

		virtual void 
		expand(warthog::search_node* current, warthog::problem_instance* pi)
        {
            reset();

            uint32_t current_id = current->get_id();
            warthog::graph::node* n = g_->get_node(current_id) ;
            warthog::graph::edge_iter begin = n->outgoing_begin();
            warthog::graph::edge_iter end = n->outgoing_end();
            
            for(warthog::graph::edge_iter it = begin; it != end; it++)
            {
                warthog::graph::edge& e = *it;
                assert(e.node_id_ < g_->get_num_nodes());
                (this->*fn_generate_successor)(e.node_id_, it - begin, e);
            }
        }

        virtual void
        get_xy(uint32_t node_id, int32_t& x, int32_t& y)
        {
            g_->get_xy(node_id, x, y);
        }

        virtual warthog::search_node* 
        generate_start_node(warthog::problem_instance* pi)
        {
            uint32_t s_graph_id = g_->to_graph_id(pi->start_id_);
            if(s_graph_id == warthog::INF) { return 0; }
            return generate(s_graph_id);
        }

        virtual warthog::search_node*
        generate_target_node(warthog::problem_instance* pi)
        {
            uint32_t t_graph_id = g_->to_graph_id(pi->target_id_);
            if(t_graph_id == warthog::INF) { return 0; }

            // also update the filter with the new target location
            filter_->set_target(t_graph_id);

            return generate(t_graph_id);
        }

        virtual size_t
		mem()
        {
            return 
                expansion_policy::mem() + 
                g_->mem() +
                sizeof(this);
        }

	private:
        warthog::graph::planar_graph* g_;
        FILTER* filter_;

        typedef void(warthog::graph_expansion_policy<FILTER>::*generate_fn)
            (uint32_t current_id, uint32_t edge_idx, warthog::graph::edge& e);

        generate_fn fn_generate_successor;

        inline void 
        fn_generate_with_filter(
                uint32_t current_id, 
                uint32_t edge_idx, 
                warthog::graph::edge& e)
        {
            if(!filter_->filter(current_id, edge_idx))
            {
                this->add_neighbour(this->generate(e.node_id_), e.wt_);
            }
        }

        inline void
        fn_generate_no_filter(
                uint32_t current_id, 
                uint32_t edge_idx, 
                warthog::graph::edge& e)
        {
            this->add_neighbour(this->generate(e.node_id_), e.wt_);
        }

};

typedef 
warthog::graph_expansion_policy<> simple_graph_expansion_policy;

}

#endif

