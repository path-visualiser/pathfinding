#ifndef WARTHOG_GRAPH_EXPANSION_POLICY_H
#define WARTHOG_GRAPH_EXPANSION_POLICY_H

#include "constants.h"
#include "expansion_policy.h"

// search/graph_expansion_policy.cpp
//
// an expansion policy for planar graphs. includes support for a node 
// filtering mechanism (i.e. never generate nodes from a specified set)
//
// @author: dharabor
// @created: 2016-02-13
// 


namespace warthog
{

namespace graph
{

class planar_graph;

}

class problem_instance;
class search_node;

class graph_expansion_policy : public expansion_policy
{
    public:
        graph_expansion_policy(warthog::graph::planar_graph* g);
        virtual ~graph_expansion_policy();

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

	private:
        warthog::graph::planar_graph* g_;
};

}

#endif

