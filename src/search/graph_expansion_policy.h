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

class search_node;
class problem_instance;
class graph_expansion_policy : public expansion_policy
{

    public:
        graph_expansion_policy(warthog::graph::planar_graph* g);
        virtual ~graph_expansion_policy();

		virtual void 
		expand(warthog::search_node*, warthog::problem_instance*);

        // never generate filtered nodes
        void 
        filter(uint32_t node_id);

        void
        unfilter(uint32_t node_id);

        void 
        reset_filter();

        virtual void
        get_xy(warthog::search_node* n, int32_t& x, int32_t& y);

		uint32_t
		mem();

	private:
        warthog::graph::planar_graph* g_;
        warthog::dbword* filter_;
        uint32_t filter_sz_;
};

}

#endif

