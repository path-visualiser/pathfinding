#ifndef WARTHOG_GRAPH_EXPANSION_POLICY_H
#define WARTHOG_GRAPH_EXPANSION_POLICY_H

#include "expansion_policy.h"

namespace warthog
{

class planar_graph;
class search_node;
class problem_instance;
class graph_expansion_policy : public expansion_policy
{

    public:
        graph_expansion_policy(warthog::planar_graph* g);
        virtual ~graph_expansion_policy();

		virtual void 
		expand(warthog::search_node*, warthog::problem_instance*);

        virtual void
        get_xy(warthog::search_node* n, int32_t& x, int32_t& y);

		uint32_t
		mem();

	private:
        warthog::planar_graph* g_;
};

}

#endif

