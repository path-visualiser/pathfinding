#ifndef CORNER_GRAPH_EXPANSION_POLICY_H
#define CORNER_GRAPH_EXPANSION_POLICY_H

// corner_graph_expansion_policy.h
//
// This expansion policy creates a graph of corner points
// and searches for paths in that graph. Problems where
// the start location or the target location (or both)
// are not corner points are handled by an insertion 
// procedure.
//
// @author: dharabor
// @created: 2016-09-21
//

#include "expansion_policy.h"
#include "gridmap.h"
#include "planar_graph.h"
#include <unordered_map>

namespace warthog
{
class search_node;

namespace graph
{
class corner_graph;
}

namespace jps
{

class corner_point_locator;
class corner_graph_expansion_policy : public expansion_policy
{
    public:
        corner_graph_expansion_policy(warthog::graph::corner_graph*);

        virtual 
        ~corner_graph_expansion_policy();

		virtual void 
		expand(warthog::search_node*, warthog::problem_instance*);

        virtual void
        get_xy(warthog::search_node* n, int32_t& x, int32_t& y);

        inline warthog::graph::corner_graph*
        get_graph() { return this->g_; }

        virtual size_t
        mem();

    private:
        warthog::graph::corner_graph* g_;
};

}

}

#endif

