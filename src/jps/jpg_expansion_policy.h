#ifndef JPG_EXPANSION_POLICY_H
#define JPG_EXPANSION_POLICY_H

// jpg_expansion_policy.h
//
// This expansion policy applies jps-style pruning to search a corner point
// graph created by the function warthog::jps::create_corner_graph
// (roughly, this is a jump point graph plus some intermediate vertices)
//
// @author: dharabor
// @created: 2016-09-21
//

#include "expansion_policy.h"
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

class jpg_expansion_policy : public expansion_policy
{

    public:
        jpg_expansion_policy(warthog::graph::corner_graph*);

        virtual 
        ~jpg_expansion_policy();

		virtual void 
		expand(warthog::search_node*, warthog::problem_instance*);

        virtual void
        get_xy(warthog::search_node* n, int32_t& x, int32_t& y);

        inline warthog::graph::corner_graph*
        get_graph() { return this->g_; }

        virtual size_t
        mem();

    private:
        // a graph of jump points
        warthog::graph::corner_graph* g_;
};

}

}

#endif

