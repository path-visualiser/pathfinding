#ifndef WARTHOG_CH_EXPANSION_POLICY_H
#define WARTHOG_CH_EXPANSION_POLICY_H

// contraction/ch_expansion_policy.h
//
// an expansion policy for contraction hierarchies. during
// search the only successors generated are those which
// are lexically larger than the current node. 
//
// For more details see:
// [Geisbergerger, Sanders, Schultes and Delling. 
// Contraction Hierarchies: Faster and Simpler Hierarchical 
// Routing in Road Networks. In Proceedings of the 2008
// Workshop on Experimental Algorithms (WEA)]
//
// @author: dharabor
// @created: 2016-05-10
//

#include "expansion_policy.h"
#include "planar_graph.h"

#include <vector>

namespace warthog{

class problem_instance;
class search_node;

class ch_expansion_policy : public  expansion_policy
{
    public:
        ch_expansion_policy(warthog::graph::planar_graph* g, 
                std::vector<uint32_t>* order)
            : expansion_policy(g->get_num_nodes()), g_(g), order_(order) { }

        virtual 
        ~ch_expansion_policy() { }

		virtual void 
		expand(warthog::search_node*, warthog::problem_instance*);

        virtual void
        get_xy(warthog::search_node* n, int32_t& x, int32_t& y)
        {
            g_->get_xy(n->get_id(), x, y);
        }

        uint32_t
        mem();

    private:
        warthog::graph::planar_graph* g_;
        std::vector<uint32_t>* order_;
};

}
#endif

