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
                std::vector<uint32_t>* ooc, bool backward=false)
            : expansion_policy(g->get_num_nodes()), g_(g)
        {
            backward_ = backward;

            // we convert the order of contraction into a ranked 
            // list s.t. rank[i] gives the contraction index of
            // node i
            assert(g->get_num_nodes() - g->id_offset() == ooc->size());

            // because DIMACS graphs are stupid and 1-indexed the 
            // planar_graph adds padding to the front of the array.
            // we need to account for the padding when assigning
            // ranks to indexes
            rank_.resize(g->get_num_nodes());
            for(uint32_t i = 0; i < g->id_offset(); i++)
            {
                rank_.at(i) = 0;
            }

            for(uint32_t i = 0; i < ooc->size(); i++)
            {
                rank_.at(ooc->at(i)) = i + g->id_offset();
            }
        }

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
        bool backward_;
        warthog::graph::planar_graph* g_;
        std::vector<uint32_t> rank_;
};

}
#endif

