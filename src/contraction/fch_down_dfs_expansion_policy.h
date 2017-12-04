#ifndef WARTHOG_FCH_DOWN_DFS_EXPANSION_POLICY_H
#define WARTHOG_FCH_DOWN_DFS_EXPANSION_POLICY_H

// contraction/fch_down_dfs_expansion_policy.h
//
// An expansion policy that combines FCH with a "down pruning" scheme.
// The idea is to label each down edge with a node-id range such that if
// the target is in the id-range it means the edge appears on an optimal
// down-path to the target.
//
// @author: dharabor
// @created: 2017-12-02
//

#include "planar_graph.h"
#include "expansion_policy.h"
#include <vector>

namespace warthog
{

namespace graph
{
class planar_graph;
}

class problem_instance;
class search_node;
class fch_down_dfs_expansion_policy : public expansion_policy
{
    public:
        fch_down_dfs_expansion_policy(
                warthog::graph::planar_graph* graph,
                std::vector<uint32_t>* rank);

        ~fch_down_dfs_expansion_policy();

		virtual void 
		expand(warthog::search_node*, warthog::problem_instance*);

        virtual void
        get_xy(uint32_t node_id, int32_t& x, int32_t& y);

        virtual warthog::search_node* 
        generate_start_node(warthog::problem_instance* pi);

        virtual warthog::search_node*
        generate_target_node(warthog::problem_instance* pi);

        virtual inline size_t
        mem()
        {
            size_t retval = sizeof(this);
            for(uint32_t i = 0; i < g_->get_num_nodes(); i++)
            {
                retval += 
                    sizeof(fch_interval)*
                    edge_labels_->at(i).size();
            }
            retval += sizeof(int32_t) * node_labels_->size();
            retval += expansion_policy::mem();
            return retval;
        }

    private:
        std::vector<uint32_t>* rank_;
        warthog::graph::planar_graph* g_;
        uint8_t* down_heads_;

        struct fch_interval
        {
            fch_interval() { left = INT32_MAX; right = INT32_MIN; } 

            fch_interval&
            operator=(const fch_interval& other)
            {
                left = other.left; right = other.right; 
                return *this;
            }

            void
            grow(int32_t val)
            {
                left = val < left  ? val : left;
                right = val > right ? val : right;
            }

            void
            merge(fch_interval& other)
            {
                left = left < other.left ? left : other.left;
                right = right > other.right ? right : other.right;
            }

            bool
            contains(int32_t val)
            {
                return (right - val) + left >= left;
            }

            int32_t left, right;
        };

        std::vector< std::vector< fch_interval >>* edge_labels_;
        std::vector<int32_t>* node_labels_;

        uint32_t s_label, t_label;

        inline bool
        filter(uint32_t node_idx, uint32_t edge_idx)
        {
            assert(edge_idx < edge_labels_->at(node_idx).size());
            return !edge_labels_->at(node_idx).at(edge_idx).contains(t_label);
        }

        inline uint32_t
        get_rank(uint32_t id)
        {
            return rank_->at(id);
        }

        // TODO: stick this in a labelling class
        void 
        compute_down_dijkstra_postorder();
};
}

#endif
