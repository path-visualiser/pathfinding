#ifndef WARTHOG_FCH_BB_JPG_EXPANSION_POLICY_H
#define WARTHOG_FCH_BB_JPG_EXPANSION_POLICY_H

// contraction/fch_bb_jpg_expansion_policy.h
//
// Forward-driven with contraction hiearchies applied to 
// corner point graphs. 
//  + jump point pruning.
//  + geometric containers.
//
// @author: dharabor
// @created: 2017-01-25
//

#include "expansion_policy.h"
#include "geom.h"
#include <vector>

namespace warthog
{

namespace graph
{
class corner_point_graph;
}

class problem_instance;
class search_node;

class bb_filter;

class fch_bb_jpg_expansion_policy : public expansion_policy
{
    public:
        fch_bb_jpg_expansion_policy(
                warthog::graph::corner_point_graph* graph,
                std::vector<uint32_t>* rank,
                warthog::bb_filter* nf);

        ~fch_bb_jpg_expansion_policy();

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
            // TODO: implement me
            return 0; 
        }

    private:
        uint32_t search_id_at_last_insert_;
        std::vector<uint32_t>* rank_;
        warthog::graph::corner_point_graph* g_;

        // geometric containers stuff
        warthog::bb_filter* nf_;
        warthog::geom::rectangle r_;
        bool
        filter(uint32_t node_id, uint32_t edge_id);

        // *** JPS STUFF ***
        
        // we label each edge of the CH with two pieces of information:
        // (1) a FIRST label that specifies the direction of the first
        // step on the diagonal-first path that begins at the edge's 
        // tail node and finishes at its head node. 
        // (2) a LAST label that is similar and which specifies the 
        // direction of the last step on the path
        enum step_type
        {
            FIRST = 0,
            LAST = 1
        };

        inline uint32_t
        get_rank(uint32_t id)
        {
            return rank_->at(id);
        }

        void
        init();

        void
        compute_direction_labels();

        void
        label_edge(
                warthog::graph::edge* e, uint32_t e_tail_id, 
                warthog::graph::planar_graph* pg);

        void
        process_edge(warthog::graph::edge* e, uint32_t e_tail_id,
                     warthog::graph::planar_graph* g);

        inline warthog::jps::direction
        get_dir(warthog::graph::edge* e, step_type which)
        {
            return (warthog::jps::direction)
                ((uint8_t*)(&e->label_))[which && which];
        }

        inline void
        set_dir(warthog::graph::edge* e, step_type which, 
                warthog::jps::direction d)
        {
            ((uint8_t*)(&e->label_))[which && which] = d;
        }
};
}

#endif
