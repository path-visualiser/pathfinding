#ifndef WARTHOG_CH_CPG_EXPANSION_POLICY_H
#define WARTHOG_CH_CPG_EXPANSION_POLICY_H

// contraction/ch_cpg_expansion_policy.h
//
// An expansion policy for contraction hierarchies on 
// top of corner point graphs.
//
// @author: dharabor
// @created: 2016-11-17
//

#include "contraction.h"
#include "expansion_policy.h"
#include "corner_point_graph.h"

#include <vector>

namespace warthog{

class problem_instance;
class search_node;

namespace ch
{

class ch_cpg_expansion_policy : public  expansion_policy
{
    public:
        // NB: @param rank: the contraction ordering used to create
        // the graph. this is a total order given in the form of 
        // an array with elements v[x] = i where x is the id of the 
        // node and i is its contraction rank
        //
        // @param backward: when true successors are generated by following 
        // incoming arcs rather than outgoing arcs (default is outgoing)
        //
        // @param sd: determines which implicit graph of the hierarchy 
        // is being searched: the up graph or the down graph.
        // in the up graph all successors have a rank larger than the
        // current node (this is the default). in the down graph, all 
        // successors have a rank smaller than the current node 
        ch_cpg_expansion_policy(warthog::graph::corner_point_graph* g, 
                std::vector<uint32_t>* rank, 
                bool backward=false,
                warthog::ch::search_direction sd = warthog::ch::UP);

        virtual 
        ~ch_cpg_expansion_policy() { }

		virtual void 
		expand(warthog::search_node*, warthog::problem_instance*);

        virtual void
        get_xy(uint32_t node_id, int32_t& x, int32_t& y)
        {
            g_->get_xy(node_id, x, y);
        }

        inline uint32_t
        get_rank(uint32_t id)
        {
            return rank_->at(id);
        }

        inline uint32_t
        get_num_nodes() { return g_->get_num_nodes(); }

        virtual warthog::search_node* 
        generate_start_node(warthog::problem_instance* pi)
        {
            if(pi->get_search_id() != search_id_at_last_insert_)
            {
                g_->insert(pi->get_start_id(), pi->get_target_id());
                search_id_at_last_insert_ = pi->get_search_id();
            }
            return this->generate(g_->get_inserted_start_id());
        }

        virtual warthog::search_node*
        generate_target_node(warthog::problem_instance* pi)
        {
            if(pi->get_search_id() != search_id_at_last_insert_)
            {
                g_->insert(pi->get_start_id(), pi->get_target_id());
                search_id_at_last_insert_ = pi->get_search_id();
            }
            return this->generate(g_->get_inserted_target_id());
        }

        virtual size_t
        mem();

    private:
        uint32_t search_id_at_last_insert_;
        bool backward_;
        warthog::graph::corner_point_graph* g_;
        std::vector<uint32_t>* rank_;
        warthog::ch::search_direction sd_;


};

}

}
#endif

