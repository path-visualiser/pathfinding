#ifndef JPG_EXPANSION_POLICY_H
#define JPG_EXPANSION_POLICY_H

// jpg_expansion_policy.h
//
// This expansion policy applies jps-style pruning to search
// a graph of corner points.
// Problems where the start location or the target location (or both)
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

namespace jps
{

class corner_point_locator;
class jpg_expansion_policy : public expansion_policy
{
    typedef std::unordered_map<uint32_t, uint32_t>::iterator id_map_iter;

    public:
        jpg_expansion_policy(
                warthog::graph::planar_graph*, 
                warthog::gridmap* gm);

        virtual 
        ~jpg_expansion_policy();

		virtual void 
		expand(warthog::search_node*, warthog::problem_instance*);

        virtual void
        get_xy(warthog::search_node* n, int32_t& x, int32_t& y);

        inline warthog::graph::planar_graph*
        get_graph() { return this->g_; }

        // insert into ::g_ start and target locations that
        // are not jump points
        void
        insert(uint32_t start_grid_id, uint32_t target_grid_id,
               uint32_t& start_graph_id, uint32_t& target_graph_id);

        // inverse of insert operations
        void
        uninsert();

        virtual size_t
        mem();

    private:
        // a graph of jump points
        warthog::graph::planar_graph* g_;

        // we use this to insert into ::g_ start and 
        // target locations that are not jump points
        warthog::jps::corner_point_locator* cpl_;

        warthog::problem_instance* problem_;

        uint32_t s_graph_id_, s_grid_id_;
        int32_t start_x_, start_y_;
        uint32_t t_graph_id_, t_grid_id_;
        std::vector<warthog::graph::edge> s_outgoing_;
        std::vector<uint32_t> t_incoming_;
        std::vector<warthog::graph::edge_iter> t_incoming_iters_;

        // maps gridmap ids to graph ids
        std::unordered_map<uint32_t, uint32_t> id_map_;

        // keep track of which nodes we are inserting
        std::vector<uint32_t> inserted_ids_;

        uint32_t
        insert_start(uint32_t start_grid_id, uint32_t target_grid_id);
        uint32_t
        insert_target(uint32_t start_grid_id, uint32_t target_grid_id);
};

}

}

#endif

