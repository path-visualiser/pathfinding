#ifndef WARTHOG_GRAPH_CORNER_GRAPH_H
#define WARTHOG_GRAPH_CORNER_GRAPH_H

// domains/corner_graph.h
//
// A corner_graph is wrapper around a planar graph created from a grid map.
// Typically the nodes of the wrapped graph are grid locations having 
// two adjacent neighbours that are not co-visible. We refer to all such 
// points as corner points.
//
// The edges of the wrapped graph can be created in different ways. For 
// example, connecting all pairs of co-visible corner points yields a
// visibility graph. Other possibilities include jump point graphs
// and subgoal graph.
//
// Common to all of the above is an insertion operation which is used 
// to handle cases when the start or target location is not a corner point.
// 
// @author: dharabor
// @created: 2016-11-06
//

#include "corner_point_locator.h"
#include "gridmap.h"
#include "planar_graph.h"

#include <unordered_map>

namespace warthog
{

namespace jps
{

class corner_point_locator;

}

namespace graph
{

class corner_graph
{
    public:
        corner_graph( warthog::gridmap* map );
                

        virtual ~corner_graph();

        inline void
        print_dimacs_gr(std::ostream& oss)
        {
            g_->print_dimacs_gr(oss);
        }

        inline void
        print_dimacs_co(std::ostream& oss)
        {
            g_->print_dimacs_co(oss);
        }

        // returns number of nodes + ID_OFFSET number of
        // dummy elements located at the start of the nodes array
        inline uint32_t
        get_num_nodes()
        {
            return g_->get_num_nodes();
        }

        inline uint32_t 
        get_num_edges()
        {
            return g_->get_num_edges();
        }

        void
        get_xy(uint32_t graph_id, int32_t& x, int32_t& y);

        inline warthog::graph::node* 
        get_node(uint32_t id)
        {
            return g_->get_node(id);
        }

        inline uint32_t
        add_node(int32_t x, int32_t y)
        {
            return g_->add_node(x, y);
        }

        inline void 
        set_verbose(bool verbose) 
        { 
            verbose_ = verbose; 
            g_->set_verbose(verbose);
        } 

        inline bool
        get_verbose() { return verbose_; }

        inline const char* 
        get_filename() { return g_->get_filename(); }

        inline warthog::gridmap* 
        get_gridmap() { return cpl_->get_map(); } 

        inline bool
        reserve(uint32_t new_cap) { return g_->reserve(new_cap); }

        virtual size_t
        mem();

        // insert the start and target locations into the graph.
        void
        insert( uint32_t start_grid_id, uint32_t target_grid_id);

        uint32_t 
        get_inserted_start_id() { return s_graph_id_; }

        uint32_t 
        get_inserted_target_id() { return t_graph_id_; }

        warthog::graph::planar_graph*
        get_planar_graph() { return g_; }

    private:
        bool verbose_;
        warthog::graph::planar_graph* g_;
        warthog::jps::corner_point_locator* cpl_;
        
        // no copy ctor
        corner_graph(warthog::graph::corner_graph& other) : S_DUMMY_ID(0), T_DUMMY_ID(0) { }

        // below the line: variables and data structures and functions 
        // for connecting the start and target nodes to the rest of the graph
        // -------------------------------------------------------------------
        uint32_t s_graph_id_, s_grid_id_;
        uint32_t t_graph_id_, t_grid_id_;
        uint32_t S_DUMMY_ID, T_DUMMY_ID; 

        std::vector<uint32_t> t_incoming_;
        std::vector<warthog::graph::edge_iter> t_incoming_iters_;

        // maps gridmap ids to graph ids
        std::unordered_map<uint32_t, uint32_t> id_map_;

        void
        construct();

        void
        insert_start(uint32_t start_grid_id, uint32_t target_grid_id);

        void
        insert_target(uint32_t start_grid_id, uint32_t target_grid_id);

        void
        uninsert();
};

}
}

#endif
