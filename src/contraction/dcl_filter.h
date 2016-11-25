#ifndef WARTHOG_DCL_FILTER_H
#define WARTHOG_DCL_FILTER_H

// contraction/dcl_filter.h
//
// Collects a variety of different edge labels that we can compute for 
// the down-closure of a node. Intended for use with forward-driven 
// contraction hierarchies.
//
// Current edge labels:
//  - down distance: the max optimal distance from the source to any
//  node in its down closure
//
//  - down level: the minimum level of any node in the down closure
//
//  - bounding box: a rectangular bounding box that contains every
//  node in the down closure
//
// @author: dharabor
// @created: 2016-08-6
// 

#include "geom.h"
#include <vector>
#include <cstdint>

namespace warthog
{

namespace graph
{

class planar_graph;
class edge;

} 

class euclidean_heuristic;
class problem_instance;
class search_node;
class dcl_filter 
{
    struct dcl_label
    {
        dcl_label()
        {
//            min_lvl_ = 0;
//            ddist_ = 0;
        }

        //std::vector<warthog::geom::rectangle> edgebb_;
        warthog::geom::rectangle bbox_;
//        uint32_t min_lvl_;
//        double ddist_;
    };

    public:
        dcl_filter(
                const char* ddfile, 
                warthog::graph::planar_graph* g, 
                std::vector<uint32_t>* rank);

        dcl_filter(
                warthog::graph::planar_graph* g,
                std::vector<uint32_t>* rank);

        ~dcl_filter();

        // return true if the ith successor (@param succ_index)
        // of node @param n does not appear on any optimal path
        // to the target location
        bool 
        filter(warthog::search_node* n, warthog::graph::edge* e);

        uint32_t
        get_level(uint32_t);

        void
        print(std::ostream& out);

        void
        compute(uint32_t startid, uint32_t endid);
                
        void
        compute();

        void
        init(warthog::problem_instance* instance);

        bool
        load_labels(const char* filename);


    private:
        // vectors to store each label
        std::vector<dcl_label> labels_;
        warthog::euclidean_heuristic* h_;

        warthog::graph::planar_graph* g_;
        std::vector<uint32_t>* rank_;
        uint32_t start_id_, last_id_;
        int32_t tx, ty;
        uint32_t target_lvl;


//        bool 
//        dd_filter(uint32_t node_id, warthog::graph::edge* e, double g_value);
//
//        bool
//        bb_filter(uint32_t node_id, warthog::graph::edge* e);
//
//        bool
//        lvl_filter(uint32_t node_id, warthog::graph::edge* e);

        bool
        edgebb_filter(uint32_t node_id, warthog::graph::edge* e);
};

}

#endif

