#ifndef WARTHOG_BB_FILTER_H
#define WARTHOG_BB_FILTER_H

// contraction/bb_filter.h
//
// A node filter making use of rectangular geometric 
// containers (i.e. bounding boxes)
//
// For theory see 
//
// [Wager & Willhalm, 2005, Geometric Containers for 
// Efficient Shortest Path Computation, 
// Journal of Experimental Algorithms, vol 10, pp 1-30]

#include "geom.h"
#include <vector>

namespace warthog
{

namespace graph
{

class planar_graph;

} 

class search_node;
class problem_instance;
class bb_filter
{

    public:
        bb_filter();

        ~bb_filter();

        inline void
        set_target_xy(int32_t tx, int32_t ty) { tx_ = tx; ty_ = ty; }

        inline bool 
        filter(uint32_t node_id, uint32_t edge_id)
        {
            // prune any node whose bounding-box 
            // does not contain the target
            if(labels_.at(node_id).at(edge_id).contains(tx_, ty_))
            {
                return false;
            }
            return true;
        }

        void
        print(std::ostream& out);

        bool
        load_labels(const char* filename, warthog::graph::planar_graph* g);

        inline warthog::geom::rectangle
        get_label(uint32_t node_id, uint32_t edge_id)
        {
            return labels_.at(node_id).at(edge_id);
        }

    private:
        int32_t tx_, ty_;
        std::vector<std::vector<warthog::geom::rectangle>> labels_;

        warthog::graph::planar_graph* g_;

};

}

#endif

