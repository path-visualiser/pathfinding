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
        bb_filter(
                const char* ddfile, 
                warthog::graph::planar_graph* g);

        bb_filter(
                warthog::graph::planar_graph* g);

        ~bb_filter();

        void
        init(warthog::problem_instance* instance);

        bool 
        filter(warthog::search_node* n, uint32_t edge_id);

        bool 
        filter__(uint32_t node_id, uint32_t edge_id);

        void
        print(std::ostream& out);

        void
        compute();

        void
        compute(uint32_t startid, uint32_t endid);

        void
        compute_ch(std::vector<uint32_t>* rank);

        void
        compute_ch( uint32_t startid, uint32_t endid, 
                    std::vector<uint32_t>* rank );
        bool
        load_labels(const char* filename);

    private:
        warthog::graph::planar_graph* g_;
        uint32_t start_id_, last_id_;
        int32_t tx, ty;
        std::vector<std::vector<warthog::geom::rectangle>> labels_;
};

}

#endif

