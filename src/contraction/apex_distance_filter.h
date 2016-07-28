#ifndef WARTHOG_APEX_DISTANCE_FILTER_H
#define WARTHOG_APEX_DISTANCE_FILTER_H

// contraction/apex_distance_filter.h
//
// @author: dharabor
// @created: 28-07-2016
//

#include <cstdint>
#include <ostream>
#include <vector>

namespace warthog
{

namespace graph
{
class planar_graph;
}

class search_node;
class apex_distance_filter
{
    public: 
        apex_distance_filter(
                const char* ddfile, 
                warthog::graph::planar_graph* g, 
                std::vector<uint32_t>* rank);

        apex_distance_filter(
                warthog::graph::planar_graph* g,
                std::vector<uint32_t>* rank);

        ~apex_distance_filter();

        bool 
        filter(warthog::search_node* n);

        double
        get_apex_distance(uint32_t);

        double
        get_down_distance(uint32_t);

        void
        print(std::ostream& out);

        void
        compute_apex_distance(uint32_t startid, uint32_t endid);

        inline void
        compute_apex_distance();


    private:
        std::vector<double>* apex_dist_;
        std::vector<double>* ddist_;

        warthog::graph::planar_graph* g_;
        std::vector<uint32_t>* rank_;
        uint32_t start_id_, last_id_;

        void
        load_apex_distance_values(const char* filename);
};
}

#endif

