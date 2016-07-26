#ifndef WARTHOG_DOWN_DISTANCE_FILTER_H
#define WARTHOG_DOWN_DISTANCE_FILTER_H

#include <vector>
#include <cstdint>

namespace warthog
{

namespace graph
{

class planar_graph;

} 

class search_node;
class down_distance_filter 
{
    public:
        down_distance_filter(
                const char* ddfile, 
                warthog::graph::planar_graph* g, 
                std::vector<uint32_t>* rank);

        down_distance_filter(
                warthog::graph::planar_graph* g,
                std::vector<uint32_t>* rank);

        ~down_distance_filter();

        bool 
        filter(warthog::search_node* n);

        void
        print(std::ostream& out);

        void
        compute_down_distance(uint32_t startid, uint32_t endid);

        inline void
        compute_down_distance()
        {
            if(g_ && rank_)
            {
                compute_down_distance(0, g_->get_num_nodes());
            }
        }


    private:
        std::vector<double>* ddist_;
        warthog::graph::planar_graph* g_;
        std::vector<uint32_t>* rank_;
        uint32_t start_id_, last_id_;

        void
        load_down_distance_values(const char* filename);
};

}

#endif

