#ifndef WARTHOG_BB_AF_FILTER_H
#define WARTHOG_BB_AF_FILTER_H

// bb_af_filter.h
//
// combined arcflags and bounding box filter
//
// @author: dharabor
// @created: 2016-08-25
//

#include <vector>
#include <cassert>
#include <cstdint>
#include <iostream>

#include "geom.h"

namespace warthog
{

namespace graph
{

class planar_graph;

}

class problem_instance;
class search_node;

class bb_af_filter
{
    struct bbaf_label
    {
        bbaf_label()
        {
            flags_ = 0;
        }

        bbaf_label(const bbaf_label& other)
        {
            flags_ = other.flags_;
            bbox_ = other.bbox_;
        }

        uint8_t* flags_;
        warthog::geom::rectangle bbox_;

    };

    public:
        // @param g: the search graph
        // @param part: a partitioning of the graph nodes
        bb_af_filter(
                warthog::graph::planar_graph* g, 
                std::vector<uint32_t>* rank,
                std::vector<uint32_t>* part);

        bb_af_filter(
                warthog::graph::planar_graph* g, 
                std::vector<uint32_t>* rank,
                std::vector<uint32_t>* part,
                const char* bb_af_file);

        ~bb_af_filter();

        // return true if the ith edge of @param node_id
        // (as specified by @param edge_index) cannot possibly 
        // appear on any optimal path to the current target;
        // return false otherwise.
        inline bool 
        filter(uint32_t node_id, uint32_t edge_idx, bool down)
        {
            bbaf_label& label = labels_.at(node_id).at(edge_idx);
            bool retval = (label.flags_[t_byte_] & t_bitmask_);
            if(down)
            {
                retval = retval && label.bbox_.contains(tx_, ty_);
            }
//            std::cerr << "bbaf filter; af " << (label.flags_[t_byte_] & t_bitmask_)
//                << " bbox " << label.bbox_.contains(tx_, ty_) << "\n";
            return !retval; 
        }

        void
        set_instance(warthog::problem_instance* instance);

        void
        compute();

        void
        compute(uint32_t firstid, uint32_t lastid);

        void
        print(std::ostream& out);

    private:    
        uint32_t nparts_;
        warthog::graph::planar_graph* g_;
        std::vector<uint32_t>* part_;
        std::vector<uint32_t>* rank_;
        std::vector<std::vector<bbaf_label>> labels_;

        // sometimes we want to compute bb_af for just 
        // a subset of nodes; those in the range [firstid, lastid)
        uint32_t firstid_;
        uint32_t lastid_;

        // values to quickly extract the flag bit for the target at hand
        uint32_t bytes_per_label_;
        uint32_t t_byte_;
        uint32_t t_bitmask_;
        int32_t tx_, ty_;
        
        void 
        init(warthog::graph::planar_graph* g, 
                std::vector<uint32_t>* rank,
                std::vector<uint32_t>* part);

        void
        load_bb_af_file(const char* filename);
};

}

#endif

