#ifndef WARTHOG_NODE_FILTER_H
#define WARTHOG_NODE_FILTER_H

// node_filter.h
//
// Sometimes during search it is desriable to skip over
// certain nodes that would otherwise need to be generated. 
//
// A node filter helps to achieve this goal. The idea is to push all 
// the decision-making into a filter class and combine the filter with a 
// generic expansion policy and search strategy.
// 
// @author: dharabor
// @created: 2016-01-26
//

#include "constants.h"
#include <cstdint>
#include <cstring>

namespace warthog
{

class node_filter
{
    public:
        node_filter(uint32_t num_nodes);
        ~node_filter();

        // filter a single node
        void
        filter(uint32_t node_id);

        // stop filtering a single node
        void 
        unfilter(uint32_t node_id);

        // returns true if @param node_id is currently
        // being filtered; false otherwise.
        bool
        get_filter_flag(uint32_t node_id);

        // clear all filter flags
        void
        reset_filter();

        inline size_t
        mem()
        {
            return 
                sizeof(*filter_)*filter_sz_ +
                sizeof(this);
        }

    private:
        warthog::dbword* filter_;
        uint32_t filter_sz_;
};

}

#endif
