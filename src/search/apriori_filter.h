#ifndef WARTHOG_APRIORI_FILTER_H
#define WARTHOG_APRIORI_FILTER_H

// search/apriori_filter.h
//
// Sometimes during search it is desriable to ignore certain nodes that
// are identified apriori as not useful for any query. This filter helps to 
// achieve the goal by keeping track of a single bit for each node in a 
// discrete graph. During search the filter can be queried about the 
// bit-state of each node and those nodes whose bit is set can be ignored 
// (not generatedor expanded).
// 
// @author: dharabor
// @created: 2016-07-19

#include "node_filter.h"
#include "constants.h"
#include <cstdint>
#include <cstring>

namespace warthog
{

class apriori_filter 
{
    public:
        apriori_filter(uint32_t num_nodes);
        virtual ~apriori_filter();

        // returns true if the node with identifier @param node_id 
        // is being filtered; otherwise false
        bool
        filter(warthog::search_node* n);

        void
        set_flag_true(uint32_t node_id);

        void
        set_flag_false(uint32_t node_id);

        bool
        get_flag(uint32_t node_id);

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
