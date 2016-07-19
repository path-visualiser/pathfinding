#ifndef WARTHOG_DUMMY_FILTER_H
#define WARTHOG_DUMMY_FILTER_H

// search/dummy_filter.h
//
// A node filter that doesn't filter anything. 
//
// @author: dharabor
// @created: 2016-07-19
//

#include "node_filter.h"
#include <cstdint>

namespace warthog
{

class dummy_filter : public node_filter
{
    public:
        dummy_filter() { } 
        virtual ~dummy_filter() { } 

        virtual bool
        filter(uint32_t node_id) { return false; }
};

}

#endif
