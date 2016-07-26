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

class search_node;
class dummy_filter 
{
    public:
        dummy_filter() { } 
        ~dummy_filter() { } 

        bool
        filter(warthog::search_node*) { return false; }
};

}

#endif
