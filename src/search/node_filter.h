#ifndef WARTHOG_NODE_FILTER_H
#define WARTHOG_NODE_FILTER_H

// search/node_filter.h
// 
// Implement this class to filter or prune nodes during
// search.
//
// @author: dharabor
// @created: 2016-01-26

#include <cstdint>
namespace warthog
{

class search_node;
class node_filter
{
    public:
        virtual 
        bool filter(warthog::search_node*) = 0;

        uint32_t
        mem() { return 0; }
};

}

#endif
