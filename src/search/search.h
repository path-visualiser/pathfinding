#ifndef WARTHOG_SEARCH_H
#define WARTHOG_SEARCH_H

// search.h
//
// A base class for functionality common to all
// search algorithms.
//
// @author: dharabor
// @created: 2016-02-24
//

#include <stdint.h>
#include <stdlib.h>

namespace warthog
{

class search
{
    public:

        search()
        {
            searchid_ = 0;
            reset_metrics();
        }

        virtual ~search() { }

		inline uint32_t 
		get_nodes_expanded() { return nodes_expanded_; }

		inline uint32_t
		get_nodes_generated() { return nodes_generated_; }

		inline uint32_t
		get_nodes_touched() { return nodes_touched_; }

		inline double
		get_search_time() { return search_time_; }

		inline bool
		get_verbose() { return verbose_; }

		inline void
		set_verbose(bool verbose) { verbose_ = verbose; } 

        inline void
        reset_metrics()
        {
			nodes_expanded_ = nodes_generated_ = nodes_touched_ = 0;
			search_time_ = 0;
        }
        
        virtual size_t
        mem() = 0;

    protected:
		uint32_t searchid_;
		uint32_t nodes_expanded_;
		uint32_t nodes_generated_;
		uint32_t nodes_touched_;
		double search_time_;
        bool verbose_;
};


}

#endif

