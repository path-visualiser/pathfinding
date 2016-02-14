#ifndef WARTHOG_BIDIRECTIONAL_SEARCH_H
#define WARTHOG_BIDIRECTIONAL_SEARCH_H

// bidirectional_search.h
//
// A customisable variant of best-first bidirectional search.
// Users can pass in any heuristic and any (domain-specific) expansion policy.
//
// @author: dharabor
// @created: 2016-02-14
//

#include <cstdlib>
#include <stack>
#include <stdint.h>

namespace warthog
{

class expansion_policy;
class pqueue;
class search_node;

typedef double (* heuristicFn)(warthog::search_node* first, warthog::search_node* second);
class bidirectional_search
{
    public:
        bidirectional_search(
                warthog::expansion_policy* expander, 
                warthog::heuristicFn heuristic);
        ~bidirectional_search();

        double
        get_length(uint32_t startid, uint32_t goal_id);
            
		inline std::stack<uint32_t>
		get_path(uint32_t startid, uint32_t goalid);

		inline bool
		get_verbose() { return verbose_; }

		inline void
		set_verbose(bool verbose) { verbose_ = verbose; } 

        size_t
        mem();

    private:
        warthog::pqueue* fwd_open_;
        warthog::pqueue* bwd_open_;
        warthog::expansion_policy* expander_;
        warthog::heuristicFn heuristic_;
        bool verbose_;

        warthog::search_node*
        search(uint32_t startid, uint32_t goalid);

		void
		cleanup();

};

}

#endif

