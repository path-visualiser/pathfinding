#ifndef WARTHOG_APEX_FILTER_H
#define WARTHOG_APEX_FILTER_H

// contraction/apex_filter.h
//
// a filter for forward contraction hierarchies in cases where the
// apex of the path is known (or when its rank can be upper-bounded):
//
// @author: dharabor
// @created: 2017-05-10
//

#include <vector>

namespace warthog
{
// allows passing apex-related data to certain expanders that support it
class apex_filter
{
    public:

        apex_filter(std::vector<uint32_t>* order)
        {
            last_search_id_ = warthog::INF;
            apex_id_ = warthog::INF;
            apex_reached_ = false;
            order_ = order;
        }

        inline bool
        filter(warthog::search_node* n)
        {
            assert(n->get_parent());

            if(n->get_search_id() != last_search_id_)
            {
                // TODO: put in search algo and avoid one branch?
                apex_reached_ = false;
                last_search_id_ = n->get_search_id();
            }
            
            // have we reached the apex yet?
            if(n->get_parent()->get_id() == apex_id_) 
            { apex_reached_ = true; }

            uint32_t n_rank = order_->at(n->get_id());
            uint32_t p_rank = order_->at(n->get_parent()->get_id());
            if(apex_reached_)
            {
                // never follow up edges after the apex is reached
                if(n_rank > p_rank) { return true; }
                return false;
            }
            
            // never follow down edges before the apex is reached
            if(n_rank < p_rank) { return true; }


            // always follow up edges before the apex is reached
            return false;
        }

        inline void
        set_apex(uint32_t node_id)
        {
            apex_id_ = node_id;
        }

        inline uint32_t
        get_apex() { return apex_id_; }
    
    private:
        uint32_t last_search_id_;
        uint32_t apex_id_;
        bool apex_reached_;
        std::vector<uint32_t>* order_;
};
}

#endif
