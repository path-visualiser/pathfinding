#ifndef WARTHOG_EXPANSION_POLICY_H
#define WARTHOG_EXPANSION_POLICY_H

// expansion_policy.h
//
// An (abstract) expansion policy for searching explicit
// domains. It is assumed the number of nodes in the 
// search space is known apriori and a description of each node can be
// generated in constant time and independent of any other node.
//
// @author: dharabor
// @created: 2016-01-26
//

#include "blocklist.h"
#include "search_node.h"
#include "problem_instance.h"

namespace warthog
{

class expansion_policy
{
    public:
        expansion_policy(uint32_t num_nodes);
        virtual ~expansion_policy();
       
		// create a warthog::search_node object from a state description
		// (in this case, an id)
		inline warthog::search_node*
		generate(uint32_t node_id)
		{
			return nodepool_->generate(node_id);
		}

        inline void
        clear()
        {
            reset();
            nodepool_->clear();
        }        

		inline void
		reset()
		{
			current_ = 0;
            neis_->clear();
            costs_->clear();
		}

		inline void
		first(warthog::search_node*& ret, warthog::cost_t& cost)
		{
            current_ = 0;
            n(ret, cost);
		}

		inline void
		n(warthog::search_node*& ret, warthog::cost_t& cost)
		{
            if(current_ < neis_->size())
            {
                ret = (*neis_)[current_];
                cost = (*costs_)[current_];
            }
            else
            {
                ret = 0;
                cost = 0;
            }
		}

		inline void
		next(warthog::search_node*& ret, warthog::cost_t& cost)
		{
            current_++;
            n(ret, cost);
		}

        inline warthog::search_node* 
        last()
        {
            return neis_->back();
        }

		virtual uint32_t
		mem()
		{
			return sizeof(*this) + 
            sizeof(warthog::search_node*) * neis_->capacity() + 
            sizeof(warthog::cost_t) * costs_->capacity() +
            nodepool_->mem();
		}

		virtual void 
		expand(warthog::search_node*, warthog::problem_instance*) = 0;

        virtual void
        get_xy(warthog::search_node*, int32_t& x, int32_t& y) = 0;

    protected:
        inline void 
        add_neighbour(warthog::search_node* nei, warthog::cost_t cost)
        {
            neis_->push_back(nei);
            costs_->push_back(cost);
            //nei->print(std::cout);
            //std::cout << " neis_.size() == " << neis_->size() << std::endl;
        }

    private:
        warthog::blocklist* nodepool_;
        std::vector<warthog::search_node*>* neis_;
        std::vector<warthog::cost_t>* costs_;
        uint32_t current_;
};

}

#endif
