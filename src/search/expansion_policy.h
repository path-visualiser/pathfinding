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
        expansion_policy(uint32_t nodes_pool_size);
        virtual ~expansion_policy();

        uint32_t
        get_nodes_pool_size() { return nodes_pool_size_; } 
       
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
		}

		inline void
		first(warthog::search_node*& ret, double& cost)
		{
            current_ = 0;
            n(ret, cost);
		}

		inline void
		n(warthog::search_node*& ret, double& cost)
		{
            if(current_ < neis_->size())
            {
                ret = (*neis_)[current_].node_;
                cost = (*neis_)[current_].cost_;
            }
            else
            {
                ret = 0;
                cost = 0;
            }
		}

		inline void
		next(warthog::search_node*& ret, double& cost)
		{
            current_++;
            n(ret, cost);
		}

        inline warthog::search_node* 
        last()
        {
            return (neis_->back()).node_;
        }

		virtual size_t
		mem()
		{
			return sizeof(*this) + 
            sizeof(neighbour_record) * neis_->capacity() + 
            nodepool_->mem();
		}

		virtual void 
		expand(warthog::search_node*, warthog::problem_instance*) = 0;

        virtual void
        get_xy(warthog::search_node*, int32_t& x, int32_t& y) = 0;

        // return a pointer to the memory allocated for node @param node_id
        // during search @param search_id. if the node was not generated 
        // during this search, or if the id is invalid, return 0
        warthog::search_node*
        get_ptr(uint32_t node_id, uint32_t search_id)
        {
            warthog::search_node* tmp = nodepool_->get_ptr(node_id);
            if(tmp && tmp->get_searchid() == search_id) 
            {
                return tmp;
            }
            return 0;
        }

    protected:
        inline void 
        add_neighbour(warthog::search_node* nei, double cost)
        {
            neis_->push_back(neighbour_record(nei, cost));
            //std::cout << " neis_.size() == " << neis_->size() << std::endl;
        }

        inline uint32_t 
        get_num_neighbours() { return neis_->size(); } 

    private:

        struct neighbour_record
        {
            neighbour_record(warthog::search_node* node, double cost)
            {
                node_ = node;
                cost_ = cost;
            }
            warthog::search_node* node_;
            double cost_;
        };

        warthog::blocklist* nodepool_;
        std::vector<neighbour_record>* neis_;
        uint32_t current_;
        uint32_t nodes_pool_size_;
};

}

#endif
