#ifndef WARTHOG_GRAPH_H
#define WARTHOG_GRAPH_H

// graph.h
//
// Common data structures for graphs.
//
// @author: dharabor
// @created: 2016-05-02
//

#include <cassert>
#include <cstdint>
#include <cfloat>
#include <iostream>

#include "constants.h"

namespace warthog
{
namespace graph
{

// adjustable limit on the maximum edge capacity (== degree) of a node. 
// smaller values make node manipulation operations more cache friendly
typedef uint16_t ECAP_T;
const uint16_t ECAP_MAX = UINT16_MAX;

class edge
{
    public:
        edge() { node_id_ = UINT32_MAX; wt_ = DBL_MAX; }

        edge(uint32_t node_id, double wt)
        {
            node_id_ = node_id;
            wt_ = wt;
        }

        edge(const warthog::graph::edge& other) 
        { node_id_ = other.node_id_; wt_ = other.wt_; }

        edge&
        operator=(const warthog::graph::edge& other) 
        { 
            node_id_ = other.node_id_; wt_ = other.wt_; 
            return *this;
        }

        // we store only one id, typically the head node
        // the id of the other node is derived contextually
        uint32_t node_id_;
        double wt_;
};
typedef edge* edge_iter;

class node
{
    public:
        node()
        {
           init(0, 0);
        };

        node(const node& other)
        {
            *this = other;
        }

        ~node()
        {
            delete [] incoming_;
            incoming_ = 0;
            in_cap_ = in_deg_ = 0;

            delete [] outgoing_;
            outgoing_ = 0;
            out_cap_ = out_deg_ = 0;
        }

        warthog::graph::node&
        operator=(warthog::graph::node&& other)
        {
            this->in_deg_ = other.in_deg_;
            this->out_deg_ = other.out_deg_;

            other.out_deg_ = 0;
            other.in_deg_ = 0;

            this->incoming_ = other.incoming_;
            other.incoming_ = 0;

            this->outgoing_ = other.outgoing_;
            other.outgoing_ = 0;
            return *this;
        }

        warthog::graph::node&
        operator=(const warthog::graph::node& other)
        {
            in_deg_ = other.in_deg_;
            in_cap_ = other.in_cap_;
            incoming_ = new edge[other.in_cap_];
            for(ECAP_T i = 0; i < other.in_deg_; i++)
            {
                incoming_[i] = other.incoming_[i];
            }

            out_deg_ = other.out_deg_;
            out_cap_ = other.out_cap_;
            outgoing_ = new edge[other.out_cap_];
            for(ECAP_T i = 0; i < other.out_deg_; i++)
            {
                outgoing_[i] = other.outgoing_[i];
            }
            return *this;
        }

        inline void
        add_incoming(warthog::graph::edge e)
        {
            add_edge(e, in_cap_, in_deg_, incoming_);
        }

        inline void
        add_outgoing(warthog::graph::edge e)
        {
            add_edge(e, out_cap_, out_deg_, outgoing_);
        }

        inline void
        del_incoming(warthog::graph::edge_iter iter)
        {
            del_edge(iter, in_deg_, incoming_);
        }

        inline void
        del_outgoing(warthog::graph::edge_iter iter)
        {
            del_edge(iter, out_deg_, outgoing_);
        }

        inline edge_iter
        incoming_begin()
        {
            return incoming_;
        }

        inline edge_iter
        incoming_end()
        {
            return &incoming_[in_deg_];
        }

        inline edge_iter
        outgoing_begin()
        {
            return outgoing_;
        }

        inline edge_iter
        outgoing_end()
        {
            return (out_deg_ == 0 ? 0 : &outgoing_[out_deg_]);
        }

        inline ECAP_T
        in_degree()
        {
            return in_deg_;
        }

        inline ECAP_T
        out_degree()
        {
            return out_deg_;
        }

        // find the index of an outgoing edge whose head 
        // node is @param to_id. the search begins from 
        // the edge iterator specified in @param it
        // (by default, the first outgoing edge).
        // the function returns an iterator for the
        // target edge or ::outgoing_end() if no suh edge
        warthog::graph::edge_iter
        find_edge(uint32_t to_id, warthog::graph::edge_iter it = 0)
        {
            if(it == 0) { it = outgoing_begin(); }
            if(it < outgoing_begin() || it >= outgoing_end())
            {
                return outgoing_end();
            }

            for( ; it < outgoing_end(); it++)
            {
                if((*it).node_id_ == to_id)
                {
                    return it;
                }
            }
            return outgoing_end();
        }

        inline size_t
        mem()
        {
            return 
                sizeof(warthog::graph::edge) * in_cap_ + 
                sizeof(warthog::graph::edge) * out_cap_ +
                sizeof(this);
        }

    private: 
        void
        init(ECAP_T in_capacity, ECAP_T out_capacity)
        {
            outgoing_ = 0;
            out_deg_ = out_cap_ = 0;
            incoming_ = 0;
            in_deg_ = in_cap_ = 0;

            if(in_capacity > 0)
            {
                incoming_ = new edge[in_capacity];
                in_cap_ = in_capacity;
            }

            if(out_capacity > 0)
            {
                outgoing_ = new edge[out_capacity];
                out_cap_ = out_capacity;
            }
        }

        // increase max (incoming or outgoing) edges that can be 
        // stored along with this node
        ECAP_T
        increase_capacity(ECAP_T newcap, ECAP_T oldcap, edge*& collection)
        {
            newcap = std::max<int>(1, newcap);
            assert(newcap > oldcap);
            if(newcap <= oldcap) { return oldcap; }

            edge* newcollection = new edge[newcap];
            for(int i = 0; i < oldcap; i++)
            {
                newcollection[i] = collection[i];
            }
            delete [] collection;
            collection = newcollection;
            return newcap;
        }

        // delete an edge and shift the remaining edges
        // to plug the hole
        inline void
        del_edge(warthog::graph::edge_iter elt, ECAP_T& deg, edge*& elts)
        {
            uint32_t index = elt - elts;
            if(index < deg)
            {
                deg--;
                for( ; index < deg; index++)
                {
                    elts[index] = elts[index+1];
                }
            }
        }

        inline void
        add_edge(warthog::graph::edge& e, 
                ECAP_T& max_elts, ECAP_T& deg, edge*& elts)
        {
            if(deg == ECAP_MAX)
            {
                std::cerr << "warthog::graph::node edge-capacity reached!\n";
                return;
            }
         
            // don't add redundant edges; we only want one:
            // the one with lowest cost
            for(uint32_t i = 0; i < deg; i++)
            {
                if(elts[i].node_id_ == e.node_id_)
                {
                    if(e.wt_ < elts[i].wt_)
                    {
                        elts[i].wt_ = e.wt_;
                    }
                    return;
                }
            }

            if(deg == max_elts)
            {
                max_elts = increase_capacity(max_elts*2, max_elts, elts);
            }
            elts[deg] = e; 
            deg++;
        }


    private:
        edge* incoming_;
        ECAP_T in_deg_;
        ECAP_T in_cap_;

        edge* outgoing_;
        ECAP_T out_deg_;
        ECAP_T out_cap_;
};
typedef node* node_iter;

}
}

#endif
