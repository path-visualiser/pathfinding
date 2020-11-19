#ifndef __ORACLE_LISTENER_H_
#define __ORACLE_LISTENER_H_

#include "cpd.h"
#include "graph_oracle.h"

namespace warthog
{

namespace cpd
{

class oracle_listener
{
  public:
    oracle_listener(warthog::cpd::graph_oracle* oracle,
                    warthog::sn_id_t* source_id,
                    std::vector<warthog::cpd::fm_coll>* s_row)
        : oracle_(oracle), source_id_(source_id), s_row_(s_row) {}

    virtual ~oracle_listener() {};

    virtual void
    generate_node(warthog::search_node *from, warthog::search_node *succ,
                  warthog::cost_t edge_cost, uint32_t edge_id) = 0;

    inline void
    expand_node(warthog::search_node* current) { }

    inline void
    relax_node(warthog::search_node* current) { }

  protected:
    graph_oracle* oracle_;
    warthog::sn_id_t* source_id_;
    std::vector<warthog::cpd::fm_coll>* s_row_;
};

// helps to precompute first-move data
class graph_oracle_listener : public oracle_listener
{
    using oracle_listener::oracle_listener;

    inline void
    generate_node(warthog::search_node *from, warthog::search_node *succ,
                  warthog::cost_t edge_cost, uint32_t edge_id) final
    {
        if(from == 0) { return; } // start node

        if(from->get_id() == *source_id_) // start node successors
        {
            //assert(s_row_.at(succ->get_id()) == 0);
            assert(edge_id <
                        oracle_->get_graph()->get_node(
                        (uint32_t)*source_id_)->out_degree());
            s_row_->at(succ->get_id()) = (1 << edge_id);
            assert(s_row_->at(succ->get_id()));
        }
        else // all other nodes
        {
            warthog::sn_id_t succ_id = succ->get_id();
            warthog::sn_id_t from_id = from->get_id();
            double alt_g = from->get_g() + edge_cost;
            double g_val =
                succ->get_search_number() == from->get_search_number() ?
                succ->get_g() : DBL_MAX;

            //  update first move
            if(alt_g < g_val)
            {
                s_row_->at(succ_id) = s_row_->at(from_id);
                assert(s_row_->at(succ_id) == s_row_->at(from_id));
            }

            // add to the list of optimal first moves
            if(alt_g == g_val)
            {
                s_row_->at(succ_id) |= s_row_->at(from_id);
                assert(s_row_->at(succ_id) >= s_row_->at(from_id));
            }

        }
    }
};

// helps to precompute first-move data, this time we build the rows in reverse.
class reverse_oracle_listener : public oracle_listener
{
    using oracle_listener::oracle_listener;

    inline void
    generate_node(warthog::search_node *from, warthog::search_node *succ,
                  warthog::cost_t edge_cost, uint32_t edge_id) final
    {
        if(from == nullptr) { return; } // start node has no predecessor

        double alt_g = from->get_g() + edge_cost;
        double g_val =
            succ->get_search_number() == from->get_search_number() ?
            succ->get_g() : DBL_MAX;
        // We record the optimal move towards a node which is the id of the
        // predecessor's edge
        graph::node* pred = oracle_->get_graph()->get_node(succ->get_id());
        graph::edge_iter eit = pred->find_edge(from->get_id());
        uint8_t fm = 1 << (eit - pred->outgoing_begin());

        assert(eit != pred->outgoing_end());

        //  update first move
        if(alt_g < g_val)
        {
            // Pointer arithmetic ftw
            s_row_->at(succ->get_id()) = fm;
            assert(s_row_->at(succ->get_id()));
        }

        // add to the list of optimal first moves
        if(alt_g == g_val)
        {
            s_row_->at(succ->get_id()) |= fm;
        }
    }

};

}

}

#endif // __ORACLE_LISTENER_H_
