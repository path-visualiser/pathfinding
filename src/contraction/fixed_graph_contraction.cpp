#include "apriori_filter.h"
#include "constants.h"
#include "euclidean_heuristic.h"
#include "fixed_graph_contraction.h"
#include "flexible_astar.h"
#include "graph_expansion_policy.h"
#include "planar_graph.h"

warthog::ch::fixed_graph_contraction::fixed_graph_contraction(
                warthog::graph::planar_graph* g)
    : graph_contraction(g)
{
    order_ = new std::vector<uint32_t>(get_graph()->get_num_nodes());
    warthog::ch::make_input_order(*g, *order_);
    init();
}

warthog::ch::fixed_graph_contraction::fixed_graph_contraction(
                warthog::graph::planar_graph* g,
                std::vector<uint32_t>* order)
    : graph_contraction(g)
{
    assert(get_graph()->get_num_nodes() == order->size());
    order_ = order;
    init();
}

warthog::ch::fixed_graph_contraction::~fixed_graph_contraction()
{
    delete alg_;
    delete filter_;
    delete heuristic_;
    delete expander_;
}


void
warthog::ch::fixed_graph_contraction::init()
{
    order_index_ = 0;

    expander_ = new warthog::graph_expansion_policy(get_graph());
    //warthog::apriori_filter* af = get_filter();
    //expander_->set_filter(af);

    heuristic_ = new warthog::euclidean_heuristic(get_graph());
    filter_ = new warthog::apriori_filter(get_graph()->get_num_nodes());
    alg_ = new flexible_astar<
                    warthog::euclidean_heuristic,
                    warthog::graph_expansion_policy,
                    warthog::apriori_filter>
                        (heuristic_, expander_, filter_);
}

void
warthog::ch::fixed_graph_contraction::preliminaries()
{
    total_searches_ = 0;
    total_expansions_ = 0;
}

uint32_t 
warthog::ch::fixed_graph_contraction::next()
{
    if(order_index_ < order_->size())
    {
        filter_->set_flag_true(order_index_); // mark as contracted
        return order_->at(order_index_++);
    }
    return warthog::INF;
}

// NB: assumes the via-node is already contracted
double
warthog::ch::fixed_graph_contraction::witness_search(
        uint32_t from_id, uint32_t to_id, double via_len)
{
    // only search for witness paths between uncontracted neighbours
    if(filter_->get_flag(from_id) || filter_->get_flag(to_id)) { return 0; }

    // run the search
    alg_->set_cost_cutoff(via_len);
    warthog::problem_instance pi(from_id, to_id);
    double witness_len = alg_->get_length(pi);
    total_expansions_ += alg_->get_nodes_expanded();
    total_searches_++;
    return witness_len;
}

size_t
warthog::ch::fixed_graph_contraction::mem()
{
    return 
        graph_contraction::mem() + 
        alg_->mem() + 
        sizeof(this);
}
