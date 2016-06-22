#include "euclidean_heuristic.h"
#include "flexible_astar.h"
#include "graph_expansion_policy.h"
#include "heap.h"
#include "lazy_graph_contraction.h"
#include "node_filter.h"

bool
warthog::ch::operator<(const ch_pair& first, const ch_pair& second)
{
    return first.edval_ < second.edval_;
}

warthog::ch::lazy_graph_contraction::lazy_graph_contraction(
        warthog::graph::planar_graph* g)
    : graph_contraction(g)

{
    expander_ = new warthog::graph_expansion_policy(g);
    expander_->set_filter(get_filter());

    heuristic_ = new warthog::euclidean_heuristic(g);
    alg_ = new euc_astar(heuristic_, expander_);

    uint32_t sz_g = get_graph()->get_num_nodes();
    heap_ = new warthog::heap<ch_pair>(sz_g, true);
    hn_pool_ = new heap_node<ch_pair>[sz_g];
}

warthog::ch::lazy_graph_contraction::~lazy_graph_contraction()
{
    delete expander_;
    expander_ = 0;

    delete heuristic_;
    heuristic_ = 0;

    delete alg_;
    alg_ = 0;

    delete heap_;
    heap_ = 0;

    delete [] hn_pool_;
    hn_pool_ = 0;
}

// create an initial ordering
void
warthog::ch::lazy_graph_contraction::preliminaries()
{
    if(get_verbose())
    {
        std::cerr << "creating initial contraction order" << std::endl;
    }

    total_searches_ = 0;
    total_expansions_ = 0;
    for(uint32_t i = 0; i < get_graph()->get_num_nodes(); i++)
    {
        if(get_verbose())
        {
            std::cerr << i << " / " << get_graph()->get_num_nodes() << "\r";
        }
        int32_t edval = edge_difference(i);
        hn_pool_[i] = heap_node<ch_pair>(ch_pair(i, edval));
        heap_->push(&hn_pool_[i]);
    }
    if(get_verbose())
    {
        std::cerr << "all " << get_graph()->get_num_nodes() 
            << " nodes contracted " << std::endl;
    }
}

// return the node that should be contracted next NB: nodes are ranked in a
// lazy way; an initial priority is generated (lower is better) but priorities
// are not updated until a node is popped off the heap.  at this stage we
// re-compute the node priority and possibly insert it back into the heap if
// it is no longer the ``best'' candidate (i.e. has the smallest edge
// difference value)
uint32_t
warthog::ch::lazy_graph_contraction::next()
{   
    uint32_t iter = 0;
    uint32_t max_iter = (uint32_t)(heap_->size() * 0.1);
    heap_node<ch_pair>* best = heap_->pop();
    if(heap_->size() > 0)
    {
        while(best->get_element().edval_ > 
                heap_->peek()->get_element().edval_)
        {
            uint32_t edgediff = edge_difference(best->get_element().node_id_);
            best->get_element().edval_ = edgediff;
            heap_->push(best);
            iter++;
            best = heap_->pop();
        }
    }

    if(get_verbose() && iter > max_iter)
    {
        std::cerr << "warthog::dynamic_node_order "
            << "#iters: "<<iter<<"; max_iter: " << max_iter 
            << std::endl;
    }
    order_.push_back(best->get_element().node_id_);
    return best->get_element().node_id_;
}

double
warthog::ch::lazy_graph_contraction::witness_search(
        uint32_t from_id, uint32_t to_id, double via_len)
{
    alg_->set_cost_cutoff(via_len);
    double witness_len = alg_->get_length(from_id, to_id);
    total_expansions_ += alg_->get_nodes_expanded();
    total_searches_++;
    return witness_len;
}

// calculate the net number of edges that result from contracting
// node @param n
int32_t
warthog::ch::lazy_graph_contraction::edge_difference(uint32_t node_id)
{
    int32_t ediff = 0;
    warthog::graph::node* n = get_graph()->get_node(node_id);
    warthog::node_filter* filter = get_filter();

    std::set<uint64_t> deleted;
    filter->filter(node_id);
    for(int i = 0; i < n->out_degree(); i++)
    {
        warthog::graph::edge& out = *(n->outgoing_begin() + i);

        // skip contracted neighbours
        if(filter->get_filter_flag(out.node_id_)) { continue; }

        for(int j = 0; j < n->in_degree(); j++)
        {
            warthog::graph::edge& in = *(n->incoming_begin() + j);

            // skip contracted neighbours
            if(filter->get_filter_flag(in.node_id_)) { continue; }

            // skip reflexive arcs
            if(out.node_id_ == in.node_id_) { continue; }

            deleted.insert((uint64_t)&in);
            double via_len = in.wt_ + out.wt_;

            double witness_len = 
                witness_search(in.node_id_, out.node_id_, via_len);

            // contraction will introduce a shortcut edge only if
            // the path <in, n, out> is the only shortest path
            if(witness_len > via_len)
            {
                ediff++;
            }
        }
        deleted.insert((uint64_t)&out);
    }
    ediff -= deleted.size();
    filter->unfilter(node_id);
    return ediff;
}
void
warthog::ch::lazy_graph_contraction::get_order(std::vector<uint32_t>* order)
{
    order->clear();
    for(uint32_t i = 0; i < order_.size(); i++)
    {
        order->push_back(order_.at(i));
    }
}

size_t
warthog::ch::lazy_graph_contraction::mem()
{
    return 
        graph_contraction::mem() + 
        heap_->mem() +
        alg_->mem() + 
        sizeof(*hn_pool_)*get_graph()->get_num_nodes() +
        sizeof(uint32_t)*order_.capacity() + 
        sizeof(this);
}
