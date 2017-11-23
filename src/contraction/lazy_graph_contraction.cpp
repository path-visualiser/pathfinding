#include "apriori_filter.h"
#include "euclidean_heuristic.h"
#include "flexible_astar.h"
#include "graph_expansion_policy.h"
#include "heap.h"
#include "lazy_graph_contraction.h"

#include <algorithm>
#include <functional>

bool
warthog::ch::operator<(const ch_pair& first, const ch_pair& second)
{
    return first.cval_ < second.cval_;
}

warthog::ch::lazy_graph_contraction::lazy_graph_contraction(
        warthog::graph::planar_graph* g)
    : graph_contraction(g)

{
    heuristic_ = new warthog::euclidean_heuristic(g);
    filter_ = new warthog::apriori_filter(get_graph()->get_num_nodes());
    expander_ = new 
        warthog::graph_expansion_policy<warthog::apriori_filter>(g, filter_);

    alg_ = new flexible_astar<
        warthog::euclidean_heuristic,
         warthog::graph_expansion_policy<warthog::apriori_filter>>
             (heuristic_, expander_);

    uint32_t sz_g = get_graph()->get_num_nodes();
    heap_ = new warthog::heap<ch_pair>(sz_g, true);
    hn_pool_ = new heap_node<ch_pair>[sz_g];

    terms_ = new node_order_terms[sz_g];

    ws_max_expansions_ = warthog::INF;
}

warthog::ch::lazy_graph_contraction::~lazy_graph_contraction()
{
    delete [] terms_;
    terms_ = 0;

    delete [] hn_pool_;
    hn_pool_ = 0;

    delete heap_;
    heap_ = 0;

    delete alg_;
    alg_ = 0;

    delete filter_;
    filter_ = 0;

    delete heuristic_;
    heuristic_ = 0;

    delete expander_;
    expander_ = 0;
}

// create an initial ordering
void
warthog::ch::lazy_graph_contraction::preliminaries()
{
    std::cerr << "creating initial contraction order" << std::endl;

    total_searches_ = 0;
    total_expansions_ = 0;
    total_lazy_updates_ = 0;
   
    for(uint32_t i = 0;
        i < get_graph()->get_num_nodes(); i++)
    {
        std::cerr << i << " / " << get_graph()->get_num_nodes() << "\r";
        terms_[i].ed_ = edge_difference(i,true); // contraction value
        hn_pool_[i] = heap_node<ch_pair>(ch_pair(i, terms_[i].ed_));
        heap_->push(&hn_pool_[i]);
    }
    std::cerr << "all " << get_graph()->get_num_nodes() 
        << " nodes contracted " << std::endl;
}

void
warthog::ch::lazy_graph_contraction::postliminaries()
{
    std::cout << "lazy contraction summary:" << std::endl
        << "\twitness searches: " << total_searches_ << std::endl 
        << "\tlazy updates: "<< total_lazy_updates_ << std::endl
        << "\tnode expansions " << total_expansions_ << std::endl;
}

// identifies the node that should be contracted next NB: nodes are ranked in a
// lazy way; an initial priority is generated (lower is better) but priorities
// are not updated until a node is popped off the heap.  at this stage we
// re-compute the node priority and possibly insert it back into the heap if
// it is no longer the ``best'' candidate 
//
// @return the next node id or warthog::INF if all nodes have been contracted
uint32_t
warthog::ch::lazy_graph_contraction::next()
{   
    // early abort; all nodes contracted
    if(heap_->size() == 0)
    {
        return warthog::INF;
    }

    uint64_t num_expansions = total_expansions_;
    uint64_t num_searches = total_searches_;
    uint64_t num_lazy = total_lazy_updates_;

    // pop the best contraction candidate off the heap
    heap_node<ch_pair>* best;
    uint32_t best_id;
    while(true)
    {
        // verify the previously stored contraction value of
        // the best candidate (the value is lazily updated)
        best = heap_->pop();
        best_id = best->get_element().node_id_; 
        terms_[best_id].ed_ = edge_difference(best_id,false);
        int32_t cval = compute_contraction_value(best_id);
        if(heap_->size() > 0 && cval > heap_->peek()->get_element().cval_)
        {
           // hmm, no longer the best node; re-insert and try again
           hn_pool_[best_id].get_element().cval_ = cval;
           heap_->push(&hn_pool_[best_id]);
           total_lazy_updates_++;
           continue;
        }
        break;
    }
    
    warthog::graph::node* bn = get_graph()->get_node(best_id);
    std::set<uint32_t> updateset;

    // figure out which un-contracted nodes are adjacent to the best node
    for(warthog::graph::edge_iter it = bn->outgoing_begin();
            it != bn->outgoing_end(); it++)
    {
        uint32_t neighbour_id = (*it).node_id_;
        if(filter_->get_flag(neighbour_id)) { continue; }
        updateset.insert(neighbour_id);
    }

    for(warthog::graph::edge_iter it = bn->incoming_begin();
            it != bn->incoming_end(); it++)
    {
        uint32_t neighbour_id = (*it).node_id_;
        if(filter_->get_flag(neighbour_id)) { continue; }
        updateset.insert(neighbour_id);
    }

    // update un-contracted neighbour priority terms
    for(std::set<uint32_t>::iterator it = updateset.begin();
            it != updateset.end(); it++)
    {
        uint32_t neighbour_id = *(it);
        node_order_terms& nb_not = terms_[neighbour_id];

        // re-compute the edge difference estimate
        filter_->set_flag_true(best_id);
        nb_not.ed_ = edge_difference(neighbour_id,true);
        filter_->set_flag_false(best_id);

        // update the "search space size" estimate
        terms_[neighbour_id].depth_ = 
            std::max(nb_not.depth_, terms_[best_id].depth_ + 1);

        // increase the contracted neighbours count (== "deleted" neighbours)
        nb_not.nc_++;

        // recompute priority and update heap 
        int32_t cval = compute_contraction_value(neighbour_id);
        if(cval < hn_pool_[neighbour_id].get_element().cval_)
        {
           hn_pool_[neighbour_id].get_element().cval_ = cval;
           heap_->decrease_key(&hn_pool_[neighbour_id]);
        }
    }
    num_expansions = total_expansions_ - num_expansions;

    std::cerr << 
        "; neis: " << updateset.size() << 
        "; num_witness: " << (total_searches_-num_searches) << 
        "; num_lazy: " << (total_lazy_updates_-num_lazy) << 
        "; num_exps: " << num_expansions << 
        "; edv " << terms_[best_id].ed_ <<
        "; nc " << terms_[best_id].nc_ << 
        "; depth " << terms_[best_id].depth_ <<
        "; PRIORITY " << compute_contraction_value(best_id) <<
        std::endl;

    filter_->set_flag_true(best_id); // mark as contracted
    order_.push_back(best_id);
    return best_id;
}

// NB: assumes the via-node is already marked as contracted
// (and will thus not be expanded)
double
warthog::ch::lazy_graph_contraction::witness_search(
        uint32_t from_id, uint32_t to_id, double via_len)
{
    // only search for witness paths between uncontracted neighbours
    if(filter_->get_flag(from_id) || filter_->get_flag(to_id)) { return 0; }

    alg_->set_cost_cutoff(via_len);
    alg_->set_max_expansions_cutoff(ws_max_expansions_);
    warthog::graph::planar_graph* g = this->get_graph();

    // need to specify start + target ids using the identifier
    // that appears in the input file
    uint32_t ext_from_id = g->to_external_id(from_id);
    uint32_t ext_to_id = g->to_external_id(to_id);

    warthog::problem_instance pi(ext_from_id, ext_to_id);
    warthog::solution sol;
    alg_->get_path(pi, sol);
    total_expansions_ += sol.nodes_expanded_;
    total_searches_++;
    return sol.sum_of_edge_costs_;
}

// calculate the net number of edges that result from contracting
// a given node 
//
// @param n: the id of the node to contract
// @param estimate: when true, witness searches are truncated after
// a fixed number of node expansions. when false, witness searches
// are run in full. 
// 
int32_t
warthog::ch::lazy_graph_contraction::edge_difference(uint32_t node_id, bool estimate)
{
    ws_max_expansions_ = (estimate ? 1000 : warthog::INF);
    int32_t eadd = 0;
    warthog::graph::node* n = get_graph()->get_node(node_id);

    std::set<uint64_t> deleted;
    filter_->set_flag_true(node_id);
    for(int i = 0; i < n->out_degree(); i++)
    {
        warthog::graph::edge& out = *(n->outgoing_begin() + i);

        // skip contracted neighbours
        if(filter_->get_flag(out.node_id_)) { continue; }

        for(int j = 0; j < n->in_degree(); j++)
        {
            warthog::graph::edge& in = *(n->incoming_begin() + j);

            // skip contracted neighbours
            if(filter_->get_flag(in.node_id_)) { continue; }

            // skip reflexive arcs
            if(out.node_id_ == in.node_id_) { continue; }

            deleted.insert((uint64_t)&in);
            double via_len = in.wt_ + out.wt_;
            double witness_len  =
                    witness_search(in.node_id_, out.node_id_, via_len);

            // contraction will introduce a shortcut edge only if
            // the path <in, n, out> is the only shortest path
            if(witness_len > via_len)
            {
                eadd++;
            }
        }
        deleted.insert((uint64_t)&out);
    }
    filter_->set_flag_false(node_id);
    int32_t retval = eadd - deleted.size();
    return retval;
}

void
warthog::ch::lazy_graph_contraction::get_order(std::vector<uint32_t>& order)
{
    order = order_;

}

size_t
warthog::ch::lazy_graph_contraction::mem()
{
    return 
        graph_contraction::mem() + 
        heap_->mem() +
        alg_->mem() + 
        sizeof(*hn_pool_)*get_graph()->get_num_nodes() +
        sizeof(this);
}

// compute the contraction value of node @param nid
// this value represents how attractive a node is to contract
// (lower is better)
//
// it is computed by taking the weighted sum of a number of
// different metrics including edge difference, number of
// contracted neighbours and the maximum depth a search must
// reach to encounter and expand @param nid.
//
// the weights are magic numbers taken from Robert Geisberger's 
// 2008 thesis.
int32_t 
warthog::ch::lazy_graph_contraction::compute_contraction_value(uint32_t nid)
{
        return
            terms_[nid].ed_*190 + 
            terms_[nid].nc_*120 + 
            terms_[nid].depth_;
}
