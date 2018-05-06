#include "apriori_filter.h"
#include "constants.h"
#include "euclidean_heuristic.h"
#include "fixed_graph_contraction.h"
#include "flexible_astar.h"
#include "graph_expansion_policy.h"
#include "planar_graph.h"

warthog::ch::fixed_graph_contraction::fixed_graph_contraction(
                warthog::graph::planar_graph* g)
    : g_(g)
{
    order_ = new std::vector<uint32_t>(get_graph()->get_num_nodes());
    warthog::ch::make_input_order(*g, *order_);
    init();
}

warthog::ch::fixed_graph_contraction::fixed_graph_contraction(
                warthog::graph::planar_graph* g,
                std::vector<uint32_t>* order)
    : g_(g)
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
    done_ = false;
    verbose_ = false;
    c_pct_ = 100;
    order_index_ = 0;

    filter_ = new warthog::apriori_filter(get_graph()->get_num_nodes());
    expander_ = new warthog::graph_expansion_policy< warthog::apriori_filter >
        (get_graph(), filter_);

    heuristic_ = new warthog::euclidean_heuristic(get_graph());
    alg_ = new flexible_astar<
                    warthog::euclidean_heuristic,
                    warthog::graph_expansion_policy<warthog::apriori_filter>>
                        (heuristic_, expander_);
}

void
warthog::ch::fixed_graph_contraction::preliminaries()
{
    total_searches_ = 0;
    total_expansions_ = 0;
}
void

warthog::ch::fixed_graph_contraction::contract()
{
    if(done_) { return; }
    done_ = true;

    if(c_pct_ < 100)
    {
        std::cerr << "partially "
                  << "("<<c_pct_<<"% of nodes) ";
    }
    std::cerr << "contracting graph " << g_->get_filename() << std::endl;
    uint32_t edges_before = g_->get_num_edges_out();

    preliminaries();

    uint32_t total_nodes = g_->get_num_nodes();
    uint32_t num_contractions = 0;
//    warthog::apriori_filter* filter = get_filter();
    for(uint32_t cid = next(); cid != warthog::INF; cid = next())
    {
        
        uint32_t pct = (num_contractions / (double)g_->get_num_nodes()) * 100;
        if(pct >= c_pct_)
        { 
            std::cerr << "\npartial contraction finished " 
                      << "(processed "<< pct << "% of all nodes)";
            break; 
        }

        std::cerr << "\r " << pct << "%; " << ++num_contractions << " /  " << total_nodes;
        if(verbose_)
        {
            std::cerr << "; current: " << cid;
        }

        warthog::graph::node* n = g_->get_node(cid);
        //filter->filter(cid); // never expand this node again

        for(int i = 0; i < n->out_degree(); i++)
        {
            warthog::graph::edge& out = *(n->outgoing_begin() + i);

            // skip already-contracted neighbours
//            if(filter->filter(out.node_id_)) { continue; }

            for(int j = 0; j < n->in_degree(); j++)
            {
                warthog::graph::edge& in = *(n->incoming_begin() + j);

                // skip already-contracted neighbours
//                if(filter->filter(in.node_id_)) { continue; }

                // no reflexive arcs please
                if(out.node_id_ == in.node_id_) { continue; }

                // terminate when we prove there is no witness 
                // path with len <= via_len
                double via_len = in.wt_ + out.wt_;
                double witness_len = 
                    witness_search(in.node_id_, out.node_id_, via_len);

                if(witness_len > via_len)
                {
                    if(verbose_)
                    {
                        std::cerr << "\tshortcut " << in.node_id_ << " -> "
                            << cid << " -> " << out.node_id_;
                        std::cerr << " via-len " << via_len;
                        std::cerr << " witness-len " << witness_len << std::endl;
                    }

                    warthog::graph::node* tail = g_->get_node(in.node_id_);
                    tail->add_outgoing(
                            warthog::graph::edge(out.node_id_, via_len));
                    warthog::graph::node* head = g_->get_node(out.node_id_);
                    head->add_incoming(
                            warthog::graph::edge(in.node_id_, via_len));
                }
            }
        }
    }

    std::cerr << "\ngraph, contracted. ";
    std::cerr << "edges before " << edges_before 
        << "; edges after " << g_->get_num_edges_out() << std::endl;
    postliminaries();
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

    // pathfinding queries must specify an external start and target id
    // (i.e. as they appear in the input file)
    warthog::graph::planar_graph* g = this->get_graph();
    uint32_t ext_from_id = g->to_external_id(from_id);
    uint32_t ext_to_id = g->to_external_id(to_id);

    // run the search
    alg_->set_cost_cutoff(via_len);
    alg_->set_max_expansions_cutoff(1000);
    warthog::problem_instance pi(ext_from_id, ext_to_id);
    warthog::solution sol;
    alg_->get_distance(pi, sol);
    total_expansions_ += sol.nodes_expanded_;
    total_searches_++;
    return sol.sum_of_edge_costs_;
}

size_t
warthog::ch::fixed_graph_contraction::mem()
{
    return 
        alg_->mem() + 
        sizeof(this);
}
