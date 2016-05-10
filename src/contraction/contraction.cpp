#include "contraction.h"
#include "euclidean_heuristic.h"
#include "flexible_astar.h"
#include "graph_expansion_policy.h"
#include "planar_graph.h"

void
warthog::contraction::contract(warthog::graph::planar_graph& g)
{
    std::vector<uint32_t> input_order;
    for(uint32_t i=0; i < g.get_num_nodes(); i++)
    {
        input_order.push_back(i);
    }
    warthog::contraction::contract(g, input_order);
}

void
warthog::contraction::contract(warthog::graph::planar_graph& g,
        std::vector<uint32_t>& order, bool verbose)
{
    assert(g.get_num_nodes() == order.size());

    warthog::graph_expansion_policy expander(&g);
    warthog::euclidean_heuristic h(&g);
    warthog::flexible_astar<warthog::euclidean_heuristic, 
        warthog::graph_expansion_policy> astar(&h, &expander);

    std::cerr << "contracting graph " << g.get_filename() << std::endl;
    uint32_t edges_before = g.get_num_edges();
    for(std::vector<uint32_t>::iterator it = order.begin(); 
            it != order.end(); it++)
    {
        uint32_t cid = (*it);
        assert(cid < g.get_num_nodes());
        warthog::graph::node* n = g.get_node(cid);
        std::cerr << "processing node " << (it - order.begin() + 1) 
            << "/" << order.size() << "\r";
        for(int i = 0; i < n->out_degree(); i++)
        {
            warthog::graph::edge out = *(n->outgoing_begin() + i);

            for(int j = 0; j < n->in_degree(); j++)
            {
                warthog::graph::edge in = *(n->incoming_begin() + j);

                if(out.node_id_ == in.node_id_) { continue; }
                double via_len = in.wt_ + out.wt_;
                double witness_len = warthog::INF;
                if(out.node_id_ > cid && in.node_id_ > cid)
                {
                    expander.filter(cid);
                    witness_len = astar.get_length(
                            in.node_id_, out.node_id_);
                    // also: if <= is used below, and two alternatives exist
                    // two shortcuts will be added (??) -- only if we allow 
                    // redundant edges though... so
                    // do not add redundant (>= existing cost) edges
                    //std::cout << "witness_len " << witness_len << " via-len " << via_len << std::endl;
                    if(witness_len >= via_len)
                    {
                        if(verbose)
                        {
                            std::cerr << "shortcut " << in.node_id_ << " -> "
                                << cid << " -> " << out.node_id_;
                            std::cerr << " via-len " << via_len;
                            std::cerr << " witness-len " << witness_len << std::endl;
                        }
                        
                        warthog::graph::node* tail = g.get_node(in.node_id_);
                        tail->add_outgoing(
                                warthog::graph::edge(out.node_id_, witness_len));
                        warthog::graph::node* head = g.get_node(out.node_id_);
                        head->add_incoming(
                                warthog::graph::edge(in.node_id_, witness_len));
                    }
                }
            }
        }
    }
    std::cerr << "graph, contracted";
    std::cerr << "edges before " << edges_before 
        << "; edges after " << g.get_num_edges() << std::endl;
}

