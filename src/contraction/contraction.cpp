#include "contraction.h"
#include "euclidean_heuristic.h"
#include "flexible_astar.h"
#include "graph_expansion_policy.h"
#include "planar_graph.h"

void
warthog::contraction::make_input_order(warthog::graph::planar_graph& g, std::vector<uint32_t>& order)
{
    order.clear();
    order.reserve(g.get_num_nodes());
    for(uint32_t i=0; i < g.get_num_nodes(); i++)
    {
        order.push_back(i);
    }
}

void 
warthog::contraction::write_node_order(const char* filename, std::vector<uint32_t>& order)
{
    std::ofstream ofs(filename, std::ios_base::out | std::ios_base::trunc);
    if(ofs.good())
    {
        for(uint32_t i = 0; i < order.size(); i++)
        {
            ofs << order.at(i) << "\n";
        }
    }
    ofs.close();
}

void 
warthog::contraction::load_node_order(const char* filename, std::vector<uint32_t>& order)
{
    order.clear();
    std::ifstream ifs(filename, std::ios_base::in);
    if(!ifs.good())
    {
        std::cerr << "\nerror trying to load node order from file " 
            << filename << std::endl;
        ifs.close();
        return;
    }

    while(ifs.good())
    {
        uint32_t tmp;
        ifs >> tmp;
        order.push_back(tmp);
    }
    ifs.close();
}

void
warthog::contraction::contract(warthog::graph::planar_graph& g,
        std::vector<uint32_t>& lexord, bool verbose)
{
    assert(g.get_num_nodes() == lexord.size());

    // sort the nodes for contraction using their
    // associated lexical order
    std::vector<uint32_t> conord;
    conord.reserve(lexord.size());
    //FIXME!!! NO SORTING HERE!!!
    for(uint32_t i = 0; i < lexord.size(); i++)
    {
        conord.push_back(lexord[i]);
    }

    warthog::graph_expansion_policy expander(&g);
    warthog::euclidean_heuristic h(&g);
    warthog::flexible_astar<warthog::euclidean_heuristic, 
        warthog::graph_expansion_policy> astar(&h, &expander);

    std::cerr << "contracting graph " << g.get_filename() << std::endl;
    uint32_t total_searches = 0;
    uint32_t total_expansions = 0;
    uint32_t edges_before = g.get_num_edges();
    for(std::vector<uint32_t>::iterator it = conord.begin(); 
            it != conord.end(); it++)
    {
        uint32_t cid = (*it);
        assert(cid < g.get_num_nodes());
        warthog::graph::node* n = g.get_node(cid);
        std::cerr << "processing node " << (it - conord.begin()) 
            << "/" << conord.size() << " local searches: " << total_searches
            << " total expansions: " << total_expansions << "\r";

        expander.filter(cid); // never expand this node again
        for(int i = 0; i < n->out_degree(); i++)
        {
            warthog::graph::edge out = *(n->outgoing_begin() + i);
            if(lexord[out.node_id_] < lexord[cid]) { continue; }

            for(int j = 0; j < n->in_degree(); j++)
            {
                warthog::graph::edge in = *(n->incoming_begin() + j);
                if(lexord[out.node_id_] < lexord[cid]) { continue; }

                if(out.node_id_ == in.node_id_) { continue; }

                double via_len = in.wt_ + out.wt_;
                double witness_len = warthog::INF;

                // terminate when we prove no witness with len < via_len
                astar.set_cost_cutoff(via_len);
                witness_len = astar.get_length(
                        in.node_id_, out.node_id_);
                total_expansions += astar.get_nodes_expanded();
                total_searches++;

//                std::cerr << " node " << cid << " in-deg " << n->in_degree() << " out-deg " << n->out_degree() << " neis " << in.node_id_ 
//                    << " " << out.node_id_ << " via_len " << via_len <<
//                    " expansions " << astar.get_nodes_expanded() << std::endl;

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
                    if(via_len > 10000 && astar.get_nodes_expanded() <= 2)
                    {
                        witness_len = astar.get_length(in.node_id_, out.node_id_);
                    }

//                    std::cerr << " shortcut! node " << cid << " in-deg " << n->in_degree() << " out-deg " 
//                        << n->out_degree() << " shortcut " << in.node_id_ 
//                        << " " << out.node_id_ << " via_len " << via_len <<
//                        " expansions " << astar.get_nodes_expanded() << std::endl;
//                    
                    warthog::graph::node* tail = g.get_node(in.node_id_);
                    tail->add_outgoing(
                            warthog::graph::edge(out.node_id_, via_len));
                    warthog::graph::node* head = g.get_node(out.node_id_);
                    head->add_incoming(
                            warthog::graph::edge(in.node_id_, via_len));
                }
            }
        }
    }
    std::cerr << "graph, contracted";
    std::cerr << "edges before " << edges_before 
        << "; edges after " << g.get_num_edges() << std::endl;
}

