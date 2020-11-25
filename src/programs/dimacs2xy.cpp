#include "cfg.h"
#include "contraction.h"
#include "dimacs_parser.h"
#include "geography.h"

#include <cmath>
#include <numeric>
#include <iostream>

// Reorder the edges of a node to match their bearing wrt North
//
// This function could be improved by not having to allocate a two vectors at
// each iteration -- e.g., pre-allocate them in 'orient_edges()' and handle
// bounds.
std::vector<size_t>
orient_node_edges(warthog::graph::xy_graph &g, warthog::graph::ECAP_T degree,
                  warthog::graph::edge_iter start)
{
    std::vector<double> angles(degree);

    for (uint32_t edge_idx = 0; edge_idx < degree; edge_idx++)
    {
        warthog::graph::edge* e = start + edge_idx;
        // Longitude and latitude of the point
        int32_t lng, lat;
        g.get_xy(e->node_id_, lat, lng);
        // Bearing wrt North
        //
        // Don't forget to convert DIMACS ids to degrees
        angles.at(edge_idx) = warthog::geo::true_bearing(
            lng / warthog::geo::DIMACS_RATIO, lat / warthog::geo::DIMACS_RATIO);
    }

    // Compute the index ordering given all edges' angles
    std::vector<size_t> idx(degree);
    std::iota(idx.begin(), idx.end(), 0);

    stable_sort(idx.begin(), idx.end(), [&](size_t i1, size_t i2)
    {
        return angles[i1] < angles[i2];
    });

    return idx;
}

// Orient graph edges where 0 is the first clockwise edge from the North.
void
orient_edges(warthog::graph::xy_graph& g)
{
    for (uint32_t i = 0; i < g.get_num_nodes(); i++)
    {
        warthog::graph::node* n = g.get_node(i);

        // So, clear does not free memory, so we can tag the old starting
        // memory, clear the node, and re-insert edges in the right order.
        warthog::graph::edge_iter head = n->outgoing_begin();
        warthog::graph::ECAP_T out_deg = n->out_degree();
        std::vector<warthog::graph::edge> save(out_deg);
        n->clear();

        for(warthog::graph::ECAP_T i = 0; i < out_deg; i++)
        {
            save.at(i) = *(head + i);
        }

        std::vector<size_t> order = orient_node_edges(g, out_deg, head);
        for(auto id : order)
        {
            n->add_outgoing(save.at(id));
        }
        assert(n->out_degree() == out_deg);
    }
}

void
help()
{
    std::cerr 
       << "Converts between the graph format used at the 9th DIMACS Implementation\n"
       << "Challenge and the xy graph format used by Warthog.\n"
       << "Usage: ./dimacs2xy --input [dimacs .co file] [dimacs .gr file]\n"
       << "\t--order enforce azimuth ordering of edges in the resulting graph\n";
}

int 
main(int argc, char** argv)
{
    int order = 0;
	// parse arguments
	warthog::util::param valid_args[] = 
	{
		{"core",  required_argument, 0, 2},
		{"input",  required_argument, 0, 2},
		{"order", no_argument, &order, 1},
		{0,  0, 0, 0}
	};

    warthog::util::cfg cfg;
	cfg.parse_args(argc, argv, "-hc:", valid_args);
    
    if(argc < 2)
    {
		help();
        exit(0);
    }

    std::string gr_file = cfg.get_param_value("input");
    std::string co_file = cfg.get_param_value("input");
    if(gr_file == "")
    {
        std::cerr << "err; missing --input [co file] [gr file]\n";
        return EINVAL;
    }

    if(gr_file == "")
    {
        std::cerr << "err; missing --input [co file] [gr file]\n";
        return EINVAL;
    }

    // load 
    warthog::dimacs_parser parser(co_file.c_str(), gr_file.c_str());

    // convert
    warthog::graph::xy_graph g_xy;
    warthog::graph::dimacs_to_xy_graph(parser, g_xy);

    // Order edges by azimuth
    if(order)
    {
        orient_edges(g_xy);
    }

    // dump
    std::cout << g_xy;
    return 0;
}

