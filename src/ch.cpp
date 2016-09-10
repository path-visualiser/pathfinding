#include "cfg.h"
#include "ch_expansion_policy.h"
#include "dimacs_parser.h"
#include "fixed_graph_contraction.h"
#include "graph.h"
#include "lazy_graph_contraction.h"
#include "planar_graph.h"

#include <iostream>
#include <string>

bool verbose=false;
warthog::util::cfg cfg;

void
help()
{
    std::cerr << 
        "create a contraction hierarchy from " <<
        "a given (currently, DIMACS-format only) input graph\n";
	std::cerr << "valid parameters:\n"
    << "\t--dimacs [gr file] [co file] (IN THIS ORDER!!)\n"
	<< "\t--order [order-of-contraction file]\n"
	<< "\t--verbose (optional)\n";
}

void contract_graph()
{
    std::string grfile = cfg.get_param_value("dimacs");
    std::string cofile = cfg.get_param_value("dimacs");
    std::cerr << "param values " << std::endl;
    std::string orderfile = cfg.get_param_value("order");
    cfg.print_values("dimacs", std::cerr);

    std::cerr << "grfile: "<< grfile << " cofile " << cofile << std::endl;

    // contract a graph and create a hierarchy
    warthog::graph::planar_graph g;
    std::vector<uint32_t> order;

    if(orderfile.compare("") != 0)
    {
        // use an existing node order
        g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true);
        warthog::ch::load_node_order(orderfile.c_str(), order);
        warthog::ch::fixed_graph_contraction contractor(&g, &order);
        contractor.set_verbose(verbose);
        contractor.contract();

        // save the order of contraction
        contractor.get_order(order);
        std::string orderfile = grfile + ".ooc";
        std::cerr << "saving order to file " << orderfile << std::endl;
        warthog::ch::write_node_order(orderfile.c_str(), order);

        // save the result
        grfile.append(".ch");
        std::cerr << "saving contracted graph to file " << grfile << std::endl;
        std::fstream ch_out(grfile.c_str(), std::ios_base::out | std::ios_base::trunc);
        if(!ch_out.good())
        {
            std::cerr << "\nerror exporting ch to file " << grfile << std::endl;
        }
        g.print_dimacs_gr(ch_out);
        ch_out.close();
    }
    else
    {
        // create a new contraction hierarchy with dynamic node ordering
        g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true);
        warthog::ch::lazy_graph_contraction contractor(&g);
        contractor.set_verbose(verbose);
        contractor.contract();

        // save the order of contraction
        contractor.get_order(order);
        std::string orderfile = grfile + ".ooc";
        std::cerr << "saving order to file " << orderfile << std::endl;
        warthog::ch::write_node_order(orderfile.c_str(), order);

        // save the result
        grfile.append(".ch");
        std::cerr << "saving contracted graph to file " << grfile << std::endl;
        std::fstream ch_out(grfile.c_str(), std::ios_base::out | std::ios_base::trunc);
        if(!ch_out.good())
        {
            std::cerr << "\nerror exporting ch to file " << grfile << std::endl;
        }
        g.print_dimacs_gr(ch_out);
        ch_out.close();
    }
    std::cerr << "all done!\n";
}

int main(int argc, char** argv)
{

	// parse arguments
    int print_help=false;
	warthog::util::param valid_args[] = 
	{
		{"help", no_argument, &print_help, 1},
		{"verbose", no_argument, &print_help, 1},
		{"dimacs",  required_argument, 0, 2},
		{"order",  required_argument, 0, 1}
	};
	cfg.parse_args(argc, argv, "-hvd:o:", valid_args);

    if(argc == 1 || print_help)
    {
		help();
        exit(0);
    }

    if(cfg.get_num_values("dimacs") != 2)
    {
        std::cerr << "insufficient values for param --dimacs (need gr and co files)\n";
        exit(0);
    }

    contract_graph();
}
