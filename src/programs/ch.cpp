#include "cfg.h"
#include "bch_expansion_policy.h"
#include "dimacs_parser.h"
#include "fixed_graph_contraction.h"
#include "graph.h"
#include "lazy_graph_contraction.h"
#include "xy_graph.h"

#include <iostream>
#include <string>

int verbose=false;
int verify=false;
warthog::util::cfg cfg;

void
help()
{
    std::cerr << 
        "create a contraction hierarchy from " <<
        "a input graph (in the 9th DIMACS format) or a gridmap (in GPPC format)\n";
	std::cerr << "valid parameters:\n"
	<< "\t--order [ fixed | lazy ]\n"
    << "\t--partial [1-100] (optional; percentage of nodes to contract)\n"
    << "\t--dimacs [gr file] [co file] (IN THIS ORDER!!)\n"
    << "\t--gridmap [map file] \n"
	<< "\t--verbose (optional)\n"
	<< "\t--verify (verify lazy node priorities before contraction)\n";
}

void 
contract_graph()
{
    warthog::ch::ch_data chd(0, warthog::ch::UP_DOWN);
    std::string outfile;

    // load up the input graph/grid
    warthog::graph::xy_graph g;
    if(cfg.get_num_values("dimacs") == 2)
    {
        std::string grfile = cfg.get_param_value("dimacs");
        std::string cofile = cfg.get_param_value("dimacs");

        if(!chd.g_->load_from_dimacs(grfile.c_str(), cofile.c_str(), false, true))
        {
            std::cerr 
                << "err; could not load gr or co input files (one or both)\n";
            return;
        }
        outfile = grfile + ".chd";
    }
    else if(cfg.get_num_values("gridmap") == 1)
    {
        std::string mapfile = cfg.get_param_value("gridmap");
        warthog::gridmap gm(mapfile.c_str());
        if(!chd.g_->load_from_grid(&gm, true))
        {
            std::cerr << "err; could not load gridmap file\n";
            return;
        }
        outfile = mapfile + ".chd";
    }
    else
    {
        std::cerr << "err; input graph or gridmap not specified\n";
        return;
    }

    // partial contraction means we process only the first k% of nodes
    std::string partial = cfg.get_param_value("partial");
    int32_t pct_nodes_to_contract = 100;
    if(partial != "")
    {
        pct_nodes_to_contract = atoi(partial.c_str());
    }

    // use an existing node order for contraction operations
    std::string order_type = cfg.get_param_value("order");
    if(order_type == "fixed")
    {
        std::string orderfile = cfg.get_param_value("order");
        if(orderfile == "")
        {
            std::cerr 
                  << "err; fixed order requires input file."
                  << " syntax: --order fixed [node order file]\n";
            return;
        }

        std::vector<uint32_t> node_order;
        if(!warthog::ch::load_node_order(orderfile.c_str(), node_order))
        {
            std::cerr << "err; could not load node order input file\n";
            return;
        }


        warthog::ch::fixed_graph_contraction contractor;
        contractor.set_verbose(verbose);
        contractor.contract(chd.g_, &node_order, (uint32_t)pct_nodes_to_contract);

        // assign a level to every node based on its contraction order
        chd.level_->resize(chd.g_->get_num_nodes(), chd.g_->get_num_nodes()-1);
        for(uint32_t i = 0; i < node_order.size(); i++)
        {
            chd.level_->at(node_order.at(i)) = i;
        }

        // save the result
        std::cerr << "saving contracted graph to file " << outfile << std::endl;
        if(!warthog::ch::save_ch_data(outfile.c_str(), &chd))
        {
            std::cerr << "\nerror writing to file " << outfile << std::endl;
        }
    }
    else if(order_type == "lazy")
    {
        // create a new contraction hierarchy with dynamic node ordering
        warthog::ch::lazy_graph_contraction contractor;
        contractor.set_verbose(verbose);
        contractor.contract(&chd, verify, (uint32_t)pct_nodes_to_contract);

    }
    else
    {
        std::cerr << "unknown parameter for --order\n";
        return;
    }

    // save the result
    std::cerr << "saving contracted graph to file " << outfile << std::endl;
    if(!warthog::ch::save_ch_data(outfile.c_str(), &chd))
    {
        std::cerr << "error exporting ch to file " << outfile << std::endl;
    }
    std::cerr << "all done!\n";
}

int main(int argc, char** argv)
{

	// parse arguments
	warthog::util::param valid_args[] = 
	{
		{"verbose", no_argument, &verbose, 1},
		{"verify", no_argument, &verify, 1},
		{"dimacs",  required_argument, 0, 2},
		{"gridmap",  required_argument, 0, 1},
		{"order",  required_argument, 0, 3},
		{"partial",  required_argument, 0, 4}
	};
	cfg.parse_args(argc, argv, "-hvd:o:", valid_args);

    if(argc == 1)
    {
		help();
        exit(0);
    }

    contract_graph();
}
