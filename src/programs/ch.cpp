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
int has_input=0;
int has_order=0;
warthog::util::cfg cfg;

void
help()
{
    std::cerr 
        << "=>  manual <== \n"
        << "This program creates a contraction hierarchy from an xy input graph \n"
        << "(a custom format, similar to that used at the 9th DIMACS challenge)\n";
	std::cerr << "The following are valid program parameters:\n"
	<< "\t--order [ fixed | lazy ]\n"
    << "\t--input [xy graph file]\n"
	<< "\t--verbose (optional; prints debugging info when compiled with debug symbols)\n"
	<< "\t--verify (optional; verify lazy priorities before contraction.\n"
    << "\t          slow but can produce less shortcut edges)\n";
}

void 
contract_graph()
{
    warthog::ch::ch_data chd(0, warthog::ch::UP_DOWN);
    std::string outfile;

    // load up the input graph/grid
    if(cfg.get_num_values("input") != 1)
    {
        std::string xy_file = cfg.get_param_value("input");
        std::ifstream ifs(xy_file);
        warthog::graph::read_xy(ifs, *chd.g_, true);
    }
    else
    {
        std::cerr << "err; input graph (--input [xy graph file]) not specified\n";
        return;
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
        contractor.contract(chd.g_, &node_order);

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
        contractor.contract(&chd, verify);

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
		{"input",  required_argument, &has_input, 1},
		{"order",  required_argument, &has_order, 1},
		{0,  0, 0, 0}
	};
	cfg.parse_args(argc, argv, "abc:d:", valid_args);

    if(!(has_input && has_order))
    {
		help();
        exit(0);
    }

    contract_graph();
}
