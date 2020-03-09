#include "cfg.h"
#include "contraction.h"
#include "dimacs_parser.h"
#include <iostream>

void
help()
{
    std::cerr 
       << "Converts between the graph format used at the 9th DIMACS Implementation\n"
       << "Challenge and the xy graph format used by Warthog.\n"
       << "Usage: ./dimacs2xy --input [dimacs .co file] [dimacs .gr file]\n";
}

int 
main(int argc, char** argv)
{
	// parse arguments
	warthog::util::param valid_args[] = 
	{
		{"core",  required_argument, 0, 2},
		{"input",  required_argument, 0, 2},
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

    // dump
    warthog::graph::write_xy(std::cout, g_xy);
    return 0;
}

