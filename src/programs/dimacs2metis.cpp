#include "dimacs_parser.h"
#include <iostream>

void
help()
{
    std::cerr 
       << "Converts between the graph format used at the 9th DIMACS Implementation\n"
       << "Challenge and the graph format used by the METIS graph partitioning library.\n"
       << "Usage: ./dimacs2metis [dimacs .gr file]\n";
}

int 
main(int argc, char** argv)
{
    if(argc != 2)
    {
		help();
        exit(0);
    }

    warthog::dimacs_parser parser;
    parser.load_graph(argv[1]);
    parser.print_undirected_unweighted_metis(std::cout);
    return 0;
}

