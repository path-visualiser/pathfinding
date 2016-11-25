#include "contraction.h"
#include "planar_graph.h"
#include <set>

int
main(int argv, char** argc)
{   
    if(argv < 4)
    {
        std::cout << "Usage: ./closure [gr file] [co file] [source id]"
            << " [optional exclude ids] " << std::endl;
        return 0;
    }

    warthog::graph::planar_graph g;
    g.load_dimacs(argc[1], argc[2]);

    std::set<uint32_t> closure;
    uint32_t node_id = atoi(argc[3]);

    if(argv > 4)
    {
        for(int i = 4; i < argv; i++)
        {
            std::cout << "ignore " << (atoi(argc[i])) << std::endl;
            closure.insert(atoi(argc[i]));
        }
    }
    warthog::ch::compute_closure(node_id, &g, &closure, 1);

    std::cout << "closure: \n";
    for(std::set<uint32_t>::iterator it = closure.begin(); 
            it != closure.end(); it++)
    {
        std::cout << (*it) << " ";
    }
    std::cout << std::endl;
    return 0;
}



