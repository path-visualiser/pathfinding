#include "dimacs_parser.h"
#include "graph.h"

#include <iostream>
#include <string>

warthog::graph::graph()
{
    filename_ = new std::string("");
    nodes_ = new std::vector<warthog::node>();
    edges_ = new std::vector<warthog::edge>();
    lat_lng_ = new std::vector<uint32_t>();
}

warthog::graph::~graph()
{
    delete nodes_;
    delete edges_;
    delete lat_lng_;
    delete [] filename_;
}

bool
warthog::graph::load_dimacs(char* gr_file, char* co_file)
{
    warthog::dimacs_parser dimacs(gr_file, co_file);
    return 1;
}

bool
warthog::graph::load_grid(char* file)
{
    std::cerr << "not yet!";
    return 1;
}


