#include "arclabels.h"
#include "ch_expansion_policy.h"
#include "bb_filter.h"
#include "flexible_astar.h"
#include "graph_expansion_policy.h"
#include "planar_graph.h"
#include "problem_instance.h"
#include "search_node.h"
#include "zero_heuristic.h"

#include <algorithm>

warthog::bb_filter::bb_filter()
{
    tx_ = ty_ = warthog::INF;
    g_ = 0;
}

warthog::bb_filter::~bb_filter()
{
}

void
warthog::bb_filter::print(std::ostream& out)
{
    if(g_) { warthog::arclabels::bb_print(labels_, g_, out); }
}

bool
warthog::bb_filter::load_labels(
        const char* filename, warthog::graph::planar_graph* g)
{
    std::ifstream ifs(filename, std::ios_base::in);
    if(!ifs.good())
    {
        std::cerr << "\nerror trying to load bbox values from file " 
            << filename << std::endl;
        ifs.close();
        return false;
    }

    uint32_t lines = 0;
    labels_.resize(g->get_num_nodes());
    // skip comment lines
    while(ifs.peek() == '#')
    {
        while(ifs.good() && ifs.get() != '\n');
        lines++;
    }
    
    bool bad_header = false;
    std::string token;

    uint32_t start_id_, last_id_;
    ifs >> token;
    if(token == "firstnode") { ifs >> start_id_; }
    else { bad_header = true; }

    ifs >> token;
    if(token == "lastnode") { ifs >> last_id_; }
    else { bad_header = true; }

    if(bad_header)
    {
        std::cerr << " invalid header in file " << filename << std::endl;
        return false;
    }

    if(!(start_id_ < g->get_num_nodes() && 
        last_id_ < g->get_num_nodes() && 
        start_id_ < last_id_)) 
    {
        std::cerr << " invalid header values in file " << filename << "\n";
        return false;
    }
    lines++;

    warthog::geom::rectangle dummy;
    for(uint32_t i = start_id_; i <= last_id_; i++)
    {
        warthog::graph::node* n = g->get_node(i);
        for(uint32_t j = 0; j < n->out_degree(); j++)
        {
            if(ifs.eof())
            {
                std::cerr << "unexpected eof in file " << filename << "\n";
                return false;
            }
            if(!ifs.good())
            {
                std::cerr << "unexpected read error while loading " 
                    << filename << "\n";
                return false;
            }

            int32_t x1, y1, x2, y2;
            ifs >> x1 >> y1 >> x2 >> y2;
            warthog::geom::rectangle rect = warthog::geom::rectangle(x1, y1, x2, y2);
            if(rect != dummy && !rect.is_valid())
            {
                std::cerr << "err; invalid label on line " << lines <<"\n";
                labels_.clear();
                return false;
            }
            labels_.at(i).push_back(rect);
            lines++;
        }
    }
    ifs.close();
    g_ = g;
    return true;
}

