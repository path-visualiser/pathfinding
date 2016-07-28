#include "ch_expansion_policy.h"
#include "down_distance_filter.h"
#include "flexible_astar.h"
#include "planar_graph.h"
#include "search_node.h"
#include "zero_heuristic.h"

#include <algorithm>

warthog::down_distance_filter::down_distance_filter(
        const char* ddfile, warthog::graph::planar_graph* g, 
        std::vector<uint32_t>* rank)
{
    ddist_ = new std::vector<double>();
    g_ = g;
    rank_ = rank;

    load_down_distance_values(ddfile);
    start_id_ = 0;
    last_id_ = ddist_->size();
}

warthog::down_distance_filter::down_distance_filter(
        warthog::graph::planar_graph* g, std::vector<uint32_t>* rank)
{
    ddist_ = new std::vector<double>();
    g_ = g;
    rank_ = rank;
    start_id_ = last_id_ = warthog::INF;
}

warthog::down_distance_filter::~down_distance_filter()
{
    delete ddist_;
}

double
warthog::down_distance_filter::get_down_distance(uint32_t node_id)
{
    if(node_id < start_id_ || node_id > last_id_)
    {
        return warthog::INF;
    }
    return ddist_->at(node_id);
}

bool
warthog::down_distance_filter::filter(warthog::search_node* n)
{
    uint32_t node_id = n->get_id();
    if(node_id < start_id_ || node_id > last_id_)
    {
        return warthog::INF;
    }

    if(rank_->at(n->get_parent()->get_id()) < rank_->at(node_id))
    {
        // cannot prune up successors
        return false; 
    }

    // the optimal distance to the target cannot be further
    // than f(n). If g(n) + down-distance < f(n) then the
    // node n cannot possibly be on an optimal path to the target
    if((n->get_g() + ddist_->at(node_id - start_id_))
            < n->get_f()) { return true; } 
    return false;
}

void
warthog::down_distance_filter::compute_down_distance()
{
    if(g_ && rank_)
    {
        compute_down_distance(0, g_->get_num_nodes());
    }
}

void
warthog::down_distance_filter::compute_down_distance(
        uint32_t startid, uint32_t endid)
{
    if(!g_ || !rank_) { return; } 

    if(startid > g_->get_num_nodes())
    { return; }
    start_id_ = startid;

    if(endid > g_->get_num_nodes())
    {
        endid = g_->get_num_nodes();
    }
    last_id_ = endid;

    std::cerr << "computing downdistance labels\n";
    warthog::zero_heuristic heuristic;
    warthog::ch_expansion_policy expander(g_, rank_, false, warthog::ch::DOWN);

    warthog::flexible_astar<
        warthog::zero_heuristic,
        warthog::ch_expansion_policy> dijkstra(&heuristic, &expander);

    ddist_->clear();
    for(uint32_t i = startid; i < endid; i++)
    {
        std::cerr << "\rprogress: " << (i+1) << " / " << endid;
        dijkstra.get_length(i, warthog::INF);
        ddist_->push_back(dijkstra.max_g_on_closed());
    }
    std::cerr << "\nall done\n"<< std::endl;
}

void
warthog::down_distance_filter::print(std::ostream& out)
{
    out << "id\tdown_dist"<<std::endl;
    for(uint32_t i = 0; i < ddist_->size(); i++)
    {
        out << i + start_id_ << "\t" << ddist_->at(i) << std::endl;
    }
}

void
warthog::down_distance_filter::load_down_distance_values(const char* filename)
{
   ddist_->clear();
   std::ifstream ifs(filename, std::ios_base::in);
    
    if(!ifs.good())
    {
        std::cerr << "\nerror trying to load down-distance values from file " 
            << filename << std::endl;
        ifs.close();
        return;
    }

    uint32_t lines = 0;
    while(true)
    {
        // skip header line
        while(ifs.peek() == 'i')
        {
            while(ifs.good() && ifs.get() != '\n');
            lines++;
        }

        uint32_t id;
        double ddist;
        ifs >> id >> ddist;
        if(ifs.eof()) { break; }

        assert(id == ddist_->size());
        ddist_->push_back(ddist);
        lines++;
    }
    ifs.close();
}
