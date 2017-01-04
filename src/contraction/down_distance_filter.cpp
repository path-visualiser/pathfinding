#include "ch_expansion_policy.h"
#include "down_distance_filter.h"
#include "flexible_astar.h"
#include "planar_graph.h"
#include "search_node.h"
#include "zero_heuristic.h"

#include <algorithm>

warthog::down_distance_filter::down_distance_filter(
        const char* ddfile, warthog::graph::planar_graph* g)
{
    ddist_ = new std::vector<double>();
    g_ = g;

    load_labels(ddfile);
    start_id_ = 0;
    last_id_ = ddist_->size();
}

warthog::down_distance_filter::down_distance_filter(
        warthog::graph::planar_graph* g)
{
    ddist_ = new std::vector<double>();
    g_ = g;
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

// @param n: the node we are trying to prune
// @param edge_id: used to indicate that @param n is the nth 
// child of its parent
//
// NB: this function assumes @param n is reached via a down-edge 
bool
warthog::down_distance_filter::filter(
        warthog::search_node* n, uint32_t edge_id)

{
    uint32_t node_id = n->get_id();
    assert((last_id_ - node_id) >= start_id_);

    // the optimal distance to the target cannot be further
    // than f(n). If g(n) + down-distance < f(n) then the
    // node n cannot possibly be on an optimal path to the target
    if((n->get_g() + ddist_->at(node_id - start_id_))
            < n->get_f()) { return true; } 
    return false;
}

void
warthog::down_distance_filter::compute(std::vector<uint32_t>* rank)
{
    if(g_ && rank)
    {
        compute(0, g_->get_num_nodes()-1, rank);
    }
}

void
warthog::down_distance_filter::compute(
        uint32_t startid, uint32_t endid, std::vector<uint32_t>* rank)
{
    if(!g_ || !rank) { return; } 

    if(startid > g_->get_num_nodes())
    { return; }

    if(endid < startid || endid >= g_->get_num_nodes())
    { return; }

    start_id_ = startid;
    last_id_ = endid;

    std::cerr << "computing downdistance labels\n";
    warthog::zero_heuristic heuristic;
    warthog::ch_expansion_policy expander(g_, rank, false, warthog::ch::DOWN);

    warthog::flexible_astar<
        warthog::zero_heuristic,
        warthog::ch_expansion_policy> dijkstra(&heuristic, &expander);


    ddist_->clear();
    for(uint32_t i = start_id_; i <= last_id_; i++)
    {
        std::cerr << "\rprocessing node " << i << ". continues to " << last_id_;
        dijkstra.get_length(warthog::problem_instance(i, warthog::INF));

        // scan the closed list and get the max g-value
        double max = 0;
        std::function<void(warthog::search_node*)> fn_max_g =
            [&max](warthog::search_node* n) -> void
            {
                assert(n);
                if(n->get_g() > max)
                {
                    max = n->get_g();
                }
            };
        dijkstra.apply_to_closed(fn_max_g);
        ddist_->push_back(max);
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

bool
warthog::down_distance_filter::load_labels(const char* filename)
{
   ddist_->clear();
   std::ifstream ifs(filename, std::ios_base::in);
    
    if(!ifs.good())
    {
        std::cerr << "\nerror trying to load down-distance values from file " 
            << filename << std::endl;
        ifs.close();
        return false;
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
    return true;
}
