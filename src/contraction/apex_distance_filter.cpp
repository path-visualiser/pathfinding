#include "apex_distance_filter.h"
#include "ch_expansion_policy.h"
#include "down_distance_filter.h"
#include "flexible_astar.h"
#include "planar_graph.h"
#include "search_node.h"
#include "zero_heuristic.h"

warthog::apex_distance_filter::apex_distance_filter(
        const char* apexfile, warthog::graph::planar_graph* g, 
        std::vector<uint32_t>* rank)
{
    apex_dist_ = new std::vector<double>();
    ddist_ = new std::vector<double>();
    g_ = g;
    rank_ = rank;

    load_apex_distance_values(apexfile);
    start_id_ = 0;
    last_id_ = apex_dist_->size();
}

warthog::apex_distance_filter::apex_distance_filter(
        warthog::graph::planar_graph* g, std::vector<uint32_t>* rank)
{
    apex_dist_ = new std::vector<double>();
    ddist_ = new std::vector<double>();
    g_ = g;
    rank_ = rank;
    start_id_ = last_id_ = warthog::INF;
}

warthog::apex_distance_filter::~apex_distance_filter()
{
    delete apex_dist_;
    delete ddist_;
}

double
warthog::apex_distance_filter::get_apex_distance(uint32_t node_id)
{
    if(node_id < start_id_ || node_id > last_id_)
    {
        return warthog::INF;
    }
    return apex_dist_->at(node_id);
}

double
warthog::apex_distance_filter::get_down_distance(uint32_t node_id)
{
    if(node_id < start_id_ || node_id > last_id_)
    {
        return warthog::INF;
    }
    return ddist_->at(node_id);
}

bool
warthog::apex_distance_filter::filter(warthog::search_node* n)
{
    uint32_t node_id = n->get_id();
    if(node_id < start_id_ || node_id > last_id_)
    {
        return warthog::INF;
    }

    // try to prune up successors
    if(rank_->at(n->get_parent()->get_id()) < rank_->at(node_id))
    {
        if((n->get_g() + apex_dist_->at(node_id - start_id_))
                < n->get_f()) { return true; } 

        return false; 
    }

    // the optimal distance to the target cannot be further
    // than f(n). If g(n) + down-distance < f(n) then the
    // node n cannot possibly be on an optimal path to the target
    if((n->get_g() + ddist_->at(node_id - start_id_))
            < n->get_f()) { return true; } 

    // cannot prune this down successor
    return false;
}

void
warthog::apex_distance_filter::compute()
{
    if(g_ && rank_)
    {
        compute(0, g_->get_num_nodes());
    }
}

void
warthog::apex_distance_filter::compute(
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

    warthog::down_distance_filter ddf(g_, rank_);
    ddf.compute();

    std::cerr << "computing apex-distance labels\n";
    warthog::zero_heuristic heuristic;
    warthog::ch_expansion_policy expander(g_, rank_, false, warthog::ch::UP);

    warthog::flexible_astar<
        warthog::zero_heuristic,
        warthog::ch_expansion_policy> dijkstra(&heuristic, &expander);

    apex_dist_->clear();
    double max_apex = 0;
    std::function<void(warthog::search_node*)> fn = 
            [&ddf,&max_apex] (warthog::search_node* n) -> void
            {
               double apex_dist = n->get_g() +
               ddf.get_down_distance(n->get_id());
               if(apex_dist > max_apex)
               {
                    max_apex = apex_dist;
               }
            };

    for(uint32_t i = startid; i < endid; i++)
    {
        std::cerr << "\rprogress: " << (i+1) << " / " << endid;
        dijkstra.get_length(i, warthog::INF);
        max_apex = 0;
        dijkstra.apply_to_closed(fn);
        apex_dist_->push_back(max_apex);
        ddist_->push_back(ddf.get_down_distance(i));
    }
    std::cerr << "\nall done\n"<< std::endl;



}

void
warthog::apex_distance_filter::print(std::ostream& out)
{
    out << "id\tddist\tapex_dist"<<std::endl;
    for(uint32_t i = 0; i < apex_dist_->size(); i++)
    {
        out << i + start_id_ << "\t" << ddist_->at(i) << "\t" 
            << apex_dist_->at(i) << std::endl;
    }
}

void
warthog::apex_distance_filter::load_apex_distance_values(const char* filename)
{
   apex_dist_->clear();
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
        double apex_dist, ddist;
        ifs >> id >> ddist >> apex_dist;
        if(ifs.eof()) { break; }

        assert(id == apex_dist_->size());
        apex_dist_->push_back(apex_dist);
        ddist_->push_back(ddist);
        lines++;
    }
    ifs.close();
}
