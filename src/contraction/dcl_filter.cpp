#include "ch_expansion_policy.h"
#include "dcl_filter.h"
#include "euclidean_heuristic.h"
#include "flexible_astar.h"
#include "planar_graph.h"
#include "problem_instance.h"
#include "search_node.h"
#include "zero_heuristic.h"

#include <algorithm>

warthog::dcl_filter::dcl_filter(
        const char* ddfile, warthog::graph::planar_graph* g, 
        std::vector<uint32_t>* rank)
{
    g_ = g;
    rank_ = rank;

    load_dcl_file(ddfile);
    start_id_ = 0;
    tx = ty = warthog::INF;
    h_ = new warthog::euclidean_heuristic(g);
}

warthog::dcl_filter::dcl_filter(
        warthog::graph::planar_graph* g, std::vector<uint32_t>* rank)
{
    g_ = g;
    rank_ = rank;
    start_id_ = last_id_ = warthog::INF;
    tx = ty = warthog::INF;
}

warthog::dcl_filter::~dcl_filter()
{
    delete h_;
}

bool 
warthog::dcl_filter::filter(warthog::search_node* n)
{
    return filter(n->get_id(), 0, n->get_g());
}

bool
warthog::dcl_filter::filter(uint32_t node_id, uint32_t edge_id, double g_value)
{
    return edgebb_filter(node_id, edge_id);
}

// prune a node based on the maximum down distance
bool 
warthog::dcl_filter::dd_filter(uint32_t node_id, uint32_t edge_id, double g_value)
{
    int32_t nx, ny;
    g_->get_xy(node_id, nx, ny);
    double fvalue = g_value + h_->h(nx, ny, tx, ty);
    if(fvalue < (g_value + labels_.at(node_id).ddist_))
    { 
        return false; 
    }
    return true;
}

// prune any node whose down-bb does not contain the target
bool
warthog::dcl_filter::bb_filter(uint32_t node_id, uint32_t edge_id)
{
    if(labels_.at(node_id).bbox_.contains(tx, ty))
    {
        return false;
    }
    return true;
}

bool
warthog::dcl_filter::edgebb_filter(uint32_t node_id, uint32_t edge_id)
{
    if(labels_.at(node_id).edgebb_.at(edge_id).contains(tx, ty))
    {
        return false;
    }
    return true;
}

// prune a node by comparing the contraction level of the target node 
// against the contraction level range created from the down closure
bool
warthog::dcl_filter::lvl_filter(uint32_t node_id, uint32_t succ_num)
{
    if(target_lvl >= labels_.at(node_id).min_lvl_ && 
       target_lvl <= rank_->at(node_id))
    {
        return false;
    }
    return true;
}

void
warthog::dcl_filter::compute()
{
    if(!g_ || !rank_) { return; } 

    start_id_ = 0;
    last_id_ = g_->get_num_nodes() - 1;

    std::cerr << "computing down-closure bbox labels\n";
    labels_.resize(rank_->size());

    std::function<void(uint32_t)> top_down;
    top_down = [&top_down, this](uint32_t node_id) 
        -> void
    {
        warthog::graph::node* n = g_->get_node(node_id);
        if(labels_.at(node_id).bbox_.is_valid()) { return; }

        dcl_label& label = labels_.at(node_id);
        label.min_lvl_ = rank_->at(node_id);

        // recursively compute a bounding box that contains
        // all nodes in the down-closure of the current node
        for(warthog::graph::edge_iter it = n->outgoing_begin(); 
                it != n->outgoing_end(); it++)
        {
            uint32_t nei_id = (*it).node_id_;
            if(rank_->at(nei_id) < rank_->at(node_id))
            {
                if(!labels_.at(nei_id).bbox_.is_valid()) { top_down(nei_id); }
                label.bbox_.grow(labels_.at(nei_id).bbox_);
                label.ddist_ = std::max<double>(
                        label.ddist_, labels_.at(nei_id).ddist_ + (*it).wt_);
                label.min_lvl_ = std::min<uint32_t>(
                        label.min_lvl_, labels_.at(nei_id).min_lvl_);
                label.edgebb_.push_back(labels_.at(nei_id).bbox_);
            }
        }
        int32_t x, y;
        g_->get_xy(node_id, x, y);
        label.bbox_.grow(x, y);
    };

    for(uint32_t i = 0; i < g_->get_num_nodes(); i++)
    {
        std::cerr << "progress " << (i+1) << " / " 
            << g_->get_num_nodes() << "\r";
        top_down(i);
    }
    std::cerr << std::endl;
}

void
warthog::dcl_filter::print(std::ostream& out)
{
    out << "node_id\tmin_lvl\tddist\tx1\ty1\tx2\ty2"<<std::endl;
    for(uint32_t i = 0; i <labels_.size(); i++)
    {
        dcl_label& label = labels_.at(i);
        out << i << "\t"
            << label.min_lvl_ << "\t" << label.ddist_ << "\t"
            << label.bbox_.x1 << "\t" << label.bbox_.y1 << "\t"
            << label.bbox_.x2 << "\t" << label.bbox_.y2 << "\t"
            << std::endl;
    }
}

void
warthog::dcl_filter::load_dcl_file(const char* filename)
{
   std::ifstream ifs(filename, std::ios_base::in);
    
    if(!ifs.good())
    {
        std::cerr << "\nerror trying to load down-distance values from file " 
            << filename << std::endl;
        ifs.close();
        return;
    }

    uint32_t lines = 0;
    labels_.resize(g_->get_num_nodes());
    while(true)
    {
        // skip header line
        while(ifs.peek() == 'n')
        {
            while(ifs.good() && ifs.get() != '\n');
            lines++;
        }

        uint32_t from_id, min_lvl;
        double ddist;
        int32_t x1, y1, x2, y2;
        ifs >> from_id >> min_lvl >> ddist >> x1 >> y1 >> x2 >> y2;
        if(ifs.eof()) { break; }

        assert(from_id < bbox_.size());
        labels_.at(from_id).bbox_ = warthog::geom::rectangle(x1, y1, x2, y2);
        labels_.at(from_id).ddist_ = ddist;
        labels_.at(from_id).min_lvl_ = min_lvl;
        lines++;
    }
    ifs.close();

    for(uint32_t node_id = 0; node_id < g_->get_num_nodes(); node_id++)
    {
        warthog::graph::node* n = g_->get_node(node_id);
        for(warthog::graph::edge_iter it = n->outgoing_begin(); 
                it != n->outgoing_end(); it++)
        {
            uint32_t nei_id = (*it).node_id_;
            if(rank_->at(nei_id) < rank_->at(node_id))
            {
                labels_.at(node_id).edgebb_.push_back(labels_.at(nei_id).bbox_);
            }
        }
    }
}

void
warthog::dcl_filter::init(warthog::problem_instance* instance)
{
    assert(instance);
    g_->get_xy(instance->get_goal(), tx, ty);
    target_lvl = rank_->at(instance->get_goal());
}
