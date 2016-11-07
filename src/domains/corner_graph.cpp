#include "corner_graph.h"
#include "corner_point_locator.h"
#include "jps.h"

typedef std::unordered_map<uint32_t, uint32_t>::iterator id_map_iter;

warthog::graph::corner_graph::corner_graph(
        warthog::gridmap* gm, 
        warthog::graph::planar_graph* g)
{
    g_ = g;
    cpl_ = new warthog::jps::corner_point_locator(gm);

    // we construct a hash map to quickly convert from the grid identifier
    // of a node to the corresponding graph identifier of the same node
    for(uint32_t graph_id = 0; graph_id < g_->get_num_nodes(); graph_id++)
    {
        int32_t x, y;
        g_->get_xy(graph_id, x, y);
        uint32_t grid_id = gm->to_padded_id(y * gm->header_width() + x);
        id_map_.insert(
                std::pair<uint32_t, uint32_t>(grid_id, graph_id));
    }
    s_grid_id_ = t_grid_id_ = warthog::INF;
    s_graph_id_ = t_graph_id_ = warthog::INF;
    s_ = new warthog::graph::node();
    t_ = new warthog::graph::node();
}

warthog::graph::corner_graph::~corner_graph()
{
    delete t_;
    delete s_;
    delete cpl_;
}

void
warthog::graph::corner_graph::get_xy(uint32_t graph_id, int32_t& x, int32_t& y)
{
    uint32_t grid_x, grid_y;
    if(graph_id < g_->get_num_nodes())
    {
        g_->get_xy(graph_id, x, y);
        return;
    }
    else if(graph_id == s_graph_id_)
    {
        cpl_->get_map()->to_unpadded_xy(s_grid_id_, grid_x, grid_y);
    }
    else if(graph_id  == t_graph_id_)
    {
        cpl_->get_map()->to_unpadded_xy(t_grid_id_, grid_x, grid_y);
    }
    else
    {
        grid_x = grid_y = warthog::INF;
    }
    x = (int32_t)grid_x;
    y = (int32_t)grid_y;
}


void
warthog::graph::corner_graph::insert(
        uint32_t start_grid_id, uint32_t target_grid_id)
{
    // clear any previous insertion data
    uninsert();

    // NB: order is important here 
    // (target before start, lest we jump over the target)
    insert_target(start_grid_id, target_grid_id);
    insert_start(start_grid_id, target_grid_id);
}

void
warthog::graph::corner_graph::insert_start(
        uint32_t start_grid_id, uint32_t target_grid_id)
{
    id_map_iter myit = id_map_.find(start_grid_id);
    if(myit == id_map_.end())
    {
        s_grid_id_ = start_grid_id;
        s_graph_id_ = g_->get_num_nodes();
        id_map_.insert(std::pair<uint32_t, uint32_t>(s_grid_id_, s_graph_id_));

        // compute jump point successors for the start node
        std::vector<uint32_t> corner_neis;
        std::vector<double> costs;
        for(uint32_t i = 0; i < 8; i++)
        {
            warthog::jps::direction d = (warthog::jps::direction)(1 << i);
            cpl_->jump(d, start_grid_id, target_grid_id, corner_neis, costs);
        }

        for(uint32_t i = 0; i < corner_neis.size(); i++)
        {
            uint32_t corner_id = corner_neis[i] & warthog::jps::ID_MASK;
            double cost = costs[i];
            uint32_t graph_id;
            if(corner_id == target_grid_id)
            {
                graph_id = t_graph_id_;
            }
            else
            {
               id_map_iter corner_iter = id_map_.find(corner_id);
               assert(corner_iter != id_map_.end());
               graph_id = corner_iter->second;
            }
            s_->add_outgoing(warthog::graph::edge(graph_id, cost));
        }
    }
    else
    {
        s_graph_id_ = myit->second;
        s_grid_id_ = start_grid_id;
    }
}

void
warthog::graph::corner_graph::insert_target(
        uint32_t start_grid_id, uint32_t target_grid_id)
{
    t_grid_id_ = target_grid_id;

    id_map_iter myit = id_map_.find(target_grid_id);
    if(myit == id_map_.end())
    {
        t_graph_id_ = g_->get_num_nodes() + 1;
        id_map_.insert(std::pair<uint32_t, uint32_t>(t_grid_id_, t_graph_id_));
    }
    else
    {
        t_graph_id_ = myit->second;
    }

    std::vector<uint32_t> corner_neis;
    std::vector<double> costs;
    for(uint32_t i = 0; i < 8; i++)
    {
        warthog::jps::direction d = (warthog::jps::direction)(1 << i);
        cpl_->jump(d, target_grid_id, warthog::INF, corner_neis, costs);
    }

    for(uint32_t i = 0; i < corner_neis.size(); i++)
    {
        uint32_t corner_id = corner_neis[i] & warthog::jps::ID_MASK;
        double cost = costs[i];
        uint32_t graph_id = (*id_map_.find(corner_id)).second;
        warthog::graph::node* nei = g_->get_node(graph_id);
        warthog::graph::edge_iter eit = 
            nei->add_outgoing(warthog::graph::edge(t_graph_id_, cost));
        if(eit != nei->outgoing_end())
        {
            //std::cout 
            //<< "inserted target edge " << t_graph_id_ << " " 
            //<< graph_id << "; cost " << cost << std::endl;
            t_incoming_.push_back(graph_id);
            t_incoming_iters_.push_back(eit);
        }
    }
}

void
warthog::graph::corner_graph::uninsert()
{
    if(s_graph_id_ != warthog::INF)
    {
        if(s_graph_id_ == g_->get_num_nodes())
        {
            id_map_.erase(id_map_.find(s_grid_id_));
            s_->clear();
        }
        s_graph_id_ = warthog::INF;
        s_grid_id_ = warthog::INF;
    }

    if(t_graph_id_ != warthog::INF)
    {
        if(t_graph_id_ == (g_->get_num_nodes() + 1))
        {
            id_map_.erase(id_map_.find(t_grid_id_));
        }
        t_graph_id_ = warthog::INF;
        t_grid_id_ = warthog::INF;
    }
    for(uint32_t i = t_incoming_.size()-1; i < t_incoming_.size() ; i--)
    {
        warthog::graph::node* nei = g_->get_node(t_incoming_.at(i));
        nei->del_outgoing(t_incoming_iters_.at(i));
    }
    t_incoming_iters_.clear();
    t_incoming_.clear();
}

size_t
warthog::graph::corner_graph::mem()
{
    return g_->mem() + cpl_->mem() + sizeof(*this);
}

