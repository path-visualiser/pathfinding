#include "jpg_expansion_policy.h"
#include "corner_point_locator.h"
#include "jps.h"
#include "planar_graph.h"

const uint32_t JP_ID_MASK = (1 << 24) - 1;
warthog::jps::jpg_expansion_policy::jpg_expansion_policy(
        warthog::graph::planar_graph* g, warthog::gridmap* gm) 
    : expansion_policy(g->get_num_nodes() + 2) // +2 because we insert s & t
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

}

warthog::jps::jpg_expansion_policy::~jpg_expansion_policy()
{
    delete cpl_;
}

void
warthog::jps::jpg_expansion_policy::get_xy(
        warthog::search_node* n, int32_t& x, int32_t& y)
{
    uint32_t graph_id = n->get_id();
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
warthog::jps::jpg_expansion_policy::expand(
        warthog::search_node* current, warthog::problem_instance* problem)
{
    reset();
    problem_ = problem;

    uint32_t graph_id = current->get_id();

    // target has no outgoing arcs, so no point expanding it
    if(graph_id == t_graph_id_) { return; }

    // special case for the start node
    if(graph_id == s_graph_id_)
    {
        // start node is part of the graph
        if(graph_id < g_->get_num_nodes())
        {
            warthog::graph::node* n = g_->get_node(graph_id);
            for(warthog::graph::edge_iter eit = n->outgoing_begin(); 
                    eit != n->outgoing_end(); eit++)
            {
                add_neighbour(generate(eit->node_id_), eit->wt_);
            }
        }
        // start node is not part of the graph and was inserted
        else
        {
            for(uint32_t i = 0; i < s_outgoing_.size(); i++)
            {
                warthog::graph::edge e = s_outgoing_.at(i);
                add_neighbour(generate(e.node_id_), e.wt_);
            }
        }

    }
    // all other nodes, all of which are in the graph
    else
    {
        // compute the travel direction from the parent
        int32_t px, py, x, y;
        this->get_xy(current->get_parent(), px, py);
        this->get_xy(current, x, y);
        warthog::jps::direction dir_c = 
            warthog::jps::compute_direction((uint32_t)px, (uint32_t)py, 
                    (uint32_t)x, (uint32_t)y);
        
        // compute the set of directions in which we should 
        // generate successors
        uint32_t tiles;
        uint32_t grid_id = cpl_->get_map()->to_padded_id(x, y);
        cpl_->get_map()->get_neighbours(grid_id, (uint8_t*)&tiles);
        uint32_t succ_dirs = warthog::jps::compute_successors(dir_c, tiles);

        warthog::graph::node* n = g_->get_node(graph_id);
        for(warthog::graph::edge_iter eit = n->outgoing_begin(); 
                eit != n->outgoing_end(); eit++)
        {
            int32_t sx, sy;
            g_->get_xy(eit->node_id_, sx, sy);
            warthog::jps::direction s_dir = 
                warthog::jps::compute_direction(
                        (uint32_t)x, (uint32_t)y, 
                        (uint32_t)sx, (uint32_t)sy);
            if(s_dir & succ_dirs || eit->node_id_ == t_graph_id_)
            {
                add_neighbour(generate(eit->node_id_), eit->wt_);
            }
        }
    }
}

void
warthog::jps::jpg_expansion_policy::insert(
        uint32_t start_grid_id, uint32_t target_grid_id, 
        uint32_t& start_graph_id, uint32_t& target_graph_id)
{
    // NB: must insert target before start in case the target
    // can be reached without traveling through any jump point
    target_graph_id = insert_target(start_grid_id, target_grid_id);
    start_graph_id = insert_start(start_grid_id, target_grid_id);
}


uint32_t
warthog::jps::jpg_expansion_policy::insert_start(
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
            uint32_t corner_id = corner_neis[i] & JP_ID_MASK;
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
            s_outgoing_.push_back(warthog::graph::edge(graph_id, cost));
        }
        return s_graph_id_;
    }
    s_graph_id_ = myit->second;
    return s_graph_id_;
}

uint32_t
warthog::jps::jpg_expansion_policy::insert_target(
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
        uint32_t corner_id = corner_neis[i] & JP_ID_MASK;
        double cost = costs[i];
        uint32_t graph_id = (*id_map_.find(corner_id)).second;
        warthog::graph::node* nei = g_->get_node(graph_id);
        warthog::graph::edge_iter eit = 
            nei->add_outgoing(warthog::graph::edge(t_graph_id_, cost));
        if(eit != nei->outgoing_end())
        {
            //std::cout << "inserted target edge " << t_graph_id_ << " " << graph_id << "; cost " << cost << std::endl;
            t_incoming_.push_back(graph_id);
            t_incoming_iters_.push_back(eit);
        }
    }
    return t_graph_id_;
}

void
warthog::jps::jpg_expansion_policy::uninsert()
{
    if(s_grid_id_ != warthog::INF)
    {
        id_map_.erase(id_map_.find(s_grid_id_));
        s_outgoing_.clear();
        s_graph_id_ = warthog::INF;
        s_grid_id_ = warthog::INF;
    }

    if(t_graph_id_ > g_->get_num_nodes())
    {
        id_map_.erase(id_map_.find(t_grid_id_));
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
warthog::jps::jpg_expansion_policy::mem()
{
    return 
        expansion_policy::mem() +
        sizeof(*this) + 
        g_->mem() +
        cpl_->mem() + 
        id_map_.size()*sizeof(uint32_t);
}

