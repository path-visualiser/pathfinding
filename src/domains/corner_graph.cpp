#include "corner_graph.h"
#include "corner_point_locator.h"
#include "jps.h"

typedef std::unordered_map<uint32_t, uint32_t>::iterator id_map_iter;

warthog::graph::corner_graph::corner_graph(std::shared_ptr<warthog::gridmap> gm)
{
    g_ = std::shared_ptr<warthog::graph::planar_graph>(
            new warthog::graph::planar_graph());
    cpl_ = new warthog::jps::corner_point_locator(gm.get());

    e_lab_index_ = new std::vector<label_index>(g_->get_num_nodes());


    // do the thing
    construct();

    // we now insert two dummy nodes into the graph which will be
    // used whenever the start and target location are not corner points
    S_DUMMY_ID = this->add_node(0, 0); 
    T_DUMMY_ID = this->add_node(0, 0);

    s_grid_id_ = t_grid_id_ = warthog::INF;
    s_graph_id_ = t_graph_id_ = warthog::INF;
}

warthog::graph::corner_graph::corner_graph(
        std::shared_ptr<warthog::gridmap> gm, 
        std::shared_ptr<warthog::graph::planar_graph> g)
{
    g_ = g;
    cpl_ = new warthog::jps::corner_point_locator(gm.get());

    s_graph_id_ = t_graph_id_ = warthog::INF;
    e_lab_index_ = new std::vector<label_index>(g_->get_num_nodes());

    // we now insert two dummy nodes into the graph which will be
    // used whenever the start and target location are not corner points
    S_DUMMY_ID = this->add_node(0, 0); 
    T_DUMMY_ID = this->add_node(0, 0);

    s_grid_id_ = t_grid_id_ = warthog::INF;
    s_graph_id_ = t_graph_id_ = warthog::INF;

    // map each graph node id to its corresponding grid id
    for(uint32_t graph_id = 0; graph_id < g_->get_num_nodes(); graph_id++)
    {
        warthog::graph::node* n = g_->get_node(graph_id);
        if(!n) { continue; } // skip padding at the start of the graph array

        int32_t x, y;
        g_->get_xy(graph_id, x, y);
        assert(x != INT32_MAX && y != INT32_MAX);
        uint32_t gm_id = gm->to_padded_id(x, y);
        id_map_.insert(std::pair<uint32_t, uint32_t>(gm_id, graph_id));
    }
}

warthog::graph::corner_graph::~corner_graph()
{
    delete e_lab_index_;
    delete cpl_;
}

void
warthog::graph::corner_graph::construct()
{
    warthog::gridmap* gm = cpl_->get_map();
    warthog::gridmap* corner_map = cpl_->get_corner_map();
    uint32_t mapwidth = gm->header_width();
    uint32_t mapheight = gm->header_height();
    for(uint32_t index = 0; index < mapheight*mapwidth; index++)
    {
        uint32_t gm_id = gm->to_padded_id(index);
        if(!corner_map->get_label(gm_id)) { continue; } 

        // each corner point in the grid becomes a node in 
        // the corresponding planar graph
        uint32_t gm_x, gm_y;
        gm->to_unpadded_xy(gm_id, gm_x, gm_y);

        uint32_t from_id;
        std::unordered_map<uint32_t, uint32_t>::iterator it_from_id;
        it_from_id = id_map_.find(gm_id);
        if(it_from_id == id_map_.end())
        {
            from_id = this->add_node(gm_x, gm_y);
            id_map_.insert(std::pair<uint32_t, uint32_t>(gm_id, from_id));
        }
        else { from_id = it_from_id->second; }

        // identify other corner points by scanning the map and adding them 
        // to the graph as necessary. also added are edges between corner 
        // points to their immediate corner-point successors
        for(uint32_t i = 0; i < 8; i++)
        {
            warthog::jps::direction d = (warthog::jps::direction)(1 << i);
            std::vector<uint32_t> jpoints;
            std::vector<double> jcosts;
            cpl_->jump(d, gm_id, warthog::INF, jpoints, jcosts);
            for(uint32_t idx = 0; idx < jpoints.size(); idx++)
            {
                //uint8_t pdir = ((uint8_t*)&jpoints[idx])[3];
                uint32_t jp_id = jpoints[idx] & ((1 << 24) - 1);
                std::unordered_map<uint32_t, uint32_t>::iterator it_to_id; 
                it_to_id = id_map_.find(jp_id);

                uint32_t to_id;
                if(it_to_id == id_map_.end())
                {
                    uint32_t jp_x, jp_y;
                    gm->to_unpadded_xy(jp_id, jp_x, jp_y);
                    to_id = this->add_node(jp_x, jp_y);
                    id_map_.insert(std::pair<uint32_t, uint32_t>(jp_id, to_id));
                }
                else { to_id = it_to_id->second; }

                uint32_t label_id = __builtin_ffs(d)-1;
                this->add_labeled_outgoing_edge(
                        from_id, warthog::graph::edge(to_id, jcosts[idx]), 
                        label_id);
            }
        }
    }
}

void
warthog::graph::corner_graph::get_xy(uint32_t graph_id, int32_t& x, int32_t& y)
{
    g_->get_xy(graph_id, x, y);
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
        s_graph_id_ = S_DUMMY_ID;

        uint32_t sx, sy;
        warthog::gridmap* gm = cpl_->get_map();
        gm->to_unpadded_xy(start_grid_id, sx, sy);
        g_->set_xy(s_graph_id_, (int32_t)sx, (int32_t)sy);

        id_map_.insert(std::pair<uint32_t, uint32_t>(s_grid_id_, s_graph_id_));

        // compute jump point successors for the start node
        std::vector<uint32_t> corner_neis;
        std::vector<double> costs;
        for(uint32_t i = 0; i < 8; i++)
        {
            warthog::jps::direction d = (warthog::jps::direction)(1 << i);
            cpl_->jump(d, start_grid_id, target_grid_id, corner_neis, costs);
        }

        warthog::graph::node* s = g_->get_node(s_graph_id_);
        for(uint32_t i = 0; i < corner_neis.size(); i++)
        {
            warthog::jps::direction d = 
                (warthog::jps::direction)(corner_neis[i] >> 24);
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
            s->add_outgoing(warthog::graph::edge(graph_id, cost, d));
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
        t_graph_id_ = T_DUMMY_ID;
        uint32_t tx, ty;
        warthog::gridmap* gm = cpl_->get_map();
        gm->to_unpadded_xy(target_grid_id, tx, ty);
        g_->set_xy(t_graph_id_, (int32_t)tx, (int32_t)ty);
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

        // goal is always the last node in the edges array of nei
        warthog::graph::node* nei = g_->get_node(graph_id);
        warthog::graph::edge_iter eit = 
            nei->add_outgoing(warthog::graph::edge(t_graph_id_, cost));
        t_incoming_.push_back(graph_id);
        t_incoming_iters_.push_back(eit);
    }
}

void
warthog::graph::corner_graph::uninsert()
{
    if(s_graph_id_ == S_DUMMY_ID)
    {
        id_map_.erase(id_map_.find(s_grid_id_));
        g_->get_node(s_graph_id_)->clear();
    }
    s_graph_id_ = warthog::INF;
    s_grid_id_ = warthog::INF;

    if(t_graph_id_ == T_DUMMY_ID)
    {
        id_map_.erase(id_map_.find(t_grid_id_));
    }
    t_graph_id_ = warthog::INF;
    t_grid_id_ = warthog::INF;

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
    size_t mem_sz = g_->mem() + cpl_->mem() + sizeof(*this);
    if(e_lab_index_)
    {
        mem_sz += sizeof(label_index)*e_lab_index_->capacity();
    }
    return mem_sz;  
}

warthog::graph::edge_iter
warthog::graph::corner_graph::add_labeled_outgoing_edge(
        uint32_t node_id, warthog::graph::edge e, uint32_t label_id)
{
    assert(label_id < NUM_LABELS);
    int32_t head = e_lab_index_->at(node_id).head_[label_id];
    
    // insert the edge into the neighbours array at 
    // the position indicated by head
    warthog::graph::node* n = g_->get_node(node_id);
    warthog::graph::edge_iter ret = n->insert_outgoing(e, head);

    // increment subsequent indexes to reflect the insertion
    for(uint32_t i = label_id+1; i < NUM_LABELS; i++)
    {
        assert(e_lab_index_->at(node_id).head_[i] >= head);
        e_lab_index_->at(node_id).head_[i]++;
    }
    return ret;
}

