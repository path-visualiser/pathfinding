#include "ch_expansion_policy.h"
#include "bb_filter.h"
#include "flexible_astar.h"
#include "graph_expansion_policy.h"
#include "planar_graph.h"
#include "problem_instance.h"
#include "search_node.h"
#include "zero_heuristic.h"

#include <algorithm>

warthog::bb_filter::bb_filter(const char* ddfile, 
        warthog::graph::planar_graph* g)

{
    g_ = g;

    load_labels(ddfile);
    start_id_ = 0;
    tx = ty = warthog::INF;
}

warthog::bb_filter::bb_filter(warthog::graph::planar_graph* g)
{
    g_ = g;
    start_id_ = last_id_ = warthog::INF;
    tx = ty = warthog::INF;
}

warthog::bb_filter::~bb_filter()
{
}

bool 
warthog::bb_filter::filter(warthog::search_node* n, uint32_t edge_id)
{
    return filter__(n->get_id(), edge_id);
}

bool
warthog::bb_filter::filter__(uint32_t node_id, uint32_t edge_id)
{
    // prune any node whose bounding-box 
    // does not contain the target
    if(labels_.at(node_id).at(edge_id).contains(tx, ty))
    {
        return false;
    }
    return true;
}

void
warthog::bb_filter::compute()
{
    compute(0, g_->get_num_nodes());
}

void
warthog::bb_filter::compute(uint32_t startid, uint32_t endid)
{
    if(!g_) { return; } 

    if(startid >= g_->get_num_nodes())
    { return; }
    start_id_ = startid;

    if(endid >= g_->get_num_nodes())
    {
        endid = g_->get_num_nodes()-1;
    }
    last_id_ = endid;

    std::cerr << "computing bbox labels\n";
    warthog::zero_heuristic heuristic;
    warthog::graph_expansion_policy expander(g_);

    warthog::flexible_astar<
        warthog::zero_heuristic,
        warthog::graph_expansion_policy> dijkstra(&heuristic, &expander);

    // need to keep track of the first edge on the way to the current node
    // (the solution is a bit hacky as we break the chain of backpointers 
    // to achieve this; it doesn't matter, we don't care about the path)
    std::function<void(warthog::search_node*)> relax_fn = 
            [] (warthog::search_node* n) -> void
            {
                // the start node and its children don't need their 
                // parent pointers updated. for all other nodes the
                // grandparent becomes the parent
                if(n->get_parent()->get_parent() != 0)
                {
                    if(n->get_parent()->get_parent()->get_parent() != 0)
                    {
                        n->set_parent(n->get_parent()->get_parent());
                    }
                }
            };
    dijkstra.apply_on_relax(relax_fn);

    for(uint32_t i = 0; i < ((last_id_ - start_id_)+1); i++)
    {
        std::cerr << "\rprogress: " << (i+start_id_) << " / " << endid;
        uint32_t source_id = i + start_id_;
        dijkstra.get_length(source_id, warthog::INF);
        
        // create a rectangle label for every outgoing edge
        warthog::graph::node* source = g_->get_node(source_id);
        labels_.push_back(
                std::vector<warthog::geom::rectangle>(
                    source->out_degree()));

        // we need an easy way to convert between the ids of nodes
        // adjacent to the source and their corresponding edge index
        std::unordered_map<uint32_t, uint32_t> idmap;
        uint32_t edge_idx = 0;
        for(warthog::graph::edge_iter it = source->outgoing_begin(); 
                it != source->outgoing_end(); it++)
        {
            idmap.insert(
                    std::pair<uint32_t, uint32_t>((*it).node_id_, edge_idx));
            edge_idx++;
        }

        // compute the extent of each rectangle bounding box
        std::function<void(warthog::search_node*)> bbox_fn = 
                [this, &source_id, &idmap](warthog::search_node* n) -> void
                {
                    // skip the source (it doesn't belong to any rectangle)
                    assert(n);
                    if(n->get_id() == source_id) { return; } 
                    assert(n->get_parent());

                    // grow the rectangle associated with each outgoing
                    // edge of the source node
                    std::unordered_map<uint32_t, uint32_t>::iterator it;
                    if(n->get_parent()->get_parent() == 0)
                    // special case for the successors of the source
                    { it = idmap.find(n->get_id()); }
                    else
                    // all the other nodes
                    { it = idmap.find(n->get_parent()->get_id()); }



                    if(it == idmap.end())
                    {
                        std::cerr << "err; invalid first edge; "
                             << " n id: "<< n->get_id() << " parent id: " 
                             << n->get_parent()->get_id() << std::endl;
                        exit(0);
                    }
                    int32_t x, y;
                    this->g_->get_xy(n->get_id(), x, y);
                    assert(x != warthog::INF && y != warthog::INF);
                    this->labels_.back().at((*it).second).grow(x, y);
                    assert(this->labels_.back().at((*it).second).is_valid());
                };
        dijkstra.apply_to_closed(bbox_fn);
    }
    std::cerr << "\nall done\n"<< std::endl;
}

void
warthog::bb_filter::compute_ch(std::vector<uint32_t>* rank)
{
    compute_ch(0, g_->get_num_nodes(), rank);
}

void
warthog::bb_filter::compute_ch(
    uint32_t startid, uint32_t endid, std::vector<uint32_t>* rank)
{
    if(!g_ || !rank) { return; } 

    if(startid >= g_->get_num_nodes())
    { return; }
    start_id_ = startid;

    if(endid >= g_->get_num_nodes())
    {
        endid = g_->get_num_nodes()-1;
    }
    last_id_ = endid;

    std::cerr << "computing bbox labels\n";
    warthog::zero_heuristic heuristic;
    warthog::ch_expansion_policy expander(g_, rank, false, warthog::ch::DOWN);

    warthog::flexible_astar<
        warthog::zero_heuristic,
        warthog::ch_expansion_policy> dijkstra(&heuristic, &expander);

    // need to keep track of the first edge on the way to the current node
    // (the solution is a bit hacky as we break the chain of backpointers 
    // to achieve this; it doesn't matter, we don't care about the path)
    std::function<void(warthog::search_node*)> relax_fn = 
            [] (warthog::search_node* n) -> void
            {
                // the start node and its children don't need their 
                // parent pointers updated. for all other nodes the
                // grandparent becomes the parent
                if(n->get_parent()->get_parent() != 0)
                {
                    if(n->get_parent()->get_parent()->get_parent() != 0)
                    {
                        n->set_parent(n->get_parent()->get_parent());
                    }
                }
            };
    dijkstra.apply_on_relax(relax_fn);

    for(uint32_t i = 0; i < ((last_id_ - start_id_)+1); i++)
    {
        std::cerr << "\rprogress: " << (i+start_id_) << " / " << endid;
        uint32_t source_id = i + start_id_;
        dijkstra.get_length(source_id, warthog::INF);
        
        // create a rectangle label for every outgoing edge
        warthog::graph::node* source = g_->get_node(source_id);
        labels_.push_back(
                std::vector<warthog::geom::rectangle>(
                    source->out_degree()));

        // we need an easy way to convert between the ids of nodes
        // adjacent to the source and their corresponding edge index
        std::unordered_map<uint32_t, uint32_t> idmap;
        uint32_t edge_idx = 0;
        for(warthog::graph::edge_iter it = source->outgoing_begin(); 
                it != source->outgoing_end(); it++)
        {
            idmap.insert(
                    std::pair<uint32_t, uint32_t>((*it).node_id_, edge_idx));
            edge_idx++;
        }

        // compute the extent of each rectangle bounding box
        std::function<void(warthog::search_node*)> bbox_fn = 
                [this, &source_id, &idmap](warthog::search_node* n) -> void
                {
                    // skip the source (it doesn't belong to any rectangle)
                    assert(n);
                    if(n->get_id() == source_id) { return; } 
                    assert(n->get_parent());

                    // grow the rectangle associated with each outgoing
                    // edge of the source node
                    std::unordered_map<uint32_t, uint32_t>::iterator it;
                    if(n->get_parent()->get_parent() == 0)
                    // special case for the successors of the source
                    { it = idmap.find(n->get_id()); }
                    else
                    // all the other nodes
                    { it = idmap.find(n->get_parent()->get_id()); }



                    if(it == idmap.end())
                    {
                        std::cerr << "err; invalid first edge; "
                             << " n id: "<< n->get_id() << " parent id: " 
                             << n->get_parent()->get_id() << std::endl;
                        exit(0);
                    }
                    int32_t x, y;
                    this->g_->get_xy(n->get_id(), x, y);
                    assert(x != warthog::INF && y != warthog::INF);
                    this->labels_.back().at((*it).second).grow(x, y);
                    assert(this->labels_.back().at((*it).second).is_valid());
                };
        dijkstra.apply_to_closed(bbox_fn);
    }
    std::cerr << "\nall done\n"<< std::endl;
}

void
warthog::bb_filter::print(std::ostream& out)
{
    out << "# Bounding box edge labels for the graph file \n"
        << g_->get_filename() << "\n" 
        << "# Before printing labels are sorted by tail index. \n"
        << "# Ties are broken using the original input ordering \n"
        << "# (lower is better)\n";
    out << "firstnode " << start_id_ << " lastnode " << last_id_ << "\n";
    for(uint32_t i = start_id_; i <= last_id_; i++)
    {
        std::vector<warthog::geom::rectangle> vrect = labels_.at(i);
        for(uint32_t j = 0; j < vrect.size(); j++)
        {
            warthog::geom::rectangle rect = vrect.at(j);
            out << rect.x1 << "\t" << rect.y1 << "\t"
                << rect.x2 << "\t" << rect.y2 << "\t"
                << std::endl;
        }
    }
}

bool
warthog::bb_filter::load_labels(const char* filename)
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
    labels_.resize(g_->get_num_nodes());
    // skip comment lines
    while(ifs.peek() == '#')
    {
        while(ifs.good() && ifs.get() != '\n');
        lines++;
    }
    
    bool bad_header = false;
    std::string token;

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

    if(!(start_id_ < g_->get_num_nodes() && 
        last_id_ < g_->get_num_nodes() && 
        start_id_ < last_id_)) 
    {
        std::cerr << " invalid header values in file " << filename << "\n";
        return false;
    }

    for(uint32_t i = start_id_; i <= last_id_; i++)
    {
        warthog::graph::node* n = g_->get_node(i);
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
            labels_.at(i).push_back(warthog::geom::rectangle(x1, y1, x2, y2));
        }
    }
    ifs.close();
    return true;
}

void
warthog::bb_filter::init(warthog::problem_instance* instance)
{
    assert(instance);
    g_->get_xy(instance->get_goal(), tx, ty);
}
