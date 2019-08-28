#include "ch_data.h"

warthog::ch::ch_data::ch_data()
{
    g_ = new warthog::graph::xy_graph();
    level_ = new std::vector<uint32_t>(0);
}

warthog::ch::ch_data::~ch_data()
{
    delete g_;
    delete level_;
}

bool
warthog::ch::save_ch_data(const char* filename, warthog::ch::ch_data* chd)
{
    std::ofstream fs_out(filename);  
    
    // comments
    fs_out << "# contraction hierarchy data generated from\n"
        << "# " << chd->g_->get_filename() << std::endl
        << "# we use 32bit integer data formatted as follows: node data first, then edge data\n"
        << "# node data format: v [id] [x] [y] [level]\n"
        << "# edge data format: e [from_node_id] [to_node_id] [cost]" << std::endl;

    // header stuff
    //fs_out << "chd 1.0" << std::endl;
    fs_out 
        << "nodes " << chd->g_->get_num_nodes() << " "
        << "edges " << chd->g_->get_num_edges_out() << std::endl;

    // node data
    for(uint32_t i = 0; i < chd->g_->get_num_nodes(); i++)
    {
        int32_t x, y; 
        chd->g_->get_xy(i, x, y);
        fs_out 
            << "v " << i << " " 
            << x << " " 
            << y << " ";

         uint32_t level;
         if(i >= chd->level_->size()) 
         { level = (uint32_t) chd->level_->size()-1; }
         else 
         { level = chd->level_->at(i); }
         fs_out << level << std::endl;
    }

    // edge data
    for(uint32_t i = 0; i < chd->g_->get_num_nodes(); i++)
    {
        warthog::graph::node* n = chd->g_->get_node(i);
        for(uint32_t edge_idx = 0; edge_idx < n->out_degree(); edge_idx++)
        {
            warthog::graph::edge* e = n->outgoing_begin()+edge_idx;
            fs_out << "e " << i << " " << e->node_id_ << " " << e->wt_ << std::endl;
        }
    }

    fs_out.close();
    return true;
}

warthog::ch::ch_data* 
warthog::ch::load_ch_data(const char* filename, bool up_only)
{
    warthog::ch::ch_data* chd = 0;
    std::ifstream fs_in(filename);

    uint32_t num_nodes=0, num_edges=0;
    while(fs_in.good())
    {
        fs_in >> std::ws;
        if(fs_in.peek() == '#') 
        { 
            while(fs_in.get() != '\n'); 
            continue;
        }

        if(fs_in.peek() == 'n') { while(fs_in.get() != ' '); } // "nodes" keyword
        fs_in >> num_nodes;
        fs_in >> std::ws;
        if(fs_in.peek() == 'e') { while(fs_in.get() != ' '); } // "edges" keyword
        fs_in >> num_edges;
        fs_in >> std::ws;
        break;
    }

    chd = new warthog::ch::ch_data();
    chd->g_ = new warthog::graph::xy_graph(num_nodes);
    chd->level_ = new std::vector<uint32_t>(num_nodes, num_nodes-1);

    uint32_t n_added = 0, e_added=0;
    while(fs_in.good())
    {
        // read nodes data
        fs_in >> std::ws;
        while(fs_in.peek() == 'v')
        {
            uint32_t id, level;
            int32_t x, y;
            fs_in.get(); // eat the 'v' char
            fs_in >> id >> x >> y >> level;
            chd->g_->set_xy(id, x, y);
            chd->level_->at(id) = level;
            fs_in >> std::ws; // trailing whitespace
            n_added++;
        }

        while(fs_in.peek() == 'e')
        {
            uint32_t from_id, to_id;
            warthog::graph::edge_cost_t cost;

            fs_in.get(); // eat the 'e' char
            fs_in >> from_id >> to_id >> cost;
            if(chd->level_->at(from_id) <= chd->level_->at(to_id))
            {
                warthog::graph::node* from = chd->g_->get_node(from_id);
                from->add_outgoing(warthog::graph::edge(to_id, cost));
                if(up_only == false)
                {
                    warthog::graph::node* to = chd->g_->get_node(to_id);
                    to->add_incoming(warthog::graph::edge(from_id, cost));
                }
            }
            else
            {
                if(up_only == false)
                {
                    warthog::graph::node* from = chd->g_->get_node(from_id);
                    from->add_outgoing(warthog::graph::edge(to_id, cost));
                }
                warthog::graph::node* to = chd->g_->get_node(to_id);
                to->add_incoming(warthog::graph::edge(from_id, cost));
            }
            e_added++;
            fs_in >> std::ws;
        }
    }
    fs_in.close();
    std::cerr << "graph, loaded.\n";
    std::cerr 
        << "read " << n_added << " nodes (total " << num_nodes << ")"
        << " and read " << e_added << " edges (total "<< num_edges << ")"
        << std::endl;
    return chd;
}
