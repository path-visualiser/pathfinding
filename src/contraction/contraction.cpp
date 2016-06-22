#include "contraction.h"
#include "planar_graph.h"
#include <stack>

void
warthog::ch::make_input_order(warthog::graph::planar_graph& g, std::vector<uint32_t>& order)
{
    order.clear();
    order.reserve(g.get_num_nodes());
    for(uint32_t i=0; i < g.get_num_nodes(); i++)
    {
        order.push_back(i);
    }
}

void 
warthog::ch::write_node_order(const char* filename, std::vector<uint32_t>& order)
{
    std::ofstream ofs(filename, std::ios_base::out | std::ios_base::trunc);
    if(ofs.good())
    {
        for(uint32_t i = 0; i < order.size(); i++)
        {
            ofs << order.at(i) << "\n";
        }
    }
    ofs.close();
}

void 
warthog::ch::load_node_order(const char* filename, std::vector<uint32_t>& order)
{
    order.clear();
    std::ifstream ifs(filename, std::ios_base::in);
    if(!ifs.good())
    {
        std::cerr << "\nerror trying to load node order from file " 
            << filename << std::endl;
        ifs.close();
        return;
    }

    while(ifs.good())
    {
        uint32_t tmp;
        ifs >> tmp;
        order.push_back(tmp);
    }
    ifs.close();
}

typedef std::set<uint32_t>::iterator set_iter;
void
warthog::ch::compute_closure(uint32_t source, 
        warthog::graph::planar_graph* g, std::set<uint32_t>* closure)
{
    std::stack<uint32_t> stack; // stack of node ids
    stack.push(source);
    while(stack.size() != 0)
    {
        uint32_t next_id = stack.top();
        std::pair<std::set<uint32_t>::iterator, bool> ret 
            = closure->insert(next_id);
        stack.pop();
        if(ret.second == false) 
        { 
            continue;  // already in the closure
        }

        // add all outgoing neighbours to the stack
        warthog::graph::node* n = g->get_node(next_id);
        for( warthog::graph::edge_iter it = n->outgoing_begin(); 
                it != n->outgoing_end(); 
                it++)
        {
            stack.push((*it).node_id_);
        }
    }
}

