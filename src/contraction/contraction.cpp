#include "contraction.h"
#include "planar_graph.h"
#include <algorithm>
#include <stack>
#include <vector>

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
        ofs << "# node ids, in order of contraction\n";
        for(uint32_t i = 0; i < order.size(); i++)
        {
            ofs << order.at(i) << "\n";
        }
    }
    else
    {
        std::cerr << "err; cannot write to file " << filename << std::endl;
    }
    ofs.close();
}

bool
warthog::ch::load_node_order(const char* filename, std::vector<uint32_t>& order,
        bool lex_order)
{
    order.clear();
    std::ifstream ifs(filename, std::ios_base::in);
    if(!ifs.good())
    {
        std::cerr << "\nerror trying to load node order from file " 
            << filename << std::endl;
        ifs.close();
        return false;
    }

    while(true)
    {
        // skip comment lines
        while(ifs.peek() == '#')
        {
            while(ifs.get() != '\n');
        }

        uint32_t tmp;
        ifs >> tmp;
        if(!ifs.good()) { break; }
        order.push_back(tmp);
    }
    ifs.close();

    if(lex_order)
    {
        warthog::ch::value_index_swap_dimacs(order);
    }
    return true;
}

typedef std::set<uint32_t>::iterator set_iter;
void
warthog::ch::compute_closure(uint32_t source, 
        warthog::graph::planar_graph* g, std::set<uint32_t>* closure, 
        uint32_t maxdepth)
{
    std::stack<std::pair<uint32_t, uint32_t>> stack; // stack of node ids
    stack.push(std::pair<uint32_t, uint32_t>(source, 0));
    while(stack.size() != 0)
    {
        // add top element to the closure
        std::pair<uint32_t, uint32_t> top = stack.top();
        uint32_t next_id = top.first;
        uint32_t depth = top.second;
        std::pair<std::set<uint32_t>::iterator, bool> ret 
            = closure->insert(next_id);
        stack.pop();
        if(ret.second == false) 
        { 
            continue;  // already in the closure
        }

        if(depth >= maxdepth)
        { 
            continue; // limited depth closure
        }

        // add all outgoing neighbours to the stack
        warthog::graph::node* n = g->get_node(next_id);
        for( warthog::graph::edge_iter it = n->outgoing_begin(); 
                it != n->outgoing_end(); 
                it++)
        {
            stack.push(std::pair<uint32_t, uint32_t>((*it).node_id_, depth+1));
        }
    }
}

void
warthog::ch::value_index_swap_dimacs(
        std::vector<uint32_t>& vec)
{
    // re-maps @param vec s.t. for each x and i
    // v[i] = x becomes v[x] = i
    std::vector<uint32_t> tmp;
    tmp.resize(vec.size());
    for(uint32_t i = 0; i < vec.size(); i++)
    {
        tmp.at(i) = vec.at(i);
    }

    // because DIMACS graphs are stupid and 1-indexed 
    // we need to insert extra elements for padding
    uint32_t min_id = *(std::min_element(tmp.begin(), tmp.end()));
    vec.clear();
    vec.resize(tmp.size() + min_id);
    for(uint32_t i = 0; i < min_id; i++)
    {
        vec.at(i) = 0;
    }

    for(uint32_t i = 0; i < tmp.size(); i++)
    {
        vec.at(tmp.at(i)) = i + min_id;
    }
}

void
warthog::ch::partition_greedy_bottom_up(
        warthog::graph::planar_graph* g, 
        std::vector<uint32_t>* rank,
        uint32_t nparts,
        std::vector<uint32_t>* part)
{

}

void
warthog::ch::unpack(uint32_t from_id,
        warthog::graph::edge_iter it_e,
        warthog::graph::planar_graph* g,
        std::set<uint32_t>& intermediate)
{
    warthog::graph::node* from = g->get_node(from_id);
    assert(it_e >= from->outgoing_begin() && it_e < from->outgoing_end());

    warthog::graph::edge* e_ft = &*it_e;
    uint32_t to_id = e_ft->node_id_;

    for(warthog::graph::edge_iter it = from->outgoing_begin(); 
            it < from->outgoing_end(); it++)
    {
        warthog::graph::node* succ = g->get_node((*it).node_id_);
        warthog::graph::edge_iter it_e_succ = succ->outgoing_begin();
        while(true)
        {
            it_e_succ = succ->find_edge(to_id, it_e_succ);
            if(it_e_succ == succ->outgoing_end()) { break; }
            assert( it_e_succ >= succ->outgoing_begin() && 
                    it_e_succ <= succ->outgoing_end());

            warthog::graph::edge* e_st = &*it_e_succ;
            if(((*it).wt_ + e_st->wt_) == e_ft->wt_) { break; }
            it_e_succ++;
        }
        if(it_e_succ == succ->outgoing_end()) { continue; } 

        intermediate.insert(from_id);
        intermediate.insert((*it).node_id_);

        // recursively unpack the two edges being represented by
        // the single shortcut edge (from_id, to_id)
        unpack(from_id, it, g, intermediate);

        succ = g->get_node((*it).node_id_);
        unpack((*it).node_id_, it_e_succ, g, intermediate);
        break;
    }
}

void
warthog::ch::unpack_list_edges(uint32_t from_id,
        warthog::graph::edge_iter it_e,
        warthog::graph::planar_graph* g,
        std::vector<warthog::graph::edge*>& edges)
{
    warthog::graph::node* from = g->get_node(from_id);
    assert(it_e >= from->outgoing_begin() && it_e < from->outgoing_end());

    warthog::graph::edge* e_ft = &*it_e;
    uint32_t to_id = e_ft->node_id_;

    bool unpacked = false;
    for(warthog::graph::edge_iter it = from->outgoing_begin(); 
            it < from->outgoing_end(); it++)
    {
        warthog::graph::node* succ = g->get_node((*it).node_id_);
        warthog::graph::edge_iter it_e_succ = succ->outgoing_begin();
        while(true)
        {
            it_e_succ = succ->find_edge(to_id, it_e_succ);
            if(it_e_succ == succ->outgoing_end()) { break; }
            assert( it_e_succ >= succ->outgoing_begin() && 
                    it_e_succ <= succ->outgoing_end());

            warthog::graph::edge* e_st = &*it_e_succ;
            if(((*it).wt_ + e_st->wt_) == e_ft->wt_) { break; }
            it_e_succ++;
        }
        if(it_e_succ == succ->outgoing_end()) { continue; } 

        // recursively unpack the two edges being represented by
        // the single shortcut edge (from_id, to_id)
        unpack_list_edges(from_id, it, g, edges);

        succ = g->get_node((*it).node_id_);
        unpack_list_edges((*it).node_id_, it_e_succ, g, edges);
        unpacked = true;
        break;
    }

    if(!unpacked)
    {
        edges.push_back(e_ft);
    }
}
