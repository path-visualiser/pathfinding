#include "firstmove_labelling.h"
#include "planar_graph.h"

#include <algorithm>

warthog::label::firstmove_labelling::firstmove_labelling(
        warthog::graph::planar_graph* g)
    : g_(g)
{
    lab_ = new std::vector< std::vector < warthog::label::fm_run > > ();
    lab_->resize(g_->get_num_nodes());
}

warthog::label::firstmove_labelling::~firstmove_labelling()
{ 
    delete lab_;
}

std::istream&
warthog::label::operator>>(std::istream& in, fm_run& the_run)
{
    in.read((char*)(&the_run.head_), 4);
    in.read((char*)(&the_run.label_), 1);
    return in ;
}

std::ostream&
warthog::label::operator<<(std::ostream& out, fm_run& the_run)
{
    out.write((char*)(&the_run.head_), 4);
    out.write((char*)(&the_run.label_), 1);
    return out;
}

std::istream&
warthog::label::operator>>(std::istream& in, 
        warthog::label::firstmove_labelling& lab)
{
    lab.lab_->resize(lab.g_->get_num_nodes());

    uint32_t num_rows;
    in.read((char*)(&num_rows), 4);

    for(uint32_t row = 0; row < num_rows; row++)
    {
        uint32_t row_id; 
        in.read((char*)(&row_id), 4);

        if(row_id > lab.g_->get_num_nodes())
        {
            std::cerr << "err; " << row_id 
                << " is aninvalid row identifier. aborting.\n";
            break;
        }

        if(lab.lab_->at(row_id).size() != 0)
        {
            std::cerr << "err; row id "<< row_id 
                << " appears more than once. aborting.\n";
            break;
        }

        uint32_t num_runs;
        in.read((char*)(&num_runs), 4);

        // read all the runs for the current row
        for(uint32_t i = 0; i < num_runs; i++)
        {
            fm_run tmp;
            in >> tmp;
            lab.lab_->at(row_id).push_back(tmp);

            assert(in.good());

            if(!in.good())
            {
                std::cerr << "err; while reading firstmove labels\n";
                std::cerr 
                    << "[debug info] row# " << row
                    << " row_id " << row_id 
                    << " run# " << i << " of " << num_runs 
                    << ". aborting.\n";
                return in;
            }
        }
    }
    return in;
}

std::ostream&
warthog::label::operator<<(std::ostream& out,
        warthog::label::firstmove_labelling& lab)
{
    uint32_t num_rows = 0;
    for(uint32_t n_id = 0; n_id < lab.g_->get_num_nodes(); n_id++)
    {
        if(lab.lab_->size() > 0) { num_rows++; }
    }
    out.write((char*)(&num_rows), 4);

    for(uint32_t row_id = 0; row_id < lab.g_->get_num_nodes(); row_id++)
    {
        if(lab.lab_->size() == 0) { continue; }
        out.write((char*)(&row_id), 4);

        uint32_t num_runs = lab.lab_->at(row_id).size();
        out.write((char*)(&num_runs), 4);

        for(uint32_t run = 0; run < num_runs; run++)
        {
            out << lab.lab_->at(row_id).at(run);
            if(!out.good())
            {
                std::cerr << "err; while writing labels\n";
                std::cerr 
                    << "[debug info] "
                    << " row_id " << row_id 
                    << " run# " << run 
                    << ". aborting.\n";
                return out;
            }
        }
    }
    return out;
}

void
warthog::label::compute_fm_dfs_preorder(
        warthog::graph::planar_graph& g, 
        std::vector<uint32_t>& column_order)
{
    typedef std::pair<uint32_t, uint32_t> dfs_pair;

    // dfs bookkeeping (we store some data to avoid cycles)
    std::vector<uint32_t> tmp(g.get_num_nodes(), INT32_MAX); 
    std::function<void(uint32_t)> dfs_preorder =
        [&column_order, &tmp, &g](uint32_t dfs_seed) -> void
    {
        std::vector<dfs_pair> dfs_stack;
        dfs_stack.push_back(dfs_pair(dfs_seed, 0));
        dfs_stack.reserve(g.get_num_nodes()*2);
        while(dfs_stack.size())
        {
            // current node in the dfs tree 
            dfs_pair dfs_node = dfs_stack.back();
            dfs_stack.pop_back();
            warthog::graph::node* n = g.get_node(dfs_node.first);

            //std::cerr << dfs_node.first << " " << dfs_node.second;

            // pre-order traversal
            assert(dfs_node.first != UINT32_MAX);
            if(tmp.at(dfs_node.first) == INT32_MAX)
            { 
                tmp.at(dfs_node.first) = dfs_node.first;
                column_order.push_back(dfs_node.first); 
            }

            warthog::graph::edge_iter begin = n->outgoing_begin();
            warthog::graph::edge_iter left = begin + dfs_node.second;

            // right branch 
            for(uint32_t i = dfs_node.second+1; i < n->out_degree(); i++)
            {
                assert((begin+i)->node_id_ != UINT32_MAX);
                if(tmp.at((begin+i)->node_id_) == INT32_MAX)
                {
                    //std::cerr << " right " << (begin+i)->node_id_;
                    dfs_stack.push_back(dfs_pair(dfs_node.first, i));
                    break;
                }
            }

            // left branch 
            if(left == n->outgoing_end()) { return; }
            assert(left->node_id_ != UINT32_MAX);
            if(tmp.at(left->node_id_) == INT32_MAX)
            {
                //std::cerr << " left " << left->node_id_;
                dfs_stack.push_back(dfs_pair(left->node_id_, 0)); 
            }

            //std::cerr << "\n";
        }
    };

    // pick a random seed
    while(true)
    {
        std::vector<uint32_t> seeds;
        for(uint32_t i = 0; i < tmp.size(); i++)
        {
            if(tmp.at(i) == INT32_MAX) { seeds.push_back(i); } 
        }
        if(!seeds.size()) { break; }

        uint32_t dfs_seed = seeds.at(rand() % seeds.size());
        dfs_preorder(dfs_seed);
    }

    std::function<void(uint32_t, uint32_t)> dfs_order = 
    [&dfs_order, &tmp, &column_order, &g](uint32_t source_id, uint32_t depth) -> void
    {
        std::cerr << source_id << " " << depth << std::endl;
        warthog::graph::node* source = g.get_node(source_id);
        warthog::graph::edge_iter begin = source->outgoing_begin();
        for( warthog::graph::edge_iter it = begin; 
                it != source->outgoing_end();
                it++)
        {
            // recurse
            if(tmp.at(it->node_id_) == INT32_MAX)
            { 
                tmp.at(it->node_id_) = it->node_id_;
                dfs_order(it->node_id_, depth+1); 
            }
        }
        if(tmp.at(source_id) == INT32_MAX)
        { tmp.at(source_id) = source_id; }
        column_order.push_back(source_id);
    };
}
