#include "cpd.h"
#include "cpd/graph_oracle.h"
#include "helpers.h"

struct shared_data
{
    warthog::cpd::graph_oracle* cpd_;
    std::vector<uint32_t>* sources_;
};

// helps to precompute first-move data
struct graph_oracle_listener
{
    inline void
    generate_node(warthog::search_node* from, 
                  warthog::search_node* succ, 
                  warthog::cost_t edge_cost,
                  uint32_t edge_id) 
    {
        if(from == 0) { return; } // start node 

        if(from->get_id() == *source_id_) // start node successors
        { 
            //assert(s_row_.at(succ->get_id()) == 0);
            assert(edge_id < 
                     oracle_->get_graph()->get_node(
                        (uint32_t)*source_id_)->out_degree());
            s_row_->at(succ->get_id()) = (1 << edge_id);
            assert(s_row_->at(succ->get_id()));
        }
        else // all other nodes
        {
            warthog::sn_id_t succ_id = succ->get_id();
            warthog::sn_id_t from_id = from->get_id();
            double alt_g = from->get_g() + edge_cost;
            double g_val = 
                succ->get_search_number() == from->get_search_number() ? 
                succ->get_g() : DBL_MAX;

            //  update first move
            if(alt_g < g_val) 
            { 
                s_row_->at(succ_id) = s_row_->at(from_id);
                assert(s_row_->at(succ_id) == s_row_->at(from_id));
            }
            
            // add to the list of optimal first moves
            if(alt_g == g_val) 
            { 
                s_row_->at(succ_id) |= s_row_->at(from_id); 
                assert(s_row_->at(succ_id) >= s_row_->at(from_id));
            }

        }
    }

    inline void
    expand_node(warthog::search_node* current) { }

    inline void
    relax_node(warthog::search_node* current) { }

    warthog::cpd::graph_oracle* oracle_;
    warthog::sn_id_t* source_id_;
    std::vector<warthog::cpd::fm_coll>* s_row_;

};

void
warthog::cpd::graph_oracle::precompute() 
{

    // The actual precompute function. We fork threads and run this
    // function. Each thread considers a selected set of source nodes
    void*(*thread_compute_fn)(void*) = 
    [] (void* args_in) -> void*
    {
        warthog::helpers::thread_params* par = 
            (warthog::helpers::thread_params*) args_in;

        shared_data* shared = (shared_data*) par->shared_;
        warthog::cpd::graph_oracle* cpd = shared->cpd_;
        warthog::graph::xy_graph* g = cpd->get_graph();

        std::vector<warthog::cpd::fm_coll> s_row(g->get_num_nodes());
        warthog::sn_id_t source_id;

        // each thread has its own copy of Dijkstra and each
        // copy has a separate memory pool
        warthog::simple_graph_expansion_policy expander(g);
        warthog::zero_heuristic h;
        warthog::pqueue_min queue;
        graph_oracle_listener listener;

        warthog::flexible_astar 
            <warthog::zero_heuristic, 
            warthog::simple_graph_expansion_policy,
            warthog::pqueue_min,
            graph_oracle_listener>
                dijk(&h, &expander, &queue, &listener);

        listener.oracle_ = cpd;
        listener.source_id_ = &source_id;
        listener.s_row_ = &s_row;
        dijk.set_listener(&listener);

        for(uint32_t i = 0; i < shared->sources_->size(); i++)
        {
            // source nodes are evenly divided among all threads;
            // skip any source nodes not intended for current thread

            { continue; }

            source_id = shared->sources_->at(i);
            warthog::problem_instance problem(source_id);
            warthog::solution sol;

            // everything begins as a wildcard 
            s_row.clear();
            s_row.assign(shared->cpd_->get_graph()->get_num_nodes(), warthog::cpd::CPD_FM_NONE);

            // compute actual first moves to each reachable node
            // (NB: non-reachable nodes retain wildcard moves)
            dijk.get_path(problem, sol);
            cpd->add_row((uint32_t)source_id, s_row);
            par->nprocessed_++;
        }
        return 0;
    };

    // Here we do the forking. Each thread receives its own
    // separate version of Dijkstra's algorithm and each
    // has a separate memory pool

    warthog::timer t;
    t.start();

    warthog::cpd::compute_dfs_preorder(g_, &order_);

    // specify which source nodes to precompute for
    // (default: all)
    std::vector<uint32_t> source_nodes;
    source_nodes.reserve(g_->get_num_nodes());
    for(uint32_t i = 0; i < g_->get_num_nodes(); i++)
    {
        source_nodes.push_back(i);
    }

    shared_data shared;
    shared.cpd_ = this;
    shared.sources_ = &source_nodes;

    std::cerr << "computing dijkstra labels\n";
    warthog::helpers::parallel_compute(
            thread_compute_fn, &shared, 
            (uint32_t)source_nodes.size());
            
    // convert the column order into a map: from vertex id to its ordered index
    warthog::helpers::value_index_swap_array(order_);

    t.stop();
    std::cerr 
        << "total preproc time (seconds): "
        << t.elapsed_time_micro() / 1000000 << "\n";

}

void
warthog::cpd::graph_oracle::add_row(uint32_t source_id, 
                 std::vector<warthog::cpd::fm_coll>& row)
{
    // source gets a wildcard move
    row.at(source_id) = warthog::cpd::CPD_FM_NONE;

    // greedily compress the row w.r.t. the current column order
    warthog::cpd::fm_coll moveset = row.at(order_.at(0));
    uint32_t head = 0;
    for(uint32_t index = 0; index < row.size(); index++)
    {
        assert(moveset > 0);
        if((moveset & row.at(order_.at(index))) == 0)
        {
            uint32_t firstmove = __builtin_ffsl(moveset) - 1;
            assert(firstmove < warthog::cpd::CPD_FM_MAX);
            fm_.at(source_id).push_back(
                    warthog::cpd::rle_run32{ (head << 4) | firstmove} );
            moveset = row.at(order_.at(index));
            head = index;
        }
        moveset = moveset & row.at(order_.at(index));
    } 
    
    // add the last run
    uint32_t firstmove = __builtin_ffsl(moveset) - 1;
    assert(firstmove < warthog::cpd::CPD_FM_MAX);
    fm_.at(source_id).push_back(
            warthog::cpd::rle_run32{ (head << 4) | firstmove} );
    
//    std::cerr << "compressed source row " << source_id << " with " 
//        << fm_.at(source_id).size() << std::endl;
}

std::istream&
warthog::cpd::operator>>(std::istream& in, 
        warthog::cpd::graph_oracle& lab)
{
    // read the graph size data
    warthog::timer mytimer;
    mytimer.start();

    uint32_t num_nodes;
    in.read((char*)(&num_nodes), 4);
    if(num_nodes != lab.g_->get_num_nodes())
    {
        std::cerr 
            << "err; " << "input mismatch. cpd file says " << num_nodes 
            << " nodes, but graph contains " << lab.g_->get_num_nodes() << "\n";
        return in;
    }

    lab.fm_.clear();
    lab.order_.resize(lab.g_->get_num_nodes());

    // read the vertex-to-column-order mapping
    for(uint32_t i = 0; i < num_nodes; i++)
    {
        uint32_t n_id; 
        in.read((char*)(&n_id), 4);
        lab.order_.at(i) = n_id;
    }

    // read the RLE data
    uint32_t run_count = 0;
    lab.fm_.resize(lab.g_->get_num_nodes());
    for(uint32_t row_id = 0; row_id < lab.g_->get_num_nodes(); row_id++)
    {
        // number of runs for this row
        uint32_t num_runs;
        in.read((char*)(&num_runs), 4);

        // read all the runs for the current row
        for(uint32_t i = 0; i < num_runs; i++)
        {
            warthog::cpd::rle_run32 tmp;
            in >> tmp;
            lab.fm_.at(row_id).push_back(tmp);
            run_count++;

            if(!in.good())
            {
                std::cerr << "err; while reading firstmove labels\n";
                std::cerr 
                    << "[debug info] "
                    << " row_id " << row_id 
                    << " run# " << i << " of " << lab.fm_.size()
                    << ". aborting.\n";
                return in;
            }
        }
    }
    mytimer.stop();

    std::cerr 
        << "read from disk " << lab.fm_.size() 
        << " rows and "
        << run_count << " runs. "
        << " time: " << (double)mytimer.elapsed_time_nano() / 1e9 << " s\n";
    return in;
}

std::ostream&
warthog::cpd::operator<<(std::ostream& out,
        warthog::cpd::graph_oracle& lab)
{
    warthog::timer mytimer;
    mytimer.start();

    // write graph size 
    uint32_t num_nodes = lab.g_->get_num_nodes();
    out.write((char*)(&num_nodes), 4);

    // write node ordering
    assert(lab.order_.size() == num_nodes);
    for(uint32_t i = 0; i < num_nodes; i++)
    {
        uint32_t n_id = lab.order_.at(i);
        out.write((char*)(&n_id), 4);
    }

    // write the runs for each row
    uint32_t run_count = 0;
    for(uint32_t row_id = 0; row_id < lab.g_->get_num_nodes(); row_id++)
    {
        // write the number of runs
        uint32_t num_runs = (uint32_t)lab.fm_.at(row_id).size();
        out.write((char*)(&num_runs), 4);

        for(uint32_t run = 0; run < num_runs; run++)
        {
            out << lab.fm_.at(row_id).at(run);
            run_count++;
            if(!out.good())
            {
                std::cerr << "err; while writing labels\n";
                std::cerr 
                    << "[debug info] "
                    << " row_id " << row_id 
                    << " run# " << lab.fm_.at(row_id).size()
                    << ". aborting.\n";
                return out;
            }
        }
        
    }
    mytimer.stop();

    std::cerr 
        << "wrote to disk " << lab.fm_.size()
        << " rows and "
        << run_count << " runs. "
        << " time: " << (double)mytimer.elapsed_time_nano() / 1e9 << " s \n";
    return out;
}

