#include "cpd/graph_oracle.h"

struct shared_data
{
    warthog::cpd::graph_oracle* cpd_;
    warthog::simple_graph_expansion_policy* expander_;
    std::vector<uint32_t>* sources_;
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

        // callback function used to record the optimal first move 
        std::function<void(
                warthog::search_node*, 
                warthog::search_node*, double, uint32_t)>  
            on_generate_fn = [&source_id, &s_row, shared]
        (warthog::search_node* succ, warthog::search_node* from,
             double edge_cost, uint32_t edge_id) -> void
        {
            if(from == 0) { return; } // start node 

            if(from->get_id() == source_id) // start node successors
            { 
                //assert(s_row.at(succ->get_id()) == 0);
                assert(edge_id < 
                         shared->cpd_->get_graph()->get_node(
                            (uint32_t)source_id)->out_degree());
                s_row.at(succ->get_id()) = (1 << edge_id);
                assert(s_row.at(succ->get_id()));
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
                    s_row.at(succ_id) = s_row.at(from_id);
                    assert(s_row.at(succ_id) == s_row.at(from_id));
                }
                
                // add to the list of optimal first moves
                if(alt_g == g_val) 
                { 
                    s_row.at(succ_id) |= s_row.at(from_id); 
                    assert(s_row.at(succ_id) >= s_row.at(from_id));
                }
            }
        };

        // each thread has its own copy of Dijkstra and each
        // copy has a separate memory pool
        warthog::simple_graph_expansion_policy expander(g);
        warthog::zero_heuristic h;
        warthog::pqueue_min queue;
        warthog::flexible_astar 
            <warthog::zero_heuristic, warthog::simple_graph_expansion_policy> 
                dijk(&h, &expander, &queue);
        dijk.apply_on_generate(on_generate_fn);

        for(uint32_t i = 0; i < shared->sources_->size(); i++)
        {
            // source nodes are evenly divided among all threads;
            // skip any source nodes not intended for current thread
            if((i % par->max_threads_) != par->thread_id_) 
            { continue; }

            source_id = shared->sources_->at(i);
            warthog::problem_instance problem(source_id);
            warthog::solution sol;

            s_row.clear();
            s_row.resize(shared->cpd_->get_graph()->get_num_nodes());
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
    warthog::cpd::fm_coll current = row.at(order_.at(0));
    uint32_t head = 0;
    for(uint32_t index = 0; index < row.size(); index++)
    {
        assert(current > 0);
        uint32_t next_id = order_.at(index);
        if(!(current & row.at(next_id)))
        {
            uint32_t firstmove = __builtin_ffsl(current) - 1;
            assert(firstmove < warthog::cpd::CPD_FM_MAX);
            fm_.at(source_id).push_back(
                    warthog::cpd::rle_run32{ (head << 4) | firstmove} );
            current = row.at(next_id);
            head = index;
        }
        current = current & row.at(next_id);
    } 
    
    // add the last run
    uint32_t firstmove = __builtin_ffsl(current) - 1;
    assert(firstmove < warthog::cpd::CPD_FM_MAX);
    fm_.at(source_id).push_back(
            warthog::cpd::rle_run32{ (head << 4) | firstmove} );
    
//    std::cerr << "compressed source row " << source_id << " with " 
//        << fm_.at(source_id).size() << std::endl;
}

