// cpd/cpd.h
//
// Various defintions / data-structures for CPDs
//
// @author: dharabor
// @created: 2020-02-25
//

#ifndef WARTHOG_CPD_CPD_H
#define WARTHOG_CPD_CPD_H

#include <cassert>
#include <functional>
#include <iostream>
#include <vector>

#include "flexible_astar.h"
#include "helpers.h"
#include "search_node.h"
#include "timer.h"
#include "zero_heuristic.h"

namespace warthog
{
namespace cpd
{

// defines a 32bit run in a run-length encoding
// the first 4 bits tell the symbol.
// the next 28 bits tell the index where the run begins
struct rle_run32
{
    uint8_t 
    get_move() { return data_ & 0xF; } 

    uint32_t 
    get_index() { return data_ >> 4; } 


    void
    print(std::ostream& out)
    {
        out << " [" << get_index() << ", " << get_move() << "]";
    }

    uint32_t data_;
};

//  limits on the number of nodes in a graph 
//  for which we compute a CPD
static const uint32_t RLE_RUN32_MAX_INDEX = (UINT32_MAX >> 4);

// limits on the maximum number of first-move labels that need to be stored.
// this value should be greater than the maximum degree of any node plus one
// extra value for the case where a node is unreachable
static const uint32_t CPD_FM_MAX = 16;

// special value to denote that no first move exists.
static const uint32_t CPD_FM_NONE=0xF;

// a collection of optimal first moves. 
// we keep one bit for each optimal move. the maximum
// degree of any node is determined by CPD_FM_MAX
typedef uint16_t fm_coll;

// a DFS pre-order traversal of the input graph is known to produce an 
// effective node order for the purposes of compression
// @param g: the input graph
// @param column_order: a list of node ids as visited by DFS from a random
// start node
void
compute_dfs_preorder(
        warthog::graph::xy_graph* g,
        std::vector<uint32_t>* column_order);

// compute and compress labels for all nodes specified by 
// @param workload_manager 
template <typename T_EXPANDER, typename T_CPD>
void
compute_firstmoves_and_compress(
        T_CPD* cpd,
        T_EXPANDER* expander,
        std::vector<uint32_t>* source_nodes)
{
    warthog::timer t;
    t.start();

    if(cpd == 0 || source_nodes == 0 || expander == 0) { return; } 

    
    struct shared_data
    {
        T_CPD* cpd_;
        T_EXPANDER* expander_;
        std::vector<uint32_t>* sources_;
    };

    // The actual precompute function. We run Dijkstra 
    // for a selected sets of source nodes
    void*(*thread_compute_fn)(void*) = 
    [] (void* args_in) -> void*
    {
        warthog::helpers::thread_params* par = 
            (warthog::helpers::thread_params*) args_in;
        shared_data* shared = (shared_data*) par->shared_;

        // bookkeeping data: track the current row and firstmove labels 
        warthog::sn_id_t source_id;
        std::vector<warthog::cpd::fm_coll>
            s_row(shared->cpd_->get_graph()->get_num_nodes());


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
                assert(s_row.at(succ->get_id()) == 0);
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

        warthog::zero_heuristic h;
        warthog::pqueue_min queue;
        warthog::flexible_astar <warthog::zero_heuristic, T_EXPANDER> 
                dijk(&h, shared->expander_, &queue);
        dijk.apply_on_generate(on_generate_fn);

        for(uint32_t i = 0; i < shared->sources_->size(); i++)
        {
            // source nodes are evenly divided among all threads;
            // skip any source nodes not intended for current thread
//            if((i % par->max_threads_) != par->thread_id_) 
//            { continue; }

            source_id = shared->sources_->at(i);
            warthog::problem_instance problem(source_id);
            warthog::solution sol;

            s_row.clear();
            s_row.resize(shared->cpd_->get_graph()->get_num_nodes());
            dijk.get_path(problem, sol);
            shared->cpd_->add_row((uint32_t)source_id, s_row);
            par->nprocessed_++;
        }
        return 0;
    };

    shared_data shared;
    shared.cpd_ = cpd;
    shared.expander_ = expander;
    shared.sources_ = source_nodes;
    warthog::helpers::thread_params params;
    params.shared_ = &shared;


    thread_compute_fn(&params);

//
//    std::cerr << "computing dijkstra labels\n";
//    warthog::helpers::parallel_compute(
//            thread_compute_fn, &shared, 
//            (uint32_t)source_nodes->size());
//            

    t.stop();
    std::cerr 
        << "total preproc time (seconds): "
        << t.elapsed_time_micro() / 1000000 << "\n";
}


}

}

#endif

