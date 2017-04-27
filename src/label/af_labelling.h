#ifndef WARTHOG_LABEL_AF_LABELLING_H
#define WARTHOG_LABEL_AF_LABELLING_H

//  label/af_labelling.h
//
// Partition the nodes of the graph into disjoint sets and store with 
// every outgoing edge of every node a bitfield where the number of bits
// is equal to the number of disjoint sets in the partition. The individual
// bits are set to true if the edge at hand appears on any optimal path
// from the current node to any node in its associated partition. 
// The bit is set false if the edge appears on no such path.
//
// More details in the following paper:
//
// [Fast Point-to-Point Shortest Path Computations with Arc-Flags,
//  Ekkehard K{\"o}hler, Rolf H. M\{"o}hring, and Heiko Schilling,
//  In The Shortest Path Problem: Ninth DIMACS Implementation Challenge, 
//  Edited by Demetrescu, Camil and Goldberg, Andrew V. and Johnson, David S 
//  pp 41-72, American Mathematical Society, 2009]
//
//  @author: dharabor
//  @created: 2017-04-15
//

#include "flexible_astar.h"
#include "planar_graph.h"
#include "search_node.h"
#include "zero_heuristic.h"

#include <algorithm>
#include <functional>
#include <iostream>
#include <pthread.h>
#include <unistd.h>

namespace warthog
{

namespace label
{

class af_labelling
{
    public:

        virtual ~af_labelling();

        // print/serialise an arcflag labeling
        // 
        // @param out: (out) the target stream to write to
        void
        print(std::ostream& out);

        inline uint8_t*
        get_label(uint32_t node_id, uint32_t edge_index)
        { 
            assert(node_id < flags_->size() && 
                   edge_index < flags_->at(node_id).size());
            return flags_->at(node_id).at(edge_index);
        }

        inline size_t
        get_label_size() { return bytes_per_label_; }

        inline std::vector<uint32_t>*
        get_partitioning() { return part_; }
        
        // load up an arcflag labeling
        //
        // @param filename: the file to load from
        // @param g: the graph associated with the labelling
        // @param partitioning: a disjoint partition of the nodes in @param g
        // @return: a pointer to a newly created af_labelling for @param g
        static warthog::label::af_labelling*
        load(const char* filename, warthog::graph::planar_graph* g,
                std::vector<uint32_t>* partitioning);

        // create a new arcflag edge labelling 
        //
        // @param g: the input graph
        // @param part: a disjoint partitioning of the nodes in @param g
        // @param expander: defines the type of search that will be used
        // @return: a pointer to a newly created af_labelling for @param g
        template <typename t_expander>
        static warthog::label::af_labelling*
        compute(warthog::graph::planar_graph* g, std::vector<uint32_t>* part, 
                std::function<t_expander*(void)>& fn_new_expander)
        {
            if(g == 0 || part  == 0) { return 0; } 
            std::cerr << "computing arcflag labels\n";
            
            warthog::label::af_labelling* afl = 
                new warthog::label::af_labelling(g, part);

            // pthreads setup
            struct afl_params
            {
                uint32_t thread_id_;
                uint32_t max_threads_;
                uint32_t nprocessed_;
                t_expander* expander_;
                warthog::label::af_labelling* afl_;
            };

            // the actual precompute function
            void*(*thread_compute_fn)(void*) = [] (void* args_in) -> void*
            {
                afl_params* par = (afl_params*) args_in;
                std::cerr << "thread " << par->thread_id_ << std::endl;
            
                // helper function to track the first edge on the path from 
                // the source node to each node in the graph. we apply this 
                // function  following every successful node relaxation.
                // (the solution is a bit hacky as we break the chain of 
                // backpointers to achieve this; it doesn't matter, we don't 
                // care about the path)
                std::function<void(warthog::search_node*)> relax_fn = 
                    [] (warthog::search_node* n) -> void
                    {
                        // the start node and its children don't need 
                        // their parent pointers updated. for all other 
                        // nodes the grandparent becomes the parent
                        if(n->get_parent()->get_parent() != 0)
                        {
                            if(n->get_parent()->get_parent()->get_parent() 
                                    != 0)
                            {
                                n->set_parent(n->get_parent()->get_parent());
                            }
                        }
                    };

                warthog::zero_heuristic heuristic;
                warthog::flexible_astar< warthog::zero_heuristic, t_expander >
                    dijkstra(&heuristic, par->expander_);
                dijkstra.apply_on_relax(relax_fn);

                for(uint32_t i = 0; i < par->afl_->g_->get_num_nodes(); i++)
                {
                    // source nodes are evenly divided among all threads;
                    // skip any source nodes not intended for current thread
                    if((i % par->max_threads_) != par->thread_id_) 
                    { continue; }

                    // process the source
                    uint32_t source_id = i;
                    warthog::graph::node* source = 
                        par->afl_->g_->get_node(source_id);

                    // run a dijkstra search from the current source node:
                    // and then analyse the closed list to compute arc flags
                    // (there's better ways to do this but, bleh)
                    uint32_t ext_source_id = 
                        par->afl_->g_->to_external_id(source_id);
                    warthog::problem_instance pi(ext_source_id, warthog::INF);
                    dijkstra.get_length(pi);

                    // first, we need an easy way to convert between the ids 
                    // of nodes adjacent to the source and their corresponding 
                    // edge index
                    std::unordered_map<uint32_t, uint32_t> idmap;
                    uint32_t edge_idx = 0;
                    for(warthog::graph::edge_iter it = 
                            source->outgoing_begin(); 
                            it != source->outgoing_end(); 
                            it++)
                    {
                        idmap.insert(std::pair<uint32_t, uint32_t>(
                                    (*it).node_id_, edge_idx));
                        edge_idx++;
                    }

                    // analyse the nodes on the closed list and label the 
                    // edges of the source node accordingly
                    std::function<void(warthog::search_node*)> fn_arcflags =
                        [par, &source_id, &idmap](warthog::search_node* n)
                        -> void
                        {
                            // skip the source
                            assert(n);
                            if(n->get_id() == source_id) { return; } 
                            assert(n->get_parent());

                            // label the edges of the source
                            // (TODO: make this stuff faster)
                            uint32_t part_id = 
                                par->afl_->part_->at(n->get_id());
                            uint32_t e_idx  = 
                                (*idmap.find(
                                    n->get_parent()->get_id() == source_id ?  
                                    n->get_id() : 
                                    n->get_parent()->get_id())).second;
                            uint8_t* label = 
                                par->afl_->flags_->at(source_id).at(e_idx);
                            label[part_id >> 3] |= (1 << (part_id & 7));
                        };
                    dijkstra.apply_to_closed(fn_arcflags);
                    par->nprocessed_++;
                }
                return 0;
            };

            // OK, let's fork some threads
            const uint32_t NUM_THREADS = 4;
            pthread_t threads[NUM_THREADS];
            afl_params task_data[NUM_THREADS];

            for(uint32_t i = 0; i < NUM_THREADS; i++)
            {
                task_data[i].thread_id_ = i;
                task_data[i].max_threads_ = NUM_THREADS;
                task_data[i].afl_ = afl;
                task_data[i].nprocessed_ = 0;
                task_data[i].expander_ = fn_new_expander();

                pthread_create(&threads[i], NULL, 
                        thread_compute_fn, (void*) &task_data[i]);
            }
            std::cerr << "forked " << NUM_THREADS << " threads \n";

            while(true)
            {
                uint32_t nprocessed = 0;
                for(uint32_t i = 0; i < NUM_THREADS; i++)
                { nprocessed += task_data[i].nprocessed_; }

                std::cerr 
                    << "\rprogress: " << nprocessed 
                    << " / " << g->get_num_nodes();
                if(nprocessed == g->get_num_nodes()) { break; }
                else { sleep(2); }
            }

            for(uint32_t i = 0; i < NUM_THREADS; i++)
            { delete task_data[i].expander_; }

            std::cerr << "\nall done\n"<< std::endl;
            return afl;
        }

    private:
        // creation via ::compute or ::load only, please
        // @param g: the graph to be labelled
        // @param partitioning: a disjoint partitioning of nodes in @param g
        af_labelling(
                warthog::graph::planar_graph* g, 
                std::vector<uint32_t>* partitioning);

        std::vector<std::vector<uint8_t*>>* flags_;
        warthog::graph::planar_graph* g_;
        std::vector<uint32_t>* part_;
        uint32_t bytes_per_label_;

};

}

}

#endif

