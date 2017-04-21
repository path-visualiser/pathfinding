#ifndef WARTHOG_LABEL_BB_LABELLING_H
#define WARTHOG_LABEL_BB_LABELLING_H

// label/bb_labelling.h
//
// For every outgoing edge of every node in a graph we store a rectangular
// bounding box. Inside the box can be found all nodes that are reached 
// optimally by a path whose first edge is the edge at hand.
//
// For theory see:
//
// [Wager & Willhalm, 2005, 
// Geometric Containers for Efficient Shortest Path Computation, 
// Journal of Experimental Algorithms, vol 10, pp 1-30]
//
// @author: dharabor
// @created: 2017-04-18
//

#include "flexible_astar.h"
#include "geom.h"
#include "planar_graph.h"
#include "search_node.h"
#include "zero_heuristic.h"

#include <algorithm>
#include <functional>
#include <iostream>
#include <ostream>
#include <pthread.h>
#include <unistd.h>

namespace warthog
{

namespace label
{

class bb_labelling
{
    public:
        virtual ~bb_labelling();

        // print/serialise every label
        void
        print(std::ostream& out);

        inline warthog::geom::rectangle
        get_label(uint32_t node_id, uint32_t edge_id)
        {
            return labels_->at(node_id).at(edge_id);
        }

        static warthog::label::bb_labelling*
        load(const char* filename, warthog::graph::planar_graph* g);
        
        template <typename t_expander>
        static warthog::label::bb_labelling*
        compute(warthog::graph::planar_graph* g, 
                std::function<t_expander*(void)>& fn_new_expander)
        {
            if(g == 0) { return 0; } 
            std::cerr << "computing bb labels\n";

            warthog::label::bb_labelling* bbl = 
                new warthog::label::bb_labelling(g);

            // pthreads setup

            struct bbl_params
            {
                uint32_t thread_id_;
                uint32_t max_threads_;
                uint32_t nprocessed_;
                t_expander* expander_;
                warthog::label::bb_labelling* bbl_;
            };

            // the actual precompute function
            void*(*thread_compute_fn)(void*) = [] (void* args_in) -> void*
            {
                bbl_params* par = (bbl_params*) args_in;

                // need to keep track of the first edge on the way to the 
                // current node (the solution is a bit hacky as we break the 
                // chain of backpointers to achieve this; it doesn't matter, 
                // we don't care about the path)
                std::function<void(warthog::search_node*)> relax_fn = 
                    [] (warthog::search_node* n) -> void
                    {
                        // the start node and its children don't need their 
                        // parent pointers updated. for all other nodes the
                        // grandparent becomes the parent
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
                warthog::flexible_astar<warthog::zero_heuristic, t_expander>
                     dijkstra(&heuristic, par->expander_);
                dijkstra.apply_on_relax(relax_fn);

                for(uint32_t i = 0; i < par->bbl_->g_->get_num_nodes(); i++)
                {
                    // source nodes are evenly divided among all threads;
                    // skip any source nodes not intended for current thread
                    if((i % par->max_threads_) != par->thread_id_) { continue; }

                    // process the source
                    uint32_t source_id = i;
                    uint32_t ext_source_id = 
                        par->bbl_->g_->to_external_id(source_id);
                    warthog::graph::node* source = 
                        par->bbl_->g_->get_node(source_id);

                    warthog::problem_instance pi(ext_source_id, warthog::INF);
                    dijkstra.get_length(pi);
                    
                    // we need an easy way to convert between the ids of nodes
                    // adjacent to the source and the edge index to reach them
                    std::unordered_map<uint32_t, uint32_t> idmap;
                    uint32_t edge_idx = 0;
                    for(warthog::graph::edge_iter it = 
                            source->outgoing_begin(); 
                            it != source->outgoing_end(); it++)
                    {
                        idmap.insert( std::pair<uint32_t, uint32_t>(
                                    (*it).node_id_, edge_idx));
                        edge_idx++;
                    }

                    // compute the extent of each rectangle bounding box
                    std::function<void(warthog::search_node*)> bbox_fn = 
                    [par, &source_id, &idmap](warthog::search_node* n) 
                    -> void
                    {
                        // skip the source node
                        // (it doesn't belong to any rectangle)
                        assert(n);
                        if(n->get_id() == source_id) { return; } 
                        assert(n->get_parent());

                        std::unordered_map<uint32_t, uint32_t>::iterator it;

                        // successors of the source are special
                        if(n->get_parent()->get_parent() == 0)
                        { it = idmap.find(n->get_id()); }
                        else
                        // all the other nodes
                        { it = idmap.find(n->get_parent()->get_id()); }

                        if(it == idmap.end())
                        {
                            // sanity check
                            std::cerr << "err; invalid first edge; "
                                 << " n id: "<< n->get_id() << " parent id: " 
                                 << n->get_parent()->get_id() << std::endl;
                            exit(0);
                        }

                        // grow the rectangle
                        int32_t x, y;
                        par->bbl_->g_->get_xy(n->get_id(), x, y);
                        assert(x != warthog::INF && y != warthog::INF);
                        par->bbl_->labels_->at(source_id).at(
                                (*it).second).grow(x, y);
                        assert(par->bbl_->labels_->at(source_id).at(
                                    (*it).second).is_valid());
                    };
                    dijkstra.apply_to_closed(bbox_fn);
                }
                return 0;
            };

            // OK, let's fork some threads
            const uint32_t NUM_THREADS = 4;
            pthread_t threads[NUM_THREADS];
            bbl_params task_data[NUM_THREADS];

            // define workloads
            for(uint32_t i = 0; i < NUM_THREADS; i++)
            {
                task_data[i].thread_id_ = i;
                task_data[i].max_threads_ = NUM_THREADS;
                task_data[i].bbl_ = bbl;
                task_data[i].nprocessed_ = 0;
                task_data[i].expander_ = fn_new_expander();

                pthread_create(&threads[i], NULL, 
                        thread_compute_fn, (void*) &task_data[i]);
            }

            while(true)
            {
                uint32_t nprocessed = 0;
                for(uint32_t i = 0; i < NUM_THREADS; i++)
                { nprocessed += task_data[i].nprocessed_; }

                std::cerr 
                    << "\rprogress: " << nprocessed 
                    << " / " << g->get_num_nodes();
                if(nprocessed == g->get_num_nodes()) { break; }
                else { sleep(5); }
            }

            for(uint32_t i = 0; i < NUM_THREADS; i++)
            { delete task_data[i].expander_; }

            std::cerr << "\nall done\n"<< std::endl;
            return bbl;
        }

    private:
        // creation via ::compute and ::load only, please
        bb_labelling(warthog::graph::planar_graph* g);
        warthog::graph::planar_graph* g_;
        std::vector<std::vector<warthog::geom::rectangle>>* labels_;

};

}

}

#endif

