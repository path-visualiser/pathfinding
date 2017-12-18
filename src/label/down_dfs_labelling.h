#ifndef WARTHOG_DOWN_DFS_LABELLING_H
#define WARTHOG_DOWN_DFS_LABELLING_H

// label/down_dfs_labelling.h
//
// Collects a variety of different labels that we can compute for 
// the down-closure of a node. Intended for use with forward-driven 
// contraction hierarchies.
//
// Current edge labels:
//  - id range: specifying the range of ids in the down closure of each edge
//
//  - bounding box: a rectangular bounding box that contains every
//  node in the down closure of each edge
//  
//  - arc flags: a bitfield that indicates set-reachability info for every
//  node in the down closure of each edge
//
// @author: dharabor
// @created: 2017-12-06
// 

#include "bbaf_labelling.h"
#include "fch_expansion_policy.h"
#include "geom.h"
#include "solution.h"

#include <vector>
#include <cstdint>

namespace warthog
{

namespace graph
{

class planar_graph;
class edge;

} 

namespace label
{

struct id_range_label
{
    id_range_label() { left = INT32_MAX; right = INT32_MIN; } 

    id_range_label&
    operator=(const id_range_label& other)
    {
        left = other.left; right = other.right; 
        return *this;
    }

    void
    grow(int32_t val)
    {
        left = val < left  ? val : left;
        right = val > right ? val : right;
    }

    void
    grow(const id_range_label& other)
    {
        left = left < other.left ? left : other.left;
        right = right > other.right ? right : other.right;
    }

    inline bool
    contains(int32_t val)
    {
        return (right - val) + left >= left;
    }

    int32_t left, right;
};

struct down_dfs_label
{
    down_dfs_label(size_t arcflag_bytes)
    {
        for(uint32_t i = 0; i < arcflag_bytes; i++) { flags_.push_back(0); }
    }

    down_dfs_label&
    operator=(const down_dfs_label& other)
    {
        rank_ = other.rank_;
        ids_ = other.ids_;
        bbox_ = other.bbox_;
        for(uint32_t i = 0; i < other.flags_.size(); i++)
        { flags_[i] = other.flags_.at(i); }
        return *this;
    }

    void
    merge(const down_dfs_label& other)
    {
        rank_.grow(other.rank_);
        ids_.grow(other.ids_);
        bbox_.grow(other.bbox_);
        for(uint32_t i = 0; i < flags_.size(); i++)
        { flags_[i] |= other.flags_[i]; }
    }

    id_range_label rank_;
    id_range_label ids_;
    std::vector<uint8_t> flags_;
    warthog::geom::rectangle bbox_;
};


class down_dfs_labelling 
{

    public:

        ~down_dfs_labelling();

        inline std::vector<uint32_t>*
        get_partitioning()
        {
            return part_;
        }

        inline warthog::graph::planar_graph*
        get_graph() 
        { 
            return g_;
        }

        down_dfs_label&
        get_label(uint32_t node_id, uint32_t edge_idx)
        {
            assert(edge_idx < lab_->at(node_id).size());
            return lab_->at(node_id).at(edge_idx);
        }

        void
        print(std::ostream& out) { /* IMPLEMENT ME */ }

        static warthog::label::down_dfs_labelling*
        load(const char* filename, warthog::graph::planar_graph* g, 
            std::vector<uint32_t>* partitioning) 
        { return 0; /* IMPLEMENT ME */ }

        // compute labels for all nodes specified by the given workload
        static warthog::label::down_dfs_labelling*
        compute(warthog::graph::planar_graph* g, 
                std::vector<uint32_t>* part, 
                std::vector<uint32_t>* rank,
                bool sort_successors = true)
        {
            if(g == 0 || part  == 0) { return 0; } 

            // sort the edges of the graph s.t. all up successors of a node
            // appear before any down successor
            if(sort_successors)
            { warthog::ch::fch_sort_successors(g, rank); }

            warthog::label::down_dfs_labelling* lab = 
                new warthog::label::down_dfs_labelling(g, part);

            std::cerr << "computing dfs labels\n";
            lab->compute_labels(rank); // single threaded

            std::cerr << "computing dijkstra labels\n";

            // select some nodes and reset their labels
            warthog::util::workload_manager workload(g->get_num_nodes());
            for(uint32_t i = 0; i < g->get_num_nodes(); i++)
            {
                warthog::graph::node* n = g->get_node(i);

                // nodes with degree >= 100
                //if(n->out_degree() >= 100)
                // top 1% highest nodes
                if(rank->at(i) >= (uint32_t)(rank->size()*0.99))
                { 
                    workload.set_flag(i, true);

                    // throw away the dfs labels
                    warthog::graph::edge_iter begin = n->outgoing_begin();
                    warthog::graph::edge_iter end = n->outgoing_end();
                    for(warthog::graph::edge_iter it = begin; it != end; it++)
                    {
                        lab->lab_->at(i).at(it - begin) = 
                            down_dfs_label(lab->bytes_per_af_label_);
                    }
                }
            }

            std::function<fch_expansion_policy*(void)> expander_fn = 
            [g, rank] (void) -> fch_expansion_policy*
            { return new warthog::fch_expansion_policy(g, rank); };

            struct shared_data
            {
                std::function<fch_expansion_policy*(void)> fn_new_expander_;
                warthog::label::down_dfs_labelling* lab_;
                warthog::util::workload_manager* workload_;
                std::vector<uint32_t>* rank_;
            };

            // The actual precompute function. We construct a 
            // Dijkstra-based preprocessing to improve the labels for 
            // selected sets of nodes
            void*(*thread_compute_fn)(void*) = 
            [] (void* args_in) -> void*
            {
                warthog::helpers::thread_params* par = 
                    (warthog::helpers::thread_params*) args_in;
                shared_data* shared = (shared_data*) par->shared_;

                warthog::label::down_dfs_labelling* lab = shared->lab_;
                std::vector<uint32_t>* rank = shared->rank_;
                warthog::util::workload_manager* workload = shared->workload_;

                // variable used to track the node currently being processed
                uint32_t source_id;
                
                // alocate memory for the first moves
                std::vector<uint32_t> first_move(lab->g_->get_num_nodes());

                // callback function used to record the optimal first move 
                std::function<void(
                        warthog::search_node*, 
                        warthog::search_node*,
                        double, uint32_t)> on_generate_fn = 
                [&source_id, &first_move, lab]
                (warthog::search_node* succ, warthog::search_node* from,
                            double edge_cost, uint32_t edge_id) -> void
                {
                    if(from == 0) { return; } // start node 

                    if(from->get_id() == source_id) // start node successors
                    { 
                        assert(edge_id < 
                        lab->g_->get_node(source_id)->out_degree());
                        first_move.at(succ->get_id()) = edge_id; 
                    }
                    else // all other nodes
                    {
                        uint32_t s_id = succ->get_id();
                        uint32_t f_id = from->get_id();
                        double alt_g = from->get_g() + edge_cost;
                        double g_val = 
                            succ->get_search_id() == from->get_search_id() ? 
                            succ->get_g() : warthog::INF; 

                        assert(first_move.at(f_id) < 
                        lab->g_->get_node(source_id)->out_degree());

                        //  update first move
                        if(alt_g < g_val) 
                        { first_move.at(s_id) = first_move.at(f_id); }
                    }
                };


                std::function<void(warthog::search_node*)> on_expand_fn =
                [&source_id, &first_move, rank, lab]
                (warthog::search_node* current) -> void
                {
                    if(current->get_id() == source_id) { return; }

                    uint32_t node_id = current->get_id();
                    assert(node_id < first_move.size());

                    uint32_t edge_idx = first_move.at(node_id);
                    assert(edge_idx < lab->lab_->at(source_id).size());
                    down_dfs_label& s_lab = 
                        lab->lab_->at(source_id).at(edge_idx);

                    s_lab.rank_.grow(rank->at(node_id));
                    s_lab.ids_.grow(lab->dfs_order_->at(node_id));

                    int32_t x, y;
                    lab->g_->get_xy(node_id, x, y);
                    s_lab.bbox_.grow(x, y);

                    uint32_t s_part = lab->part_->at(node_id);
                    s_lab.flags_[s_part >> 3] |= (1 << (s_part & 7)); // div8, mod8
                    assert(s_lab.flags_[s_part >> 3] & (1 << (s_part & 7)));
                };

                warthog::zero_heuristic h;
                std::shared_ptr<warthog::fch_expansion_policy> 
                    expander(shared->fn_new_expander_());
                warthog::flexible_astar 
                    <warthog::zero_heuristic, warthog::fch_expansion_policy>
                        dijk(&h, expander.get());
                dijk.apply_on_generate(on_generate_fn);
                dijk.apply_on_expand(on_expand_fn);

                for(uint32_t i = 0; i < lab->g_->get_num_nodes(); i++)
                {
                    // skip any nodes not part of the precomputation workload
                    if(!workload->get_flag(i))
                    { continue; }

                    // source nodes are evenly divided among all threads;
                    // skip any source nodes not intended for current thread
                    if((i % par->max_threads_) != par->thread_id_) 
                    { continue; }

                    source_id = i;
                    uint32_t ext_source_id = 
                        lab->g_->to_external_id(source_id);
                    warthog::problem_instance problem(ext_source_id, 
                            warthog::INF);
                    //problem.verbose_ = true;
                    warthog::solution sol;
                    dijk.get_path(problem, sol);
                    par->nprocessed_++;
                }

                return 0;
            };

            shared_data shared;
            shared.fn_new_expander_ = expander_fn;
            shared.lab_ = lab;
            shared.workload_ = &workload;
            shared.lab_ = lab;
            shared.rank_ = rank;

            warthog::helpers::parallel_compute(
                    thread_compute_fn, &shared, 
                    workload.num_flags_set());
            std::cerr << "\nall done\n"<< std::endl;
            return lab;
        }

        inline int32_t
        get_dfs_index(uint32_t graph_id) { return dfs_order_->at(graph_id); }

        inline size_t
        mem()
        {
            size_t retval = sizeof(this);
            for(uint32_t i = 0; i < lab_->size(); i++)
            {
                retval += (sizeof(down_dfs_label) + bytes_per_af_label_) 
                    * lab_->at(i).size();
            }
            retval += sizeof(int32_t) * dfs_order_->size();
            return retval;
        }


    private:
        // only via ::compute or ::load please
        down_dfs_labelling(
            warthog::graph::planar_graph* g, 
            std::vector<uint32_t>* partitioning);

        // DFS-based preprocessing computes labels for every edge 
        // @param contraction order of every node in the graph
        void 
        compute_labels(std::vector<uint32_t>*);

        warthog::graph::planar_graph* g_;
        std::vector<uint32_t>* part_;
        size_t bytes_per_af_label_;

        std::vector<int32_t>* dfs_order_;
        std::vector< std::vector< down_dfs_label >>* lab_;

};

}

}

#endif

