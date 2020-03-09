#include "domains/xy_graph.h"
#include "heuristics/zero_heuristic.h"
#include "label/dfs_labelling.h"
#include "contraction/fch_expansion_policy.h"
#include "search/flexible_astar.h"
#include "search/graph_expansion_policy.h"
#include "search/problem_instance.h"
#include "search/search_node.h"
#include "search/solution.h"
#include "util/workload_manager.h"

#include <algorithm>

std::istream&
warthog::label::operator>>(std::istream& in, warthog::label::dfs_label& label)
{
    in >> label.bbox_;
    return in;
}

std::ostream&
warthog::label::operator<<(std::ostream& out, warthog::label::dfs_label& label)
{
    out << label.bbox_;
    return out;
}

std::istream&
warthog::label::operator>>(std::istream& in, warthog::label::dfs_labelling& lab)
{
    for(uint32_t n_id = 0; n_id < lab.g_->get_num_nodes(); n_id++)
    {
        warthog::graph::node* n = lab.g_->get_node(n_id);
        warthog::graph::edge_iter begin = n->outgoing_begin();
        warthog::graph::edge_iter end = n->outgoing_end();
        for(warthog::graph::edge_iter it = begin; it != end; it++)
        {
            in >> lab.lab_->at(n_id).at((size_t)(it - begin));
            assert(in.good());

            if(!in.good())
            {
                std::cerr << "unexpected error while reading labels\n";
                std::cerr 
                    << "[debug info] node: " << n_id
                    << " out-edge-index: " << (it - begin) << "\n";
                return in;
            }
        }
    }
    return in;
}

std::ostream&
warthog::label::operator<<(std::ostream& out,warthog::label::dfs_labelling& lab)
{
    for(uint32_t n_id = 0; n_id < lab.g_->get_num_nodes(); n_id++)
    {
        warthog::graph::node* n = lab.g_->get_node(n_id);
        warthog::graph::edge_iter begin = n->outgoing_begin();
        warthog::graph::edge_iter end = n->outgoing_end();
        for(warthog::graph::edge_iter it = begin; it != end; it++)
        {
            out << lab.lab_->at(n_id).at((size_t)(it - begin));
            if(!out.good())
            {
                std::cerr << "unexpected error while writing labels\n";
                std::cerr 
                    << "[debug info] node: " << n_id
                    << " out-edge-index: " << (it-begin) << "\n";
                return out;
            }
        }
    }
    return out;
}

warthog::label::dfs_labelling::dfs_labelling(
        warthog::graph::xy_graph* g, 
        std::vector<uint32_t>* rank)
    : g_(g), rank_(rank)
{
    dfs_order_ = new std::vector< uint32_t >();
    apex_id_ = compute_dfs_postorder_ids_ch(g_, rank_, dfs_order_);
   
    // allocate memory for edge labels
    lab_ = new std::vector< std::vector < dfs_label > >();
    lab_->resize(g_->get_num_nodes());
    dfs_label dummy;
    for(uint32_t n_id = 0; n_id < g_->get_num_nodes(); n_id++)
    {
        warthog::graph::node* n = this->g_->get_node(n_id);
        lab_->at(n_id).resize(n->out_degree(), dummy);
    }   
}

warthog::label::dfs_labelling::~dfs_labelling()
{ 
    delete lab_;
    delete dfs_order_;
}

void
warthog::label::dfs_labelling::compute_dfs_labels(
        warthog::util::workload_manager* workload)
{
    // Here we:
    // 1. perform a downwards dfs traversal of the CH from a given source node
    // 2. label the nodes in the subtree using a post-order visitation scheme
    // 3. store for every down edge a label that bounds the reachable subtree
    // 4. store for all nodes a label to bound its down closure
    std::vector< dfs_label > closure(
            this->g_->get_num_nodes(), dfs_label());
    std::vector< uint8_t > recurse(this->g_->get_num_nodes(), 1);

    std::function<void(uint32_t)> label_fn = 
        [this, workload, &recurse, &closure, &label_fn] 
        (uint32_t source_id) -> void
        {
            dfs_label& s_lab = closure.at(source_id);
            warthog::graph::node* source = this->g_->get_node(source_id);
            warthog::graph::edge_iter begin = source->outgoing_begin();

            for( warthog::graph::edge_iter it = begin; 
                    it != source->outgoing_end();
                    it++)
            {
                // skip up edges
                if(this->rank_->at(it->node_id_) > this->rank_->at(source_id)) 
                { continue; }

                // DFS
                if(recurse.at(it->node_id_))
                { label_fn(it->node_id_); }

                // grow the label of the down edge at hand
                if(workload->get_flag(source_id))
                {
                    lab_->at(source_id).at((size_t)(it - begin)) = 
                        closure.at(it->node_id_);
                }
                s_lab.merge(lab_->at(source_id).at((size_t)(it - begin)));
            }
            recurse.at(source_id) = 0;

            int32_t x, y;
            this->g_->get_xy(source_id, x, y);
            s_lab.bbox_.grow(x, y);

            //uint32_t s_part = this->part_->at(source_id);
            //s_lab.flags_[s_part >> 3] |= (1 << (s_part & 7)); // div8, mod8
            //assert(s_lab.flags_[s_part >> 3] & (1 << (s_part & 7)));
        };

    // Here we:
    // 1. perform an upwards dfs traversal of the CH from a given source node
    // 2. store for every up edge a label that bounds the up-reachable subtree
    // 3. store for all nodes a label to bound its up-reachable closure
    std::vector< dfs_label > up_closure(
            this->g_->get_num_nodes(), dfs_label());
    std::function<void(uint32_t)> up_label_fn = 
        [this, workload, &closure, &up_closure, &recurse, &up_label_fn] 
        (uint32_t source_id) -> void
        {
            warthog::graph::node* source = this->g_->get_node(source_id);
            warthog::graph::edge_iter begin = source->outgoing_begin();
            warthog::graph::edge_iter end = source->outgoing_end();
            dfs_label& s_lab = up_closure.at(source_id);

            // compute labels for each up edge by taking the union of:
            // 1. the up-closure of every up-edge
            // 2. the down-closure of every node in the up-closure
            for( warthog::graph::edge_iter it = begin; it != end; it++)
            {
                // skip down edges
                if(this->rank_->at(it->node_id_) < this->rank_->at(source_id)) 
                { continue; }

                // DFS
                if(recurse.at(it->node_id_))
                { up_label_fn(it->node_id_); }

                // grow the label for the up edge at hand
                if(workload->get_flag(source_id))
                {
                    dfs_label& e_lab = lab_->at(source_id).at((size_t)(it-begin));
                    e_lab.merge(up_closure.at(it->node_id_)); 
                }

                // grow the up-reachable closure for the source node
                s_lab.merge(lab_->at(source_id).at((size_t)(it-begin)));
            }

            recurse.at(source_id) = 0;
            s_lab.merge(closure.at(source_id));
        };

    // Top-down DFS from the apex to compute down-edge labels
    label_fn(apex_id_);
    recurse.clear();
    recurse.assign(this->g_->get_num_nodes(), 1);

    // bottom-up DFS to compute up-edge labels
    for(uint32_t n_id = 0; n_id < this->g_->get_num_nodes(); n_id++)
    { 
        up_label_fn(n_id); 
    }
}

// Computes a DFS post-order id for every node in a contraction
// hierarchy (i.e. a top-down traversal)
// @param id of the highest node in the contraction hierarchy
uint32_t
warthog::label::dfs_labelling::compute_dfs_postorder_ids_ch(
        warthog::graph::xy_graph* g, 
        std::vector<uint32_t>* rank, 
        std::vector<uint32_t>* dfs_ids)
{
    // find the apex of the hierarchy
    uint32_t apex_id = 0;
    for(uint32_t i = 0; i < rank->size(); i++)
    { 
        if(rank->at(i) > rank->at(apex_id)) 
        { apex_id = i; } 
    }

    uint32_t next_id = 0;
    dfs_ids->resize(g->get_num_nodes(), INT32_MAX);
    std::function<void(uint32_t)> dfs_id_fn =
    [dfs_ids, rank, &next_id, &dfs_id_fn, g] (uint32_t source_id) 
    -> void
    {
        warthog::graph::node* source = g->get_node(source_id);
        warthog::graph::edge_iter begin = source->outgoing_begin();
        for( warthog::graph::edge_iter it = begin; 
                it != source->outgoing_end();
                it++)
        {
            // skip up edges
            if(rank->at(it->node_id_) > rank->at(source_id)) 
            { continue; }

            // recurse
            if(dfs_ids->at(it->node_id_) == INT32_MAX)
            { dfs_id_fn(it->node_id_); }
        }
        if(dfs_ids->at(source_id) == INT32_MAX)
        { dfs_ids->at(source_id) = next_id++; }
    };

    // gogogo
    dfs_id_fn(apex_id);
    return apex_id;
}

struct shared_data
{
    std::function<warthog::fch_expansion_policy*(void)> fn_new_expander_;
    warthog::label::dfs_labelling* lab_;
    warthog::util::workload_manager* workload_;
    std::vector<uint32_t>* rank_;
};

struct search_listener
{
    std::vector<uint32_t>* first_move;
    warthog::label::dfs_labelling* lab; 
    uint32_t source_id;

   

    inline void
    generate_node(warthog::search_node* from, 
                  warthog::search_node* succ, 
                  warthog::cost_t edge_cost,
                  uint32_t edge_id) 
    {
        if(from == 0) { return; } // start node 

        if((uint32_t)from->get_id() == source_id) // start node successors
        { 
            assert(edge_id < 
            lab->get_graph()->get_node(source_id)->out_degree());
            first_move->at((uint32_t)succ->get_id()) = edge_id; 
        }
        else // all other nodes
        {
            uint32_t s_id = (uint32_t)succ->get_id();
            uint32_t f_id = (uint32_t)from->get_id();
            double alt_g = from->get_g() + edge_cost;
            double g_val = 
                succ->get_search_number() == from->get_search_number() ? 
                succ->get_g() : warthog::INF32; 

            assert(first_move->at(f_id) < 
            lab->get_graph()->get_node(source_id)->out_degree());

            //  update first move
            if(alt_g < g_val) 
            { first_move->at(s_id) = first_move->at(f_id); }
        }
    };

    inline void
    expand_node(warthog::search_node* current)
    {
        if((uint32_t)current->get_id() == source_id) { return; }

        uint32_t node_id = (uint32_t)current->get_id();
        assert(node_id < first_move->size());

        uint32_t edge_idx = first_move->at(node_id);
        warthog::label::dfs_label& s_lab = lab->get_label(source_id, edge_idx);

        int32_t x, y;
        lab->get_graph()->get_xy(node_id, x, y);
        s_lab.bbox_.grow(x, y);
    };

    inline void
    relax_node(warthog::search_node* current) { }
};

// compute labels for all nodes specified by the given workload
warthog::label::dfs_labelling*
warthog::label::dfs_labelling::precompute(
        warthog::graph::xy_graph* g, std::vector<uint32_t>* rank,
        warthog::util::workload_manager* workload)
{
    warthog::timer t;
    t.start();

    if(g == 0 || rank == 0 ) { return 0; } 

    std::function<fch_expansion_policy*(void)> expander_fn = 
    [g, rank] (void) -> fch_expansion_policy*
    { return new warthog::fch_expansion_policy(g, rank); };


    // The actual precompute function. We construct a 
    // Dijkstra-based preprocessing to improve the labels for 
    // selected sets of nodes
    void*(*thread_compute_fn)(void*) = 
    [] (void* args_in) -> void*
    {
        warthog::helpers::thread_params* par = 
            (warthog::helpers::thread_params*) args_in;
        shared_data* shared = (shared_data*) par->shared_;

        warthog::label::dfs_labelling* lab = shared->lab_;
        warthog::util::workload_manager* workload = shared->workload_;

        // variable used to track the node currently being processed
        uint32_t source_id;
        std::vector<uint32_t> first_move(lab->g_->get_num_nodes());


        warthog::zero_heuristic h;
        warthog::pqueue_min open;
        search_listener listener;
        std::shared_ptr<warthog::fch_expansion_policy> 
            expander(shared->fn_new_expander_());

        warthog::flexible_astar 
            <warthog::zero_heuristic, 
            warthog::fch_expansion_policy,
            warthog::pqueue_min,
            search_listener>
                dijk(&h, expander.get(), &open, &listener);

        listener.first_move = &first_move;
        listener.lab = lab;
        for(uint32_t i = 0; i < workload->num_flags_set(); i++)
        {
            // skip any nodes not part of the precomputation workload
            if(!workload->get_flag(i))
            { continue; }

            // source nodes are evenly divided among all threads;
            // skip any source nodes not intended for current thread
            if((i % par->max_threads_) != par->thread_id_) 
            { continue; }

            source_id = i;
            listener.source_id = source_id;

            uint32_t ext_source_id = 
                lab->g_->to_external_id(source_id);
            warthog::problem_instance problem(ext_source_id, 
                    warthog::INF32);
            //problem.verbose_ = true;
            warthog::solution sol;
            dijk.get_path(problem, sol);
            par->nprocessed_++;
        }
        return 0;
    };

    warthog::label::dfs_labelling* lab = 
        new warthog::label::dfs_labelling(g, rank);

    shared_data shared;
    shared.fn_new_expander_ = expander_fn;
    shared.lab_ = lab;
    shared.workload_ = workload;
    shared.lab_ = lab;

    std::cerr << "computing dijkstra labels\n";
    warthog::helpers::parallel_compute(
            thread_compute_fn, &shared, 
            workload->num_flags_set());

    std::cerr << "computing dfs labels...\n";
    workload->set_all_flags_complement();
    lab->compute_dfs_labels(workload); // single threaded
    t.stop();

    std::cerr 
        << "total preproc time (seconds): "
        << t.elapsed_time_micro() / 1000000 << "\n";

    return lab;
}
