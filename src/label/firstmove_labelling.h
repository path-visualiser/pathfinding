#ifndef WARTHOG_CPD_LABELLING_H
#define WARTHOG_CPD_LABELLING_H

// label/firstmove_labelling.h
//
// Computes a compressed first-move labelling for every node in a given
// graph. The compression scheme is run-length encoding. The column order
// is given by a depth-first search post-order traversal of the graph.
// For more details see:
//  
// @inproceedings{DBLP:conf/socs/StrasserHB14,
//   author    = {Ben Strasser and Daniel Harabor and Adi Botea},
//   title     = {Fast First-Move Queries through Run-Length Encoding},
//   booktitle = {Seventh Annual Symposium on Combinatorial Search}
//   year      = {2014} }
//
// @author: dharabor
// @created: 2017-12-21
// 

#include "bbaf_labelling.h"
#include "fch_expansion_policy.h"
#include "geom.h"
#include "solution.h"
#include "timer.h"

#include <fstream>
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

// limits on the maximum degree help to save memory 
// NB: if the max degree is not a power of two, we round up to 
// the next power
const uint32_t FM_MAX_DEGREE = 255;
const uint32_t FM_QWORD_SZ = sizeof(uint64_t);
const uint32_t FM_MAX_QWORDS =  
    std::max<int32_t>(1, 
        (FM_MAX_DEGREE & (FM_MAX_DEGREE-1) ? 
            (((FM_MAX_DEGREE & (FM_MAX_DEGREE-1))<<1) >> 6) : 
            (FM_MAX_DEGREE >> 6)) );

// we use run-length encoding to compress first-move data
struct fm_run
{
    uint32_t head_;
    uint8_t label_;

    void
    print(std::ostream& out)
    {
        out << " [" << head_ << ", " << label_ << "]";
    }
}

istream&
operator>>(fm_run& the_run, std::istream& in);

ostream&
operator<<(fm_run& the_run, std::ostream& out);

// a collection of optimal first moves. 
// this data structure is useful during preprocessing.
struct fm_coll
{
    fm()
    {
        for(uint32_t i = 0; i < FM_MAX_QWORDS; i++)
        { moves_[i] = 0; }
    }

    fm&
    operator=(const fm& other)
    {
        for(uint32_t i = 0; i < FM_MAX_QWORDS; i++)
        { moves_[i] = other.moves_[i]; }
        return *this;
    }

    fm&
    operator|=(const fm& other)
    {
        for(uint32_t i = 0; i < FM_MAX_QWORDS; i++)
        { moves_[i] |= other.moves_[i]; }
        return *this;
    }

    fm&
    operator&=(const fm& other)
    {
        for(uint32_t i = 0; i < FM_MAX_QWORDS; i++)
        { moves_[i] = moves_[i] & other.moves_[i]; }
        return retval;
    }

    void
    add_move(uint32_t move)
    {
        assert(move <= FM_MAX_DEGREE);
        uint32_t byte = move >> 6;
        uint32_t bit  = move & 63;
        moves_[byte] |= (1 << bit);
    }

    // return true if the collection has anyfirst moves, otherwise false
    bool
    eval() 
    {
        uint64_t retval = 0; 
        for(uint32_t i = 0; i < moves_[i]; i++)
        { retval |= moves_[i]; }
        return retval;
    }

    uint64_t moves_[FM_MAX_QWORDS];
}

class fm_labelling 
{
    friend std::ostream&
    operator<<(std::ostream& out, fm_labelling& lab);

    friend std::istream&
    warthog::label::operator>>(std::istream& in, fm_labelling& lab);

    public:

        ~fm_labelling();

        inline warthog::graph::planar_graph*
        get_graph() 
        { 
            return g_;
        }

        // @param node_id: the current node
        // @param target_id: the target node
        //
        // @return the index of the optimal move @param node_id to
        // @param target_id
        fm_label&
        get_label(uint32_t node_id, uint32_t uint32_t target_id)
        {
            assert(node_id < g_->get_num_nodes());
            // IMPLEMENT ME
        }

        inline size_t
        mem()
        {
            size_t retval = sizeof(this);
            retval += sizeof(lab_) * lab_->size();
            for(uint32_t i = 0; i < lab_->size(); i++)
            {
                retval += (sizeof(fm_label) * lab_->at(i).size();
            }
            return retval;
        }

        static warthog::label::fm_labelling*
        load(const char* filename, warthog::graph::planar_graph* g, 
            std::vector<uint32_t>* rank, std::vector<uint32_t>* part)
        {
            std::cerr << "loading fm_labelling from file " 
                << filename << std::endl;
            std::ifstream ifs(filename, 
                    std::ios_base::in|std::ios_base::binary);

            if(!ifs.good())
            {
                std::cerr << "load failed (no such file?)" << std::endl;
                ifs.close();
                return 0;
            }

            warthog::label::fm_labelling* lab = 
                new warthog::label::fm_labelling();

            ifs >> *lab;

            if(!ifs.good())
            {
                std::cerr << "load failed" << std::endl;
                delete lab;
                lab = 0;
            }

            ifs.close();
            return lab;
        }
        
        static void
        save(const char* filename, warthog::label::fm_labelling& lab)
        {
            std::cerr << "writing labels to file " << filename << "\n";
            std::ofstream out(filename, 
                    std::ios_base::out|std::ios_base::binary);

            out << lab;

            if(!out.good())
            {
                std::cerr << "\nerror trying to write to file " 
                    << filename << std::endl;
            }
            out.close();
        }

        // compute labels for all nodes specified by the given workload
        template <typename t_expander>
        static warthog::label::fm_labelling*
        compute(warthog::graph::planar_graph* g, 
                std::vector<uint32_t>* column_order,
                std::function<t_expander*(void)>& fn_new_expander,
                warthog::util::workload_manager* workload)
        {
            warthog::timer t;
            t.start();

            if(g == 0 || column_order == 0) { return 0; } 

            struct shared_data
            {
                std::function<t_expander*(void)> fn_new_expander_;
                warthog::label::fm_labelling* lab_;
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

                warthog::label::fm_labelling* lab = shared->lab_;
                warthog::util::workload_manager* workload = shared->workload_;

                // variable used to track the node currently being processed
                uint32_t source_id;

                // alocate memory for the first moves
                std::vector<fm> row(lab->g_->get_num_nodes());

                // callback function used to record the optimal first move 
                std::function<void(
                        warthog::search_node*, 
                        warthog::search_node*,
                        double, uint32_t)> on_generate_fn = 
                [&source_id, &row, lab]
                (warthog::search_node* succ, warthog::search_node* from,
                            double edge_cost, uint32_t edge_id) -> void
                {
                    if(from == 0) { return; } // start node 

                    if(from->get_id() == source_id) // start node successors
                    { 
                        assert(edge_id < 
                        lab->g_->get_node(source_id)->out_degree());
                        row.at(succ->get_id()).add_move(edge_id);
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
                        if(alt_g < g_val) { row.at(s_id) = row .at(f_id); }
                        
                        // add to the list of optimal first moves
                        if(alt_g == g_val) { row.at(s_id) |= row.at(f_id); }
                    }
                };

                compress_fn(std::vector<fm>& row, 
                        std::vector<fm_run>& rle_row)
                {
                    fm current = row.at(0);
                    for(uint32_t col = 0; col < row.size(); col++)
                    {
                        fm


                    } 
                }

                warthog::zero_heuristic h;
                std::shared_ptr<warthog::fch_expansion_policy> 
                    expander(shared->fn_new_expander_());
                warthog::flexible_astar 
                    <warthog::zero_heuristic, warthog::fch_expansion_policy>
                        dijk(&h, expander.get());
                dijk.apply_on_generate(on_generate_fn);

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
                    compress_fn(row, lab_->at(i));
                    par->nprocessed_++;
                }
                return 0;
            };

            compute_dfs_order(g, dfs_order, apex_id);

            warthog::label::fm_labelling* lab = 
                new warthog::label::fm_labelling(g, rank, part);

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
            lab->compute_fm_labels(workload); // single threaded
            t.stop();

            std::cerr 
                << "total preproc time (seconds): "
                << t.elapsed_time_micro() / 1000000 << "\n";

            return lab;
        }

    private:
        // only via ::compute or ::load please
        fm_labelling(
            warthog::graph::planar_graph* g, 
            std::vector<uint32_t>* rank,
            std::vector<uint32_t>* partitioning);

        // CPD-based preprocessing computes labels for every edge 
        // @param contraction order of every node in the graph
        void 
        compute_fm_labels(warthog::util::workload_manager* workload);

        // Computes a CPD post-order id for every node in the hierarchy
        // (top-down traversal)
        // @param id of the highest node in the contraction hierarchy
        void
        compute_firstmove_ids(uint32_t apex_id);

        void
        compress(std::vector< std::vector < uint32_t >>& fm);

        warthog::graph::planar_graph* g_;
        std::vector< std::vector< fm_label >>* lab_;
};

std::istream&
operator>>(std::istream& in, warthog::label::fm_labelling& lab);

std::ostream&
operator<<(std::ostream& in, warthog::label::fm_labelling& lab);


}


}

#endif

