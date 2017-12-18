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
                warthog::util::workload_manager* workload,
                bool sort_successors = true)
        {
            if(g == 0 || part  == 0) { return 0; } 

            // sort the edges of the graph s.t. all up successors of a node
            // appear before any down successor
            if(sort_successors)
            { warthog::ch::fch_sort_successors(g, rank); }

            std::function<fch_expansion_policy*(void)> get_expander_fn = 
            [g, rank] (void) -> fch_expansion_policy*
            {
                return new warthog::fch_expansion_policy(g, rank, 
                        warthog::ch::DOWN, false);
            };

            warthog::label::down_dfs_labelling* lab = 
                new warthog::label::down_dfs_labelling(g, part);

            std::cerr << "computing labels\n";
            lab->compute_labels(rank);
            //lab->improve_labels(rank);

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

        // Dijkstra-based preprocessing improves the labels for 
        // selected sets of nodes
        // @param contraction order of every node in the graph
        void 
        improve_labels(std::vector<uint32_t>*);

        warthog::graph::planar_graph* g_;
        std::vector<uint32_t>* part_;
        size_t bytes_per_af_label_;

        std::vector<int32_t>* dfs_order_;
        std::vector< std::vector< down_dfs_label >>* lab_;

};

}

}

#endif

