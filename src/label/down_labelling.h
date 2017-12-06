#ifndef WARTHOG_DOWN_LABELLING_H
#define WARTHOG_DOWN_LABELLING_H

// label/down_labelling.h
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
    merge(id_range_label& other)
    {
        left = left < other.left ? left : other.left;
        right = right > other.right ? right : other.right;
    }

    bool
    contains(int32_t val)
    {
        return (right - val) + left >= left;
    }

    int32_t left, right;
};

struct down_label
{
    down_label(bbaf_label& bbaf, id_range_label& id) 
        : bbaf_(bbaf), ids_(id) { }

    bbaf_label& bbaf_;
    id_range_label& ids_;
};


class down_labelling 
{

    public:

        ~down_labelling();

        inline std::vector<uint32_t>*
        get_partitioning()
        {
            return lab_bbaf_->get_partitioning();
        }

        inline warthog::graph::planar_graph*
        get_graph() 
        { 
            return lab_bbaf_->get_graph();
        }

        down_label
        get_label(uint32_t node_id, uint32_t edge_idx)
        {
            assert(edge_idx < lab_idr_->at(node_id).size());
            //down_label ret(
            //        lab_bbaf_->get_label(node_id, edge_idx), 
            //        lab_idr_->at(node_id).at(edge_idx));
            bbaf_label dummy;
            down_label ret( dummy, lab_idr_->at(node_id).at(edge_idx));
            return ret;
        }

        void
        print(std::ostream& out) { /* IMPLEMENT ME */ }

        static warthog::label::down_labelling*
        load(const char* filename, warthog::graph::planar_graph* g, 
            std::vector<uint32_t>* partitioning) 
        { return 0; /* IMPLEMENT ME */ }

        // compute labels for all nodes specified by the given workload
        static warthog::label::down_labelling*
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

            warthog::label::down_labelling* lab = 
                new warthog::label::down_labelling();

            std::cerr << "computing id-range labels\n";
            lab->compute_id_range_labels(g, rank);

        //    std::cerr << "computing bbaf labels\n";
        //    lab->lab_bbaf_ = 
        //        warthog::label::bbaf_labelling::compute
        //        <warthog::fch_expansion_policy>
        //        (g, part, get_expander_fn, workload);
            
            return lab;
        }

        inline int32_t
        get_dfs_index(uint32_t graph_id) { return dfs_order_->at(graph_id); }

        inline size_t
        mem()
        {
            size_t retval = sizeof(this);
            for(uint32_t i = 0; i < lab_idr_->size(); i++)
            {
                retval += sizeof(id_range_label)* lab_idr_->at(i).size();
            }
            retval += sizeof(int32_t) * dfs_order_->size();
            retval += lab_bbaf_->mem();
            return retval;
        }


    private:
        // only via ::compute or ::load please
        down_labelling();

        void 
        compute_id_range_labels(
                warthog::graph::planar_graph*, std::vector<uint32_t>*);

        std::vector<int32_t>* dfs_order_;
        warthog::label::bbaf_labelling* lab_bbaf_;
        std::vector< std::vector< id_range_label >>* lab_idr_;
};

}

}

#endif

