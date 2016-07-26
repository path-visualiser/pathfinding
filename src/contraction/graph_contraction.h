#ifndef WARTHOG_GRAPH_CONTRACTION_H
#define WARTHOG_GRAPH_CONTRACTION_H

// graph_contraction.h
//
// abstract class from which different graph contraction
// schemes can be cooked up.
// 
// @author: dharabor
// @created: 2016-06-15
//

#include <cstring> // size_t
#include <cstdint>

namespace warthog
{

//class apriori_filter;
   
namespace graph
{

class planar_graph;

} 

namespace ch
{

class graph_contraction
{
    public:
        graph_contraction(warthog::graph::planar_graph* g);
        virtual ~graph_contraction();

        // contract a graph according to the node
        // order induced by ::next
        void
        contract();

//        inline warthog::apriori_filter*
//        get_filter() { return filter_; }

        inline warthog::graph::planar_graph*
        get_graph() { return g_; }
       
        // helper functions
        void
        set_verbose(bool verbose) { this->verbose_ = verbose; }

        bool 
        get_verbose() { return verbose_; } 

        virtual size_t
        mem();

    protected:

        // the method ::contract calls this function first
        // before doing anything else
        virtual void
        preliminaries() = 0;

        virtual void
        postliminaries() = 0;
        
        // return the id of the node that should be contracted next
        // and the value warthog::INF if all nodes have been contracted
        virtual uint32_t
        next() = 0;
        
        // return the result (path length) of a witness
        // search between @param from and @param to.
        virtual double
        witness_search(uint32_t from_id, uint32_t to_id, 
                       double via_len) = 0;

    private:
        bool verbose_;
        bool done_;
//        warthog::apriori_filter* filter_;
        warthog::graph::planar_graph* g_;
};


}

}

#endif
