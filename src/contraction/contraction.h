#ifndef WARTHOG_CONTRACTION_H
#define WARTHOG_CONTRACTION_H

// contraction.h
//
// Applies a contraction operation to each node in a 
// graph.
//
// Contraction in this case refers to a localised graph
// operation where additional ``shortcut'' edges are added to 
// the graph in such a way that the ``contracted'' node can
// be bypassed entirely during search (i.e. its neighbours are
// connected directly).
//
// In a contraction hierarchy every node is assigned a specific
// ``level'' which is equal to the contraction priority of the
// node. The highest priority nodes (== lowest level) are 
// contracted first and the process continues until every node
// has been contracted.
//
// For more details see:
// [Geisbergerger, Sanders, Schultes and Delling. 
// Contraction Hierarchies: Faster and Simpler Hierarchical 
// Routing in Road Networks. In Proceedings of the 2008
// Workshop on Experimental Algorithms (WEA)]
//
// @author: dharabor
// @created: 2016-01-25
//

#include <cstdint>
#include <fstream>
#include <vector>

namespace warthog
{

namespace graph
{
    class planar_graph;
}

namespace contraction
{
    // contract a graph using a specified node ordering 
    void 
    contract(warthog::graph::planar_graph& g, std::vector<uint32_t>& order, bool verbose = false);

    // the simplest contraction node order
    void
    make_input_order(warthog::graph::planar_graph& g, std::vector<uint32_t>& order);

    // helper functions
    void 
    load_node_order(const char* filename, std::vector<uint32_t>& order);
    void
    write_node_order(const char* filename, std::vector<uint32_t>& order);
}

}

#endif

