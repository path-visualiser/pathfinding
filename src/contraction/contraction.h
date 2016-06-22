#ifndef WARTHOG_CONTRACTION_H
#define WARTHOG_CONTRACTION_H

// contraction.h
//
// helper functions for working with contraction
// hierarchies.
//
// @author: dharabor
// @created: 2016-06-14
//

#include <cstdint>
#include <fstream>
#include <set>
#include <vector>

namespace warthog
{

class euclidean_heuristic;
class graph_expansion_policy;
template<typename H, typename E>
class flexible_astar;

namespace graph
{
class planar_graph;
}

namespace ch
{

typedef warthog::flexible_astar<
            warthog::euclidean_heuristic, 
            warthog::graph_expansion_policy> euc_astar;

void 
load_node_order(const char* filename, std::vector<uint32_t>& order);

void
write_node_order(const char* filename, std::vector<uint32_t>& order);

void
make_input_order(warthog::graph::planar_graph& g, std::vector<uint32_t>& order);

// compute the out-arc closure of @param source 
void
compute_closure(uint32_t source, warthog::graph::planar_graph* g, 
        std::set<uint32_t>* closure);
}

}

#endif

