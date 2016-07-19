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

#include "constants.h"
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

// An order-of-contraction is a sequential list of node-ids that indicates
// the order in which nodes should be contracted. e.g. 
// 5, 7, 4, 1, 3, 2, 6  (etc)
//
// This function takes an OOC and turns it into a ranked list, indexed by 
// node-id. Given a node-id, x, rank[x] indicates the position of node-id 
// x in the OOC; e.g. for the OOC above, rank[1] = 5
void
convert_order_of_contraction_to_ranked_list(std::vector<uint32_t>& ooc, 
        std::vector<uint32_t>& rank);

// compute the out-arc closure of @param source 
void
compute_closure(uint32_t source, warthog::graph::planar_graph* g, 
        std::set<uint32_t>* closure, uint32_t maxdepth=warthog::INF);
}

}

#endif

