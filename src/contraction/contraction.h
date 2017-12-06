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
#include "graph.h"

#include <cstdint>
#include <fstream>
#include <set>
#include <vector>

namespace warthog
{

class dummy_filter;
class euclidean_heuristic;

template<typename FILTER>
class graph_expansion_policy;

template<typename H, typename E>
class flexible_astar;

namespace graph
{

class planar_graph;

}


namespace ch
{

// in the case of bidirectional search
enum search_direction { UP = 1, DOWN = 2, ANY = 3 };

// load a node order-of-contraction file. such files comprise a list of
// node ids in the order each node was contracted. thus, the first value
// in the file is the identifier of the first node to be contracted.
//
// @param filename: the file containing the order
// @param order: the container in which to load the node order
// @param lex_order: (default: false) when true, the contraction order is 
// converted into a lexical order. This is equivalent to calling 
// ::value_index_swap_dimacs after loading the contraction order
bool
load_node_order(const char* filename, 
        std::vector<uint32_t>& order, bool lex_order=false);

// The main algorithm for searching a contraction hierarhcy is bi-directional
// Dijkstra search (BCH). This algorithm traverses the hierarchy as follows:
// In the forward direction, it follows only up edges
// In the backward direction, it follows only incoming down edges
// In this function we prune away all other edges in order to speed up
// search (there's no longer a need to filter on the fly other types of edges)
//
// @param g: the graph to be optimised
// @param rank: the contraction order for the nodes of @param g
void
optimise_graph_for_bch(warthog::graph::planar_graph* g,
        std::vector<uint32_t>* rank);

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
//
// NB: min{vec[i]} elements are added to the front of @param vec
void
value_index_swap_dimacs(std::vector<uint32_t>& vec);

// compute the out-arc closure of @param source 
void
compute_closure(uint32_t source, warthog::graph::planar_graph* g, 
        std::set<uint32_t>* closure, uint32_t maxdepth=warthog::INF);

void
compute_down_closure(uint32_t source, 
        warthog::graph::planar_graph* g, std::vector<uint32_t>* rank,
        std::set<uint32_t>* closure);

// given an existing contraction hierarchy, create a partition
// of the nodes using a greedy bottom-up partitioning scheme.
// the result of the partitioning is a labeling of each node with
// with an integer in the range [0, @param nparts]. 
// the labeling is stored in @param part
void
partition_greedy_bottom_up(
        warthog::graph::planar_graph* g, 
        std::vector<uint32_t>* rank,
        uint32_t nparts, 
        std::vector<uint32_t>* part);

// unpack a shortcut edge by locating the two intermediate edges
// that are being bypassed
void
unpack_and_list_edges(warthog::graph::edge* scut,
        uint32_t scut_tail_id, warthog::graph::planar_graph* g,
        std::vector<warthog::graph::edge*>& unpacked, bool recurse=true);

// recursively unpack a shortcut edge and
// collect all the intermediate nodes, in order
void
unpack(uint32_t from_id,
        warthog::graph::edge_iter it_e,
        warthog::graph::planar_graph* g,
        std::set<uint32_t>& intermediate);

// sorts the list of successors such that all up successors appear first
// then all down successors. also stores the offset for the first down succ.
// this sorting improves FCH performance (no need to always check rank)
// NB: assumes max out degree (i.e. all up + down successors) is 256
void
fch_sort_successors(
        warthog::graph::planar_graph* g, 
        std::vector<uint32_t>* rank);

}

}

#endif

