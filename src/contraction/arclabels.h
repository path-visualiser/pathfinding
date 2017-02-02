#ifndef WARTHOG_ARCLABELS_H
#define WARTHOG_ARCLABELS_H

// contraction/arclabels.h
//
// A collection of functions that compute arclabels for different types of 
// graphs and using different types of expansion policies.
// 
// @author: dharabor
// @created: 2017-01-28
//

#include "geom.h"

#include <ostream>
#include <vector>

namespace warthog
{

namespace graph
{

class planar_graph;
class corner_point_graph;

}

namespace arclabels
{

////////////////////////// GEOMETRIC CONTAINERS ////////////////////////////
//                (we use rectangular bounding boxes)

typedef std::vector<std::vector<warthog::geom::rectangle>> t_arclabels_bb;

// vanilla bounding-box arclabels:
// for every outgoing edge of every node in a graph we store a rectangular
// bounding box. Inside the box can be found all nodes that are reached 
// optimally by a path whose first edge is the edge at hand.
void
bb_compute(warthog::graph::planar_graph* g, std::ostream& out);

// bb arclabels computed on a graph that is a contraction hierarchy:
// In a CH edges are either "up" or "down". Nodes reached via an up edge
// have all of their successors generated and relaxed. Nodes reached via
// a down edge only generate and relax those successors further down 
// in the hierarchy.
void
ch_bb_compute(warthog::graph::planar_graph* g, 
        std::vector<uint32_t>* rank, std::ostream& out);

// bb arclabels that are: 
//  - computed on a corner point graph
//  - the graph is also a contraction hierarchy 
//  - the set of successors are pruned ala Jump Point Search
void
ch_bb_jpg_compute( 
        warthog::graph::corner_point_graph* g, 
        std::vector<uint32_t>* rank, std::ostream& out);

// generic helper function for printing computed labels
void
bb_print( t_arclabels_bb& labels, 
        warthog::graph::planar_graph* g_, std::ostream& out);

/////////////////////////////// ARC FLAGS //////////////////////////////
// NB: the idea with all compute functions is that a labeling is computed
// then printed to file as part of a preprocessing step. the file should
// be loaded, right before commencing search

typedef std::vector<std::vector<uint8_t*>> t_arclabels_af;

struct af_params
{
    uint32_t bytes_per_label_;
    uint32_t nparts_;
};

// helper function to figure out how many partitions there are
// and how many bytes are needed to store a corresponding label
warthog::arclabels::af_params
get_af_params(std::vector<uint32_t>* part);

// helper function to print an arcflag labeling
void
af_print( t_arclabels_af& flags_, 
        warthog::arclabels::af_params params, std::ostream& out);

// vanilla arcflags:
// partition the nodes of the graph into disjoint sets and store with 
// every outgoing edge of every node a bitfield where the number of bits
// is equal to the number of disjoint sets in the partition. the individual
// bits are set to true if the edge at hand appears on any optimal path
// to any node in its associated partition and false if it fails to appear
// on any optimal path to a node in its associated partition.
warthog::arclabels::t_arclabels_af*
af_compute(warthog::graph::planar_graph* g, std::vector<uint32_t>* part,
        warthog::arclabels::af_params par);

// arcflags computed on a graph that is a contraction hierarchy.
// In a CH edges are either "up" or "down". Nodes reached via an up edge
// have all of their successors generated and relaxed. Nodes reached via
// a down edge only generate and relax those successors further down 
// in the hierarchy.
warthog::arclabels::t_arclabels_af*
ch_af_compute(warthog::graph::planar_graph* g, std::vector<uint32_t>* part,
        std::vector<uint32_t>* rank, warthog::arclabels::af_params par);

// arcflags that are:
//  - computed on a corner point graph
//  - the graph is also a contraction hierarchy 
//  - the set of successors are pruned ala Jump Point Search
warthog::arclabels::t_arclabels_af*
ch_af_jpg_compute(warthog::graph::corner_point_graph* g, 
        std::vector<uint32_t>* part, std::vector<uint32_t>* rank, 
        warthog::arclabels::af_params);



} // warthog::arcflags::

} // warthog::

#endif

