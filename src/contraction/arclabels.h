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

}

}

#endif

