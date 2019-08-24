#ifndef WARTHOG_CH_CH_DATA_H
#define WARTHOG_CH_CH_DATA_H

#include "xy_graph.h"

namespace warthog
{
namespace ch
{

class ch_data
{

public:
    warthog::graph::xy_graph* g_;
    std::vector<uint32_t>* level_;

};

// load a contraction hierarchy for use with bi-directional search (BCH).
// In BCH, only outgoing up edges and incoming up edges are stored with
// the graph. Other edges, which may be read from file, are discarded. 
bool
load_bch_data(warthog::ch::ch_data*, const char* filename);

bool
save_ch_data(const char* filename, warthog::ch::ch_data* chd);

}
}

#endif
