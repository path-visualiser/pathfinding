#ifndef WARTHOG_CH_CH_DATA_H
#define WARTHOG_CH_CH_DATA_H

#include "xy_graph.h"

namespace warthog
{
namespace ch
{

typedef enum
{
    UP_DOWN = 0,
    UP_ONLY = 1
} ch_direction_t;

class ch_data
{

public:
    ch_data(uint32_t num_nodes, warthog::ch::ch_direction_t type);
    virtual ~ch_data();

    warthog::graph::xy_graph* g_;
    std::vector<uint32_t>* level_;
    ch_direction_t type_;
};

warthog::ch::ch_data* 
load_ch_data(const char* filename, warthog::ch::ch_direction_t type);

bool
save_ch_data(const char* filename, warthog::ch::ch_data* chd);

}
}

#endif
