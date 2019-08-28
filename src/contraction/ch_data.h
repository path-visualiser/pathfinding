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
    ch_data();
    virtual ~ch_data();

    warthog::graph::xy_graph* g_;
    std::vector<uint32_t>* level_;
};

warthog::ch::ch_data* 
load_ch_data(const char* filename, bool up_only = false);

bool
save_ch_data(const char* filename, warthog::ch::ch_data* chd);

}
}

#endif
