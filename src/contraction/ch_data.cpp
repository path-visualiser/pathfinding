#include "ch_data.h"

warthog::ch::ch_data::ch_data()
{
    g_ = new warthog::graph::xy_graph();
    level_ = new std::vector<uint32_t>(0);
}

warthog::ch::ch_data::~ch_data()
{
    delete g_;
    delete level_;
}

bool
warthog::ch::load_bch_data(warthog::ch::ch_data*, const char* filename)
{
    return false;
}

bool
warthog::ch::save_ch_data(const char* filename, warthog::ch::ch_data* chd)
{
    return false;
}

