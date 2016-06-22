#include "node_filter.h"
#include <cassert>

warthog::node_filter::node_filter(uint32_t num_elements)
{
    filter_sz_ = (num_elements >> warthog::LOG2_DBWORD_BITS)+1;
    filter_ = new warthog::dbword[filter_sz_];
}

warthog::node_filter::~node_filter()
{
    delete [] filter_;
}

void 
warthog::node_filter::filter(uint32_t node_id)
{
    int index = node_id >> warthog::LOG2_DBWORD_BITS;
    int pos = node_id & DBWORD_BITS_MASK;
    filter_[index] |= 1 << pos;
}

void
warthog::node_filter::unfilter(uint32_t node_id)
{
    int index = node_id >> warthog::LOG2_DBWORD_BITS;
    int pos = node_id & DBWORD_BITS_MASK;
    filter_[index] &= ~(1 << pos);
}

void 
warthog::node_filter::reset_filter()
{
    for(uint32_t i = 0; i < filter_sz_; i++)
    {
        filter_[i] = 0;
    }
}

bool
warthog::node_filter::get_filter_flag(uint32_t id) 
{
   assert((id / warthog::DBWORD_BITS) < filter_sz_);
   uint32_t word = id / warthog::DBWORD_BITS;
   uint32_t pos = id % warthog::DBWORD_BITS;
   return this->filter_[word] & (1 << pos);
}

