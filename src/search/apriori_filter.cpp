#include "apriori_filter.h"
#include "search_node.h"
#include <cassert>

warthog::apriori_filter::apriori_filter(uint32_t num_elements) 
{
    filter_sz_ = (num_elements >> warthog::LOG2_DBWORD_BITS)+1;
    filter_ = new warthog::dbword[filter_sz_];
    reset_filter();
}

warthog::apriori_filter::~apriori_filter()
{
    delete [] filter_;
}

void 
warthog::apriori_filter::set_flag_true(uint32_t node_id)
{
    int index = node_id >> warthog::LOG2_DBWORD_BITS;
    int pos = node_id & DBWORD_BITS_MASK;
    filter_[index] |= (1 << pos);
}

void
warthog::apriori_filter::set_flag_false(uint32_t node_id)
{
    int index = node_id >> warthog::LOG2_DBWORD_BITS;
    int pos = node_id & DBWORD_BITS_MASK;
    filter_[index] &= ~(1 << pos);
}

void 
warthog::apriori_filter::reset_filter()
{
    for(uint32_t i = 0; i < filter_sz_; i++)
    {
        filter_[i] = 0;
    }
}

bool
warthog::apriori_filter::filter(uint32_t node_id, uint32_t edge_idx)
{
    return get_flag(node_id);
}

bool
warthog::apriori_filter::get_flag(uint32_t id) 
{
   assert((id / warthog::DBWORD_BITS) < filter_sz_);
   uint32_t word = id / warthog::DBWORD_BITS;
   uint32_t pos = id % warthog::DBWORD_BITS;
   return this->filter_[word] & (1 << pos);
}

