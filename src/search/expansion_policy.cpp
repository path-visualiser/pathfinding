#include "expansion_policy.h"

warthog::expansion_policy::expansion_policy(uint32_t nodes_pool_size)
{
    nodes_pool_size_ = nodes_pool_size;
    nodepool_ = new warthog::blocklist(nodes_pool_size);
    neis_ = new std::vector<neighbour_record>();
    neis_->reserve(32);
}

warthog::expansion_policy::~expansion_policy()
{
    neis_->clear();
    delete neis_;
    delete nodepool_;
}
