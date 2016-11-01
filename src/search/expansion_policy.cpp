#include "expansion_policy.h"

warthog::expansion_policy::expansion_policy(uint32_t nodes_pool_size)
{
    nodes_pool_size_ = nodes_pool_size;
    nodepool_ = new warthog::blocklist(nodes_pool_size);
    costs_ = new std::vector<double>();
    neis_ = new std::vector<warthog::search_node*>();
    costs_->reserve(8);
    neis_->reserve(8);
}

warthog::expansion_policy::~expansion_policy()
{
    neis_->clear();
    costs_->clear();
    delete neis_;
    delete costs_;
    delete nodepool_;
}
