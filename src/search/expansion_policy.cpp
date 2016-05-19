#include "expansion_policy.h"

warthog::expansion_policy::expansion_policy(uint32_t num_nodes)
{
    nodepool_ = new warthog::blocklist(num_nodes);
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
