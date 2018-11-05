#include "node_pool.h"
#include "helpers.h"
#include "search_node.h"

warthog::mem::node_pool::node_pool(uint32_t num_nodes)
	: blocks_(0), pool_(0)
{
    init(num_nodes);
}

void
warthog::mem::node_pool::init(uint32_t num_nodes)
{
	num_blocks_ = ((num_nodes) >> warthog::mem::node_pool_ns::LOG2_NBS)+1;
	blocks_ = new warthog::search_node**[num_blocks_];
	for(uint32_t i=0; i < num_blocks_; i++)
	{
		blocks_[i] = 0;
	}
	blockspool_ = 
		new warthog::mem::cpool(sizeof(void*)*warthog::mem::node_pool_ns::NBS, 1);
	pool_ = new warthog::mem::cpool(sizeof(warthog::search_node));
}

warthog::mem::node_pool::~node_pool()
{
	delete pool_;
	delete blockspool_;
    delete blocks_;
}

warthog::search_node*
warthog::mem::node_pool::generate(uint32_t node_id)
{
	uint32_t block_id = node_id >> warthog::mem::node_pool_ns::LOG2_NBS;
	uint32_t list_id = node_id &  warthog::mem::node_pool_ns::NBS_MASK;
	assert(block_id <= num_blocks_);

	if(!blocks_[block_id])
	{
		// add a new block of nodes
		//std::cerr << "generating block: "<<block_id<<std::endl;
		warthog::search_node** list = new (blockspool_->allocate())
		   	warthog::search_node*[warthog::mem::node_pool_ns::NBS];
		uint32_t i = 0;
		for( ; i < warthog::mem::node_pool_ns::NBS; i+=8)
		{
			list[i] = 0;
			list[i+1] = 0;
			list[i+2] = 0;
			list[i+3] = 0;
			list[i+4] = 0;
			list[i+5] = 0;
			list[i+6] = 0;
			list[i+7] = 0;
		}
		// generate node_id
		warthog::search_node* mynode = new (pool_->allocate())
			warthog::search_node(node_id);
		list[list_id] = mynode;
		blocks_[block_id] = list;
		return mynode;
	}

	// look for node_id in an existing block
	warthog::search_node* mynode = blocks_[block_id][list_id];
	if(mynode)
	{
		return mynode;
	}

	// not in any existing block; generate it
	mynode = new (pool_->allocate()) warthog::search_node(node_id);
	blocks_[block_id][list_id] = mynode;
	return mynode;
}

warthog::search_node*
warthog::mem::node_pool::get_ptr(uint32_t node_id)
{
	uint32_t block_id = node_id >> warthog::mem::node_pool_ns::LOG2_NBS;
	uint32_t list_id = node_id &  warthog::mem::node_pool_ns::NBS_MASK;
	assert(block_id <= num_blocks_);

	if(!blocks_[block_id])
    {
        return 0;
    }
    return blocks_[block_id][list_id];
}

uint32_t
warthog::mem::node_pool::mem()
{
	uint32_t bytes = sizeof(*this) + blockspool_->mem() +
		pool_->mem() + num_blocks_*sizeof(warthog::search_node**);

	return bytes;
}
