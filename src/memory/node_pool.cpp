#include "node_pool.h"
#include "helpers.h"
#include "search_node.h"

warthog::mem::node_pool::node_pool(uint32_t num_nodes)
	: blocks_(0)
{
    init(num_nodes);
}

void
warthog::mem::node_pool::init(uint32_t num_nodes)
{
	num_blocks_ = ((num_nodes) >> warthog::mem::node_pool_ns::LOG2_NBS)+1;
	blocks_ = new warthog::search_node*[num_blocks_];
	for(uint32_t i=0; i < num_blocks_; i++)
	{
		blocks_[i] = 0;
	}
    uint32_t block_sz = 
        warthog::mem::node_pool_ns::NBS * sizeof(warthog::search_node);
	//uint32_t mem_chunks_prealloc = 
    //    (((num_blocks_ / 10) * block_sz) / warthog::mem::DEFAULT_CHUNK_SIZE)+1;
   // blockspool_ = new warthog::mem::cpool(block_sz, mem_chunks_prealloc);
    blockspool_ = new warthog::mem::cpool(block_sz, 1);
    
    
    // we keep one bit per node to track whether memory is initialised
    //node_init_sz_ = (num_nodes >> 6) + 1;
    //node_init_ = new uint64_t[node_init_sz_];
    //for(uint32_t i = 0; i < node_init_sz_; i++)
    //{
    //    node_init_[i] = 0;
    //}
}

warthog::mem::node_pool::~node_pool()
{
    //delete [] node_init_;

	blockspool_->reclaim();
	delete blockspool_;

	for(uint32_t i=0; i < num_blocks_; i++)
	{
		if(blocks_[i] != 0)
		{
			//std::cerr << "deleting block: "<<i<<std::endl;
			blocks_[i] = 0;
		}
	}
    delete blocks_;
}

warthog::search_node*
warthog::mem::node_pool::generate(uint32_t node_id)
{
	uint32_t block_id = node_id >> warthog::mem::node_pool_ns::LOG2_NBS;
	uint32_t list_id = node_id &  warthog::mem::node_pool_ns::NBS_MASK;
	assert(block_id <= num_blocks_);

    // add a new block of nodes if necessary
	if(!blocks_[block_id])
	{
		//std::cerr << "generating block: "<<block_id<<std::endl;
        blocks_[block_id] = new (blockspool_->allocate())
		   	warthog::search_node[warthog::mem::node_pool_ns::NBS];

        // initialise memory 
        uint32_t current_id = node_id - list_id;
		for( uint32_t i  = 0; i < warthog::mem::node_pool_ns::NBS; i+=8)
		{
            new (&blocks_[block_id][i]) warthog::search_node(current_id++);
            new (&blocks_[block_id][i+1]) warthog::search_node(current_id++);
            new (&blocks_[block_id][i+2]) warthog::search_node(current_id++);
            new (&blocks_[block_id][i+3]) warthog::search_node(current_id++);
            new (&blocks_[block_id][i+4]) warthog::search_node(current_id++);
            new (&blocks_[block_id][i+5]) warthog::search_node(current_id++);
            new (&blocks_[block_id][i+6]) warthog::search_node(current_id++);
            new (&blocks_[block_id][i+7]) warthog::search_node(current_id++);
		}
	}

    // initialise memory 
//    uint64_t flag = 1;
//    if(!(node_init_[node_id >> 6] & (flag << (node_id & 63))))
//    {
//        new (&blocks_[block_id][list_id]) warthog::search_node(node_id);
//        node_init_[node_id >> 6] |= (flag << (node_id & 63));
//    }

//    warthog::search_node* blah = &(blocks_[block_id][list_id]);
//    assert(blah->get_id() == node_id);

	// return the node from its position in the assocated block 
    return &(blocks_[block_id][list_id]);
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
    return &(blocks_[block_id][list_id]);
}

uint32_t
warthog::mem::node_pool::mem()
{
	uint32_t bytes = 
        sizeof(*this) + 
        blockspool_->mem() +
		num_blocks_*sizeof(void*);

	return bytes;
}
