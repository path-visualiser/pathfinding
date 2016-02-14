#include "bidirectional_search.h"
#include "constants.h"
#include "expansion_policy.h"
#include "pqueue.h"
#include "search_node.h"

warthog::bidirectional_search::
    bidirectional_search(
            warthog::expansion_policy* expander, 
            warthog::heuristicFn heuristic)
    : expander_(expander), heuristic_(heuristic)
{
    verbose_ = false;
    fwd_open_ = new pqueue(512, true);
    bwd_open_ = new pqueue(512, true);
}

warthog::bidirectional_search::~bidirectional_search()
{
    delete fwd_open_;
    delete bwd_open_;
}


std::stack<uint32_t> 
warthog::bidirectional_search::get_path(uint32_t startid, uint32_t goalid)
{
    std::stack<uint32_t> path;
//    warthog::search_node* goal = search(startid, goalid);
//    if(goal)
//    {
//        // follow backpointers to extract the path
//        assert(goal->get_id() == goalid);
//        for(warthog::search_node* cur = goal;
//                cur != 0;
//                cur = cur->get_parent())
//        {
//            path.push(cur->get_id());
//        }
//        assert(path.top() == startid);
//    }
//    cleanup();
//    return path;
    return path;
}

double 
warthog::bidirectional_search::get_length(uint32_t startid, uint32_t goalid)
{
    warthog::search_node* goal = search(startid, goalid);
    warthog::cost_t len = warthog::INF;
    if(goal)
    {
        assert(goal->get_id() == goalid);
        len = goal->get_g();
    }

#ifndef NDEBUG
    if(verbose_)
    {
        std::stack<warthog::search_node*> path;
        warthog::search_node* current = goal;
        while(current != 0)	
        {
            path.push(current);
            current = current->get_parent();
        }

        while(!path.empty())
        {
            warthog::search_node* n = path.top();
            int32_t x, y;
            expander_->get_xy(n, x, y);
            std::cerr 
                << "final path: (" << x << ", " << y << ")...";
            n->print(std::cerr);
            std::cerr << std::endl;
            path.pop();
        }
    }
#endif
    cleanup();
    return len / (double)warthog::ONE;
        
}

warthog::search_node*
warthog::bidirectional_search::search(uint32_t startid, uint32_t goalid)
{
    return 0;
}

size_t
warthog::bidirectional_search::mem()
{
    return sizeof(*this) + 
        fwd_open_->mem() +
        bwd_open_->mem() +
        expander_->mem();
}

void
warthog::bidirectional_search::cleanup()
{
    fwd_open_->clear();
    bwd_open_->clear();
    expander_->clear();
}
