#ifndef WARTHOG_PROBLEM_INSTANCE_H
#define WARTHOG_PROBLEM_INSTANCE_H

#include "search_node.h"

namespace warthog
{

class problem_instance
{
	public:
        problem_instance(uint32_t startid, uint32_t targetid) :
            start_id_(startid), target_id_(targetid), search_id_(0) { }
        
		problem_instance() :  start_id_(0), target_id_(0), search_id_(0) { }

		~problem_instance() { } 

		problem_instance(const warthog::problem_instance& other) 
        { 
            this->start_id_ = other.start_id_;
            this->target_id_ = other.target_id_;
            this->search_id_ = other.search_id_;
        }

		warthog::problem_instance& 
		operator=(const warthog::problem_instance& other) 
        { 
            this->start_id_ = other.start_id_;
            this->target_id_ = other.target_id_;
            this->search_id_ = other.search_id_;
            return *this; 
        }

		inline void
		set_target_id(uint32_t target_id_id) { target_id_ = target_id_id; }

		inline uint32_t
		get_target_id() { return target_id_; }

		inline uint32_t
		get_start_id() { return start_id_; }

		inline void
		set_start_id(uint32_t start_id_id) { start_id_ = start_id_id; }

		inline uint32_t
		get_search_id() { return search_id_; } 

		inline void
		set_search_id(uint32_t search_id) { search_id_ = search_id; }

        void
        print(std::ostream& out)
        {
            out << "problem instance; start_id = " << start_id_ << " "
                << " target_id " << target_id_ << " " << " search_id " 
                << search_id_;
        }

	private:
		uint32_t start_id_;
		uint32_t target_id_;
		uint32_t search_id_;

};

}

#endif

