#ifndef WARTHOG_PROBLEM_INSTANCE_H
#define WARTHOG_PROBLEM_INSTANCE_H

#include "search_node.h"

namespace warthog
{

class problem_instance
{
	public:
        problem_instance(uint32_t startid, uint32_t targetid, bool verbose=0) :
            start_id_(startid), 
            target_id_(targetid), 
            instance_id_(instance_counter_++),
            verbose_(verbose)

        { }
        
		problem_instance() :  
            start_id_(0), 
            target_id_(0), 
            instance_id_(instance_counter_++),
            verbose_(0)

        { }

		~problem_instance() { } 

		problem_instance(const warthog::problem_instance& other) 
        { 
            this->start_id_ = other.start_id_;
            this->target_id_ = other.target_id_;
            this->instance_id_ = other.instance_id_;
            this->verbose_ = other.verbose_;
        }

		warthog::problem_instance& 
		operator=(const warthog::problem_instance& other) 
        { 
            this->start_id_ = other.start_id_;
            this->target_id_ = other.target_id_;
            this->instance_id_ = other.instance_id_;
            this->verbose_ = other.verbose_;
            return *this; 
        }

        void
        print(std::ostream& out)
        {
            out << "problem instance; start_id = " << start_id_ << " "
                << " target_id " << target_id_ << " " << " search_id " 
                << instance_id_;
        }

		uint32_t start_id_;
		uint32_t target_id_;
		uint32_t instance_id_;
        bool verbose_;

        private:
            static uint32_t instance_counter_;

};


}

#endif

