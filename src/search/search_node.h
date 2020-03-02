#ifndef WARTHOG_SEARCH_NODE_H
#define WARTHOG_SEARCH_NODE_H

// search_node.h
//
// @author: dharabor
// @created: 10/08/2012
//

#include "constants.h"
#include "cpool.h"
#include "jps.h"

#include <iostream>

namespace warthog
{
	
class search_node
{
	public:
		search_node(warthog::sn_id_t id = warthog::SN_ID_MAX) :
            id_(id), parent_id_(warthog::SN_ID_MAX),
            g_(warthog::COST_MAX), f_(warthog::COST_MAX),
            status_(0), priority_(warthog::INF32), search_number_(0),
            ub_(warthog::COST_MAX)
		{
            set_pdir(warthog::jps::direction::NONE);
			refcount_++;
		}

		~search_node()
		{
			refcount_--;
		}

		inline void
		init(uint32_t search_number,
             warthog::sn_id_t parent_id,
             double g,
             double f,
             warthog::cost_t ub=warthog::COST_MAX)
		{
			status_ = 0;
            parent_id_= parent_id;
            f_ = f;
            g_ = g;
            ub_ = ub;
			search_number_ = search_number;
		}

		inline uint32_t
		get_search_number()
		{
			return search_number_;
		}

		inline void
		set_search_number(uint32_t search_number)
		{
			search_number_ = search_number;
		}

		inline warthog::sn_id_t
		get_id() const { return id_; }

		inline void
		set_id(warthog::sn_id_t id) 
		{ 
			id_ = id;
		} 

		inline warthog::jps::direction
		get_pdir() const
		{
			return (warthog::jps::direction)jps_parent_direction_;
			   	
		}

		inline void
		set_pdir(warthog::jps::direction d)
		{
			jps_parent_direction_ = d;
		}

		inline bool
		get_expanded() const { return status_; }

		inline void
		set_expanded(bool expanded) 
		{ 
            status_ = expanded;
		}

		inline warthog::sn_id_t
		get_parent() const { return parent_id_; }

		inline void
		set_parent(warthog::sn_id_t parent_id) { parent_id_ = parent_id; } 

		inline uint32_t
		get_priority() const { return priority_; }

		inline void
		set_priority(uint32_t priority) { priority_ = priority; } 

		inline double
		get_g() const { return g_; }

		inline void
		set_g(double g) { g_ = g; }

		inline double 
		get_f() const { return f_; }

		inline void
		set_f(double f) { f_ = f; }

		inline warthog::cost_t
		get_ub() const { return ub_; }

		inline void
		set_ub(warthog::cost_t ub) { ub_ = ub; }

		inline void
		relax(double g, warthog::sn_id_t parent_id)
		{
			assert(g < g_);
			f_ = (f_ - g_) + g;
			g_ = g;
			parent_id_ = parent_id;
		}

		inline bool
		operator<(const warthog::search_node& other) const
		{
			if(f_ < other.f_)
			{
				return true;
			}
			if(f_ > other.f_)
			{
				return false;
			}

			// break ties in favour of larger g
			if(g_ > other.g_)
			{
				return true;
			}
			return false;
		}

		inline bool
		operator>(const warthog::search_node& other) const
		{
			if(f_ > other.f_)
			{
				return true;
			}
			if(f_ < other.f_)
			{
				return false;
			}

			// break ties in favour of larger g
			if(g_ > other.g_)
			{
				return true;
			}
			return false;
		}

		inline bool
		operator==(const warthog::search_node& other) const
		{
			if( !(*this < other) && !(*this > other))
			{
				return true;
			}
			return false;
		}

		inline bool
		operator<=(const warthog::search_node& other) const
		{
			if(*this < other)
			{
				return true;
			}
			if(!(*this > other))
			{
				return true;
			}
			return false;
		}

		inline bool
		operator>=(const warthog::search_node& other) const
		{
			if(*this > other)
			{
				return true;
			}
			if(!(*this < other))
			{
				return true;
			}
			return false;
		}

		inline void
		print(std::ostream&  out) const
		{
			out << "search_node id:" << get_id();
            out << " p_id: ";
            out << parent_id_;
            out << " g: "<<g_ <<" f: "<<this->get_f() << " ub: " << ub_
                << " expanded: " << get_expanded() << " "
                << " search_number_: " << search_number_
                << " pdir: "<< get_pdir() << " ";
		}

		uint32_t
		mem()
		{
			return sizeof(*this);
		}

        static uint32_t
        get_refcount() { return refcount_; }

	private:
		warthog::sn_id_t id_; // bit 0 is expansion status; 1-63 are id
        warthog::sn_id_t parent_id_;
		double g_;
		double f_;
        uint8_t status_; // open or closed
		uint32_t priority_; // expansion priority

		uint32_t search_number_;
        uint8_t jps_parent_direction_; // hack
        warthog::cost_t ub_;


static uint32_t refcount_;
};

struct cmp_less_search_node
{
    inline bool
    operator()(
            const warthog::search_node& first,
            const warthog::search_node& second)
    {
        return first < second;
    }
};

struct cmp_greater_search_node
{
    inline bool
    operator()(
            const warthog::search_node& first,
            const warthog::search_node& second)
    {
        return first > second;
    }
};

struct cmp_less_search_node_f_only
{
    inline bool
    operator()(
            const warthog::search_node& first,
            const warthog::search_node& second)
    {
        return first.get_f() < second.get_f();
    }
};

}

std::ostream& operator<<(std::ostream& str, const warthog::search_node& sn);

#endif
