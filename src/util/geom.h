#ifndef WARTHOG_GEOM_H
#define WARTHOG_GEOM_H

// util/geom.h
//
// miscellaneous collection of geometric things
//
// @author: dharabor
// @created 2016-07-31
//

#include <ostream>
#include <float.h>

namespace warthog
{
namespace geom
{

struct rectangle
{
    int32_t x1, y1, x2, y2;

    rectangle()
    { 
        clear();
    }

    rectangle(int32_t x1, int32_t y1, int32_t x2, int32_t y2)
    {
        this->x1 = x1; this->y1 = y1; 
        this->x2 = x2; this->y2 = y2;
    }

    rectangle(const rectangle& other)
    {
        x1 = other.x1;
        x2 = other.x2;
        y1 = other.y1;
        y2 = other.y2;
    }

    bool
    operator==(const rectangle& other)
    {
        return 
            x1 == other.x1 && x2 == other.x2 &&
            y1 == other.y1 && y2 == other.y2;
    }

    bool
    operator!=(const rectangle& other)
    {
        return !this->operator==(other);
    }

    uint32_t
    get_width() { return x2 - x1; }

    uint32_t
    get_height() { return y2 - y1; }

    uint32_t
    get_area() { return (x2 - x1) * (y2 - y1); }

    void
    grow(int32_t x, int32_t y)
    {
        if(x < x1) { x1 = x; }
        if(y < y1) { y1 = y; }
        if(x > x2) { x2 = x; }
        if(y > y2) { y2 = y; }
    }

    void
    grow(warthog::geom::rectangle& r)
    {
        if(r.x1 < x1) { x1 = r.x1; }
        if(r.y1 < y1) { y1 = r.y1; }
        if(r.x2 > x2) { x2 = r.x2; } 
        if(r.y2 > y2) { y2 = r.y2; } 
    }

    bool 
    contains(int32_t x, int32_t y)
    {
        return x >= x1 && x <= x2 && y >= y1 && y <= y2;
    }

    bool 
    intersects(warthog::geom::rectangle& r)
    {
        rectangle tmp(*this);
        tmp.grow(r);
        uint32_t bound = 
            (this->get_width() + r.get_width()) *
            (this->get_height() + r.get_height());
        if(tmp.get_area() > bound) { return false; }
        return true;
    }

    bool 
    is_valid() { return x2 >= x1 && y2 >= y1; }

    void
    print(std::ostream& out)
    {
        out << " rect " << x1 << "\t" << y1 << "\t"
            << x2 << "\t" << y2 <<std::endl;
    }

    void
    clear()
    {
        x1 = y1 = INT32_MAX;
        x2 = y2 = INT32_MIN;
    }
};

}
}

#endif

