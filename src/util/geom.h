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
        x1 = y1 = INT32_MAX;
        x2 = y2 = INT32_MIN;
    }

    rectangle(int32_t x1, int32_t y1, int32_t x2, int32_t y2)
    {
        this->x1 = x1; this->y1 = y1; 
        this->x2 = x2; this->y2 = y2;
    }

    int32_t
    get_width() { return x2 - x1; }

    int32_t
    get_height() { return y2 - y1; }

    void
    grow(int32_t x, int32_t y)
    {
        if(x < x1) { x1 = x; }
        if(y < y1) { y1 = y; }
        if(x > x2) { x2 = x; }
        if(y > y2) { y2 = y; }
    }

    bool 
    contains(int32_t x, int32_t y)
    {
        return x >= x1 && x <= x2 && y >= y1 && y <= y2;
    }

    bool 
    is_valid() { return x2 >= x1 && y2 >= y1; }

    void
    print(std::ostream& out)
    {
        out << " rect " << x1 << "\t" << y1 << "\t"
            << x2 << "\t" << y2 <<std::endl;
        }
};

}
}

#endif

