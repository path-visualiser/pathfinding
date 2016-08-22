#ifndef WARTHOG_HELPERS_H
#define WARTHOG_HELPERS_H

// helpers.h
//
// Helper functions that don't fit anywhere else.
//
// @author: dharabor
// @created: 21/08/2012
//

#include <fstream>
#include <iostream>
#include <vector>

namespace warthog
{
namespace helpers
{

// convert id into x/y coordinates on a grid of width 'mapwidth'
inline void
index_to_xy(unsigned int id, unsigned int mapwidth, 
		unsigned int& x, unsigned int& y)
{	
	y = id / mapwidth;
	x = id % mapwidth;
}

// sometimes the nodes of the search graph need to be labeled with
// integer data (e.g. during a partitioning of the graph)
// this file loads up such integer labels from a file. It assumes each
// line of the file contains a single integer.
// Comments lines can exist in the file; these need to begin with one
// of the following characters: #, %, c
inline void
load_integer_labels(const char* filename, std::vector<uint32_t>& labels)
{
    std::ifstream ifs(filename, std::ios_base::in);
    if(!ifs.good())
    {
        std::cerr << "\nerror trying to load file " << filename << std::endl;
        ifs.close();
        return;
    }

    while(true)
    {
        // skip comment lines
        while(ifs.peek() == '#' || ifs.peek() == 'c' || ifs.peek() == '%')
        {
            while(ifs.get() != '\n');
        }

        uint32_t tmp;
        ifs >> tmp;
        if(!ifs.good()) { break; }
        labels.push_back(tmp);
    }
    ifs.close();
}

inline void
load_integer_labels_dimacs(const char* filename, std::vector<uint32_t>& labels)
{
    // add a dummy to the front of the list if the labels are for use with
    // a dimacs graph. Such graphs use a 1-indexed scheme for node ids. we 
    // add the dummy so we can use the dimacs ids without conversion
    labels.push_back(0);
    load_integer_labels(filename, labels);
}


}
}

#endif


