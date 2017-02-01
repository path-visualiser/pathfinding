#include "arclabels.h"
#include "af_filter.h"
#include "constants.h"
#include "flexible_astar.h"
#include "fch_expansion_policy.h"
#include "graph_expansion_policy.h"
#include "planar_graph.h"
#include "problem_instance.h"
#include "zero_heuristic.h"

#include <fstream>
#include <iostream>
#include <set>

warthog::af_filter::af_filter(std::vector<uint32_t>* part)
{
    g_ = 0;
    t_byte_ = 0;
    t_bitmask_ = 0;

    par_ = warthog::arclabels::get_af_params(part);
    part_ = part;
}

warthog::af_filter::~af_filter()
{
    for(uint32_t i = 0; i < flags_.size(); i++)
    {
        std::vector<uint8_t*>& node_flags = flags_.at(i);
        for(uint32_t j = 0; j < node_flags.size(); j++)
        {
            uint8_t* label = node_flags.back();
            delete [] label;
            node_flags.pop_back();
        }
    }
}

void
warthog::af_filter::print(std::ostream& out)
{
    warthog::arclabels::af_print(flags_, par_, out);
}

bool
warthog::af_filter::load_labels(const char* filename, 
        warthog::graph::planar_graph* g)
{
    g_ = g;

    std::cerr << "loading arcflags file\n";
    std::ifstream ifs(filename);
    uint32_t firstid, lastid, num_parts;

    // skip comment lines
    while(ifs.peek() == '#')
    {
        while(ifs.get() != '\n');
    }

    // read header
    bool bad_header = false;
    std::string token;

    ifs >> token;
    if(token == "firstnode")
    {
        ifs >> firstid;

        ifs >> token;
        if(token == "lastnode")
        {
            ifs >> lastid;
        }
        else { bad_header = true; }

        ifs >> token;
        if(token == "partitions")
        {
            ifs >> num_parts;
            if(num_parts != par_.nparts_)
            { 
                std::cerr 
                    << "err; number of partitions in arcflags file\n"
                    << " does not match input graph (read: " << num_parts 
                    << " expected: " << par_.nparts_ << ")\n";
                return false;
            }
        }
        else { bad_header = true; }
    }
    else { bad_header = true; }

    if(bad_header)
    {
        std::cerr << "err; arcflags file has invalid header. required format: \n"
            << "firstnode [integer] lastnode [integer] partitions [integer]\n";
        return false;
    }

    // allocate enough memory to label the outgoing edge of every node 
    // and then read the labels from file
    // NB: label range is : [firstid, lastid)
    if(lastid == g->get_num_nodes()-1)
    {
        lastid++;
    }
    flags_.resize(lastid);
    const uint32_t word_sz = sizeof(uint64_t);
    uint32_t words_per_label = ceil(par_.bytes_per_label_ / (double)word_sz);

    for(uint32_t i = firstid; i < lastid; i++)
    {
        warthog::graph::node* n = g_->get_node(i);
        flags_.at(i).resize(n->out_degree());

        for(uint32_t j = 0; j < n->out_degree(); j++)
        {
            uint8_t* label = new uint8_t[par_.bytes_per_label_];
            flags_.at(i).at(j) = label;

            for(uint32_t word = 0; word < words_per_label; word++)
            {
                if(!ifs.good())
                {
                    std::cerr << "errr; unexpected read error loading arcflags; "
                        << "aborting\n";
                    return false;
                }

                uint64_t label = 0;
                ifs >> label;
                //std::cout << "node " << i << " arc " << j << 
                //" label " << label << std::endl;
                *((uint64_t*)&flags_.at(i).at(j)[word*word_sz]) = label;
                assert(*((uint64_t*)&flags_.at(i).at(j)[word*word_sz])
                            == label);
            }
        }
    }
    return true;
}

