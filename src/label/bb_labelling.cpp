#include "bb_labelling.h"

#include <algorithm>
#include <fstream>

warthog::label::bb_labelling::bb_labelling(warthog::graph::planar_graph* g)
{
    g_ = g;

    // allocate memory for every label
    labels_ = new std::vector<std::vector<warthog::geom::rectangle>>();
    labels_->resize(g->get_num_nodes());
    for(uint32_t i = 0; i < g->get_num_nodes(); i++)
    {
        warthog::graph::node* n = g->get_node(i);
        labels_->at(i).resize(n->out_degree());
    }
}

warthog::label::bb_labelling::~bb_labelling()
{
    delete labels_;
}

void
warthog::label::bb_labelling::print(
        std::ostream& out, uint32_t first_id, uint32_t last_id)
{
    out << "# Before printing labels are sorted by tail index. \n"
        << "# Ties are broken using the order edges appear in the file \n"
        << "# " << g_->get_filename() << std::endl
        << "# (lower is better)\n";

    warthog::geom::rectangle dummy;
    for(uint32_t i = first_id; i < last_id; i++)
    {
        std::vector<warthog::geom::rectangle> vrect = labels_->at(i);
        for(uint32_t j = 0; j < vrect.size(); j++)
        {
            warthog::geom::rectangle rect = vrect.at(j);
            if(rect != dummy) { assert(rect.is_valid()); }
            out 
                //<< i << " " << j << " " << 
                << rect.x1 << "\t" << rect.y1 << "\t"
                << rect.x2 << "\t" << rect.y2 << "\t"
                << std::endl;
        }
    }
}

warthog::label::bb_labelling*
warthog::label::bb_labelling::load(
        const char* filename, warthog::graph::planar_graph* g)
{
    std::ifstream ifs(filename, std::ios_base::in);
    if(!ifs.good())
    {
        std::cerr << "\nerror trying to load bbox values from file " 
            << filename << std::endl;
        ifs.close();
        return 0;
    }

    uint32_t lines = 0;
    // skip comment lines
    while(ifs.peek() == '#')
    {
        while(ifs.good() && ifs.get() != '\n');
        lines++;
    }

    warthog::label::bb_labelling* bbl = new warthog::label::bb_labelling(g);
    warthog::geom::rectangle dummy;
    for(uint32_t i = 0; i < g->get_num_nodes(); i++)
    {
        warthog::graph::node* n = g->get_node(i);
        for(uint32_t j = 0; j < n->out_degree(); j++)
        {
            int32_t x1, y1, x2, y2;
            ifs >> x1 >> y1 >> x2 >> y2;
            if(ifs.eof() || !ifs.good())
            {
                std::cerr << "unexpected error while reading" 
                    << filename << "; line=" << lines << "\n";
                delete bbl; 
                return 0;
            }

            warthog::geom::rectangle rect = 
                warthog::geom::rectangle(x1, y1, x2, y2);
            if(rect != dummy && !rect.is_valid())
            {
                std::cerr << "err; invalid label on line " << lines <<"\n";
                delete bbl;
                return 0;
            }
            bbl->labels_->at(i).at(j) = rect;
            lines++;
        }
    }
    ifs.close();
    return bbl;
}
