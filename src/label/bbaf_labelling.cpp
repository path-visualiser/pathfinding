#include "bbaf_labelling.h"

#include <fstream>
#include <pthread.h>
warthog::label::bbaf_labelling::bbaf_labelling(
        warthog::graph::planar_graph* g, 
        std::vector<uint32_t>* partitioning)
{
    g_ = g;
    part_ = partitioning;

    // figure out how many bytes are required per label
    uint32_t max_id = *(std::max_element(part_->begin(), part_->end()));
    bytes_per_af_label_ = (max_id / 8) + !!(max_id % 8);
    
    // allocate memory for each label
    labels_.resize(g_->get_num_nodes());
    for(uint32_t i = 0; i < g_->get_num_nodes(); i++)
    {
        warthog::graph::node* n = g_->get_node(i);
        for(uint32_t j = 0; j < n->out_degree(); j++)
        {
            // create one label per outgoing edge
            bbaf_label label;

            // initialise arcflag bits to zero
            label.flags_ = new uint8_t[bytes_per_af_label_];
            for(uint32_t k = 0; k < bytes_per_af_label_; k++)
            {
                label.flags_[k] = 0;
            }

            // add the label to the collection
            labels_.at(i).push_back(label);
        }
    }
}

warthog::label::bbaf_labelling::~bbaf_labelling()
{
    for(uint32_t i = 0; i < labels_.size(); i++)
    {
        std::vector<bbaf_label>& arclabs = labels_.at(i);
        for(uint32_t j = 0; j < arclabs.size(); j++)
        {
            bbaf_label& label = arclabs.back();
            delete [] label.flags_;
            arclabs.pop_back();
        }
    }
}

void
warthog::label::bbaf_labelling::print(std::ostream& out, 
        uint32_t first_id, uint32_t last_id)
{
    out 
    << "# Each line (aside from comments & header) is an edge label.\n"
    << "# A label comprises a set of arcflags (written out as 64bit words) \n"
    << "# and a rectangular bounding box which contains all nodes in the \n"
    << "# down-closure of the associated edge. The last 4 digits of each\n"
    << "# line describe the bounding box: minx, miny, maxx, maxy.\n"
    << "#\n"
    << "# NB: Before printing, labels are sorted by edge-tail index. \n"
    << "# Ties are broken using the order edges appear in the file \n"
    << "# " << g_->get_filename() << std::endl
    << "# (lower/earlier is better)\n";

    // we split labels into 64bit words for printing
    const uint32_t word_sz = sizeof(uint64_t);
    uint32_t words_per_label = bytes_per_af_label_ / word_sz;
    if((bytes_per_af_label_ % word_sz) > 0) words_per_label++;

    // iterate over the labels for each outgoing arc of each node
    warthog::geom::rectangle dummy;
    last_id = last_id > g_->get_num_nodes() ? 
        g_->get_num_nodes() : last_id;
    for(uint32_t i = first_id; i < last_id; i++)
    {
        std::vector<bbaf_label>& arclabs = labels_.at(i);
        for(uint32_t j = 0; j < arclabs.size(); j++)
        {
            // print the arcflags
            bbaf_label& label = arclabs.at(j);
            for(uint32_t word = 0; word < words_per_label; word++)
            {
                uint8_t printed_word[word_sz];
                for(uint32_t k = 0; k < word_sz; k++)
                {
                    // read the flag label, byte by byte, one word at a time
                    if((word*word_sz+k) < bytes_per_af_label_)
                    {
                        printed_word[k] = label.flags_[word*word_sz+k];
                    }
                    // pad the last word with leading zeroes if necessary
                    else
                    {
                        printed_word[k] = 0;
                    }
                }
                out << *((uint64_t*)&(printed_word)) << " ";
            }

            // print the bounding box
            assert(label.bbox_ == dummy || label.bbox_.is_valid());
            out << label.bbox_.x1 << " " << label.bbox_.y1 << " " 
                << label.bbox_.x2 << " " << label.bbox_.y2;
            out << std::endl;
        }
    }
}

warthog::label::bbaf_labelling*
warthog::label::bbaf_labelling::load(const char* filename, 
        warthog::graph::planar_graph* g, 
        std::vector<uint32_t>* partitioning)
{
    std::cerr << "loading bbaf file\n";
    std::ifstream ifs(filename);

    if(!ifs.good())
    {
        std::cerr << "\nerror trying to load bbaf labelling file " 
            << filename << std::endl;
        ifs.close();
        return 0;
    }

    warthog::label::bbaf_labelling* lab = 
        new warthog::label::bbaf_labelling(g, partitioning);

    uint32_t lines = 1;

    // skip comment lines
    while(ifs.peek() == '#')
    {
        while(ifs.get() != '\n') { lines++; }
    }

    // read labels for each outgoing arc
    std::string token;
    const uint32_t word_sz = sizeof(uint64_t);
    uint32_t words_per_label = (lab->bytes_per_af_label_ / word_sz);
    if(lab->bytes_per_af_label_ % word_sz != 0) { words_per_label++; }
    for(uint32_t i = 0; i < lab->g_->get_num_nodes(); i++)
    {
        warthog::graph::node* n = lab->g_->get_node(i);
        for(uint32_t j = 0; j < n->out_degree(); j++)
        {
            for(uint32_t word = 0; word < words_per_label; word++)
            {
                if(!ifs.good())
                {
                    std::cerr 
                        << "unexpected error reading bbaf label on line "
                        << lines << "; aborting\n";
                    delete lab;
                    return 0;
                }

                uint64_t label;
                ifs >> label;
                *((uint64_t*)
                    &lab->labels_.at(i).at(j).flags_[word*word_sz]) = label;
                assert(*((uint64_t*)
                            &lab->labels_.at(i).at(j).flags_[word*word_sz]) 
                        == label);
            }

            if(!ifs.good())
            {
                std::cerr 
                    << "unexpected error reading bbaf label on line "
                    << lines << "; aborting\n";
                delete lab;
                return 0;
            }

            uint32_t x1, y1, x2, y2;
            ifs >> x1 >> y1 >> x2 >> y2;
            if(x2 < x1 || y2 < y1)
            {
                std::cerr << "err; invalid label on line " << lines <<"\n";
                delete lab;
                return 0;
            }

            lab->labels_.at(i).at(j).bbox_.x1 = x1;
            lab->labels_.at(i).at(j).bbox_.y1 = y1;
            lab->labels_.at(i).at(j).bbox_.x2 = x2;
            lab->labels_.at(i).at(j).bbox_.y2 = y2;
            lines++;
       }
    }
    return lab;
}
