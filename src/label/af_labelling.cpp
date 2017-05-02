#include "af_labelling.h"
#include <fstream>

warthog::label::af_labelling::af_labelling(
        warthog::graph::planar_graph* g, 
        std::vector<uint32_t>* partitioning)
{
    g_ = g;
    part_ = partitioning;

    // figure out how many bytes are required per label
    uint32_t max_id = *(std::max_element(part_->begin(), part_->end()));
    bytes_per_label_ = (max_id / 8) + !!(max_id % 8);

    // allocate memory
    flags_ = new std::vector<std::vector<uint8_t*>>();
    flags_->resize(g->get_num_nodes());
    for(uint32_t i = 0; i < g->get_num_nodes(); i++)
    {
        warthog::graph::node* source = g->get_node(i);
        for(uint32_t j = 0; j < source->out_degree(); j++)
        {
            uint8_t* label = new uint8_t[bytes_per_label_];
            flags_->at(i).push_back(label);
            // zero every label
            for(uint32_t k = 0; k < bytes_per_label_; k++)
            {
                label[k] = 0;
            }
        }
    }
}

warthog::label::af_labelling::~af_labelling()
{
    for(uint32_t i = 0; i < flags_->size(); i++)
    {
        std::vector<uint8_t*>& node_flags = flags_->at(i);
        for(uint32_t j = 0; j < node_flags.size(); j++)
        {
            uint8_t* label = node_flags.back();
            delete [] label;
            node_flags.pop_back();
        }
    }
    delete flags_;
}

void
warthog::label::af_labelling::print(
        std::ostream& out, uint32_t first_id, uint32_t last_id)
{
    if(flags_->size() == 0) { return; }

    out << "# arcflag labels for the graph\n";
    out << "# " << g_->get_filename() << std::endl;
    out << "# labels are written out using 64-bit words\n";
    out << "# with spaces separating words\n";
    // we split flags into 64bit words for printing
    const uint32_t word_sz = sizeof(uint64_t);
    uint32_t words_per_label =  bytes_per_label_ / word_sz;
    if((bytes_per_label_ % word_sz) > 0) words_per_label++;

    // iterate over the flags for each outgoing arc of each node
    last_id = last_id > g_->get_num_nodes() ? 
        g_->get_num_nodes() : last_id;
    for(uint32_t i = first_id; i < last_id; i++)
    {
        std::vector<uint8_t*>& node_arcflags_ = flags_->at(i);
        for(uint32_t j = 0; j < node_arcflags_.size(); j++)
        {
            uint8_t* label = node_arcflags_.at(j);
            for(uint32_t word = 0; word < words_per_label; word++)
            {
                uint8_t printed_word[word_sz];
                for(uint32_t k = 0; k < word_sz; k++)
                {
                    // read the label, byte by byte, one word at a time
                    if((word*word_sz+k) < bytes_per_label_) 
                    {
                        printed_word[k] = label[word*word_sz+k];
                    }
                    // pad the last word with leading zeroes if necessary
                    else
                    {
                        printed_word[k] = 0;
                    }
                }
                out << *((uint64_t*)&(printed_word))
                    << (word < words_per_label ? " " : "");
            }
            out << std::endl;
        }
    }
}

warthog::label::af_labelling*
warthog::label::af_labelling::load(
        const char* filename, warthog::graph::planar_graph* g,
        std::vector<uint32_t>* partitioning)
{
    std::cerr << "loading arcflags labelling file\n";
    std::ifstream ifs(filename, std::ios_base::in);
    if(!ifs.good())
    {
        std::cerr << "\nerror trying to load arcflags file " 
            << filename << std::endl;
        ifs.close();
        return 0;
    }

    // skip comment lines
    uint32_t lines = 1;
    while(ifs.peek() == '#')
    {
        while(ifs.get() != '\n');
        lines++;
    }

    warthog::label::af_labelling* afl = 
        new warthog::label::af_labelling(g, partitioning);

    // allocate enough memory to label the outgoing edge of every node 
    // and then read the labels from file
    // NB: label range is : [firstid, lastid)
    const uint32_t word_sz = sizeof(uint64_t);
    uint32_t words_per_label = ceil(afl->bytes_per_label_ / (double)word_sz);
    for(uint32_t i = 0; i < afl->g_->get_num_nodes(); i++)
    {
        warthog::graph::node* n = afl->g_->get_node(i);
        for(uint32_t j = 0; j < n->out_degree(); j++)
        {
            for(uint32_t word = 0; word < words_per_label; word++)
            {
                uint64_t label = 0;
                ifs >> label;
                if(ifs.eof() || !ifs.good())
                {
                    std::cerr << "unexpected error while reading " 
                        << filename << "; line=" << lines << "\n";
                    std::cerr 
                        << "[debug info] node: " << i 
                        << " out-edge-index: " << j << "\n";
                    delete afl;
                    return 0;
                }

                //std::cout << "node " << i << " arc " << j << 
                //" label " << label << std::endl;
                *((uint64_t*)&afl->flags_->at(i).at(j)[word*word_sz]) = label;
                assert(*((uint64_t*)&afl->flags_->at(i).at(j)[word*word_sz])
                            == label);
            }
            lines++;
        }
    }
    return afl;
}
