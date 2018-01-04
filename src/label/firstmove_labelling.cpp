#include "firstmove_labelling.h"
#include "planar_graph.h"

#include <algorithm>

warthog::label::firstmove_labelling::firstmove_labelling(
        warthog::graph::planar_graph* g)
{
    lab_ = new std::vector< std::vector < warthog::label::fm_run > > ();
    lab_->resize(g_->get_num_nodes());
}

warthog::label::firstmove_labelling::~firstmove_labelling()
{ 
    delete lab_;
}

std::istream&
warthog::label::operator>>(std::istream& in, fm_run& the_run)
{
    in.read((char*)(&the_run.head_), 4);
    in.read((char*)(&the_run.label_), 1);
    return in ;
}

std::ostream&
warthog::label::operator<<(std::ostream& out, fm_run& the_run)
{
    out.write((char*)(&the_run.head_), 4);
    out.write((char*)(&the_run.label_), 1);
    return out;
}

std::istream&
warthog::label::operator>>(std::istream& in, 
        warthog::label::firstmove_labelling& lab)
{
    lab.lab_->resize(lab.g_->get_num_nodes());

    uint32_t num_rows;
    in.read((char*)(&num_rows), 4);

    for(uint32_t row = 0; row < num_rows; row++)
    {
        uint32_t row_id; 
        in.read((char*)(&row_id), 4);

        if(row_id > lab.g_->get_num_nodes())
        {
            std::cerr << "err; " << row_id 
                << " is aninvalid row identifier. aborting.\n";
            break;
        }

        if(lab.lab_->at(row_id).size() != 0)
        {
            std::cerr << "err; row id "<< row_id 
                << " appears more than once. aborting.\n";
            break;
        }

        uint32_t num_runs;
        in.read((char*)(&num_runs), 4);

        // read all the runs for the current row
        for(uint32_t i = 0; i < num_runs; i++)
        {
            fm_run tmp;
            in >> tmp;
            lab.lab_->at(row_id).push_back(tmp);

            assert(in.good());

            if(!in.good())
            {
                std::cerr << "err; while reading firstmove labels\n";
                std::cerr 
                    << "[debug info] row# " << row
                    << " row_id " << row_id 
                    << " run# " << i << " of " << num_runs 
                    << ". aborting.\n";
                return in;
            }
        }
    }
    return in;
}

std::ostream&
warthog::label::operator<<(std::ostream& out,
        warthog::label::firstmove_labelling& lab)
{
    uint32_t num_rows = 0;
    for(uint32_t n_id = 0; n_id < lab.g_->get_num_nodes(); n_id++)
    {
        if(lab.lab_->size() > 0) { num_rows++; }
    }
    out.write((char*)(&num_rows), 4);

    for(uint32_t row_id = 0; row_id < lab.g_->get_num_nodes(); row_id++)
    {
        if(lab.lab_->size() == 0) { continue; }
        out.write((char*)(&row_id), 4);

        uint32_t num_runs = lab.lab_->at(row_id).size();
        out.write((char*)(&num_runs), 4);

        for(uint32_t run = 0; run < num_runs; run++)
        {
            out << lab.lab_->at(row_id).at(run);
            if(!out.good())
            {
                std::cerr << "err; while writing labels\n";
                std::cerr 
                    << "[debug info] "
                    << " row_id " << row_id 
                    << " run# " << run 
                    << ". aborting.\n";
                return out;
            }
        }
    }
    return out;
}

