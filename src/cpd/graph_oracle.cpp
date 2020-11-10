#include "cpd.h"
#include "graph_oracle.h"
#include "helpers.h"

void
warthog::cpd::compute_row(uint32_t source_id,
                          warthog::cpd::graph_oracle* cpd,
                          warthog::search* dijk,
                          std::vector<warthog::cpd::fm_coll> &s_row)
{
    warthog::problem_instance problem(source_id);
    warthog::solution sol;

    std::fill(s_row.begin(), s_row.end(), warthog::cpd::CPD_FM_NONE);
    dijk->get_path(problem, sol);
    cpd->add_row(source_id, s_row);
}

void
warthog::cpd::graph_oracle::add_row(uint32_t source_id,
                 std::vector<warthog::cpd::fm_coll>& row)
{
    // source gets a wildcard move
    row.at(source_id) = warthog::cpd::CPD_FM_NONE;

    // greedily compress the row w.r.t. the current column order
    warthog::cpd::fm_coll moveset = row.at(order_.at(0));
    uint32_t head = 0;
    for(uint32_t index = 0; index < row.size(); index++)
    {
        assert(moveset > 0);
        if((moveset & row.at(order_.at(index))) == 0)
        {
            uint32_t firstmove = __builtin_ffsl(moveset) - 1;
            assert(firstmove < warthog::cpd::CPD_FM_MAX);
            fm_.at(source_id).push_back(
                    warthog::cpd::rle_run32{ (head << 4) | firstmove} );
            moveset = row.at(order_.at(index));
            head = index;
        }
        moveset = moveset & row.at(order_.at(index));
    }

    // add the last run
    uint32_t firstmove = __builtin_ffsl(moveset) - 1;
    assert(firstmove < warthog::cpd::CPD_FM_MAX);
    fm_.at(source_id).push_back(
            warthog::cpd::rle_run32{ (head << 4) | firstmove} );

//    std::cerr << "compressed source row " << source_id << " with "
//        << fm_.at(source_id).size() << std::endl;
}

std::istream&
warthog::cpd::operator>>(std::istream& in,
        warthog::cpd::graph_oracle& lab)
{
    // read the graph size data
    warthog::timer mytimer;
    mytimer.start();

    uint32_t num_nodes;
    in.read((char*)(&num_nodes), 4);
    // Need to check whether we have initialized the graph as serialising
    // removes the internal pointer.
    if(lab.g_ != nullptr && num_nodes != lab.g_->get_num_nodes())
    {
        std::cerr
            << "err; " << "input mismatch. cpd file says " << num_nodes
            << " nodes, but graph contains " << lab.g_->get_num_nodes() << "\n";
        return in;
    }

    lab.fm_.clear();
    lab.order_.resize(num_nodes);

    // read the vertex-to-column-order mapping
    for(uint32_t i = 0; i < num_nodes; i++)
    {
        uint32_t n_id;
        in.read((char*)(&n_id), 4);
        lab.order_.at(i) = n_id;
    }

    // read the RLE data
    uint32_t run_count = 0;
    lab.fm_.resize(num_nodes);
    for(uint32_t row_id = 0; row_id < num_nodes; row_id++)
    {
        // Check if we have a partial CPD file
        if (in.peek() == EOF)
        {
            lab.fm_.resize(row_id);
            std::cerr << "early stop; ";
            break;
        }
        // number of runs for this row
        uint32_t num_runs;
        in.read((char*)(&num_runs), 4);

        // read all the runs for the current row
        for(uint32_t i = 0; i < num_runs; i++)
        {
            warthog::cpd::rle_run32 tmp;
            in >> tmp;
            lab.fm_.at(row_id).push_back(tmp);
            run_count++;

            if(!in.good())
            {
                std::cerr << "err; while reading firstmove labels\n";
                std::cerr
                    << "[debug info] "
                    << " row_id " << row_id
                    << " run# " << i << " of " << lab.fm_.size()
                    << ". aborting.\n";
                return in;
            }
        }
    }
    mytimer.stop();

    std::cerr
        << "read from disk " << lab.fm_.size()
        << " rows and "
        << run_count << " runs. "
        << " time: " << (double)mytimer.elapsed_time_nano() / 1e9 << " s\n";
    return in;
}

std::ostream&
warthog::cpd::operator<<(std::ostream& out,
        warthog::cpd::graph_oracle& lab)
{
    warthog::timer mytimer;
    mytimer.start();

    // write graph size
    uint32_t num_nodes = lab.g_->get_num_nodes();
    out.write((char*)(&num_nodes), 4);

    // write node ordering
    assert(lab.order_.size() == num_nodes);
    for(uint32_t i = 0; i < num_nodes; i++)
    {
        uint32_t n_id = lab.order_.at(i);
        out.write((char*)(&n_id), 4);
    }

    // write the runs for each row
    uint32_t row_count = 0;
    uint32_t run_count = 0;
    for(uint32_t row_id = 0; row_id < lab.g_->get_num_nodes(); row_id++)
    {
        // write the number of runs
        uint32_t num_runs = (uint32_t)lab.fm_.at(row_id).size();
        // Skip empty runs
        if (num_runs == 0) { continue; }

        out.write((char*)(&num_runs), 4);
        row_count++;

        for(uint32_t run = 0; run < num_runs; run++)
        {
            out << lab.fm_.at(row_id).at(run);
            run_count++;
            if(!out.good())
            {
                std::cerr << "err; while writing labels\n";
                std::cerr
                    << "[debug info] "
                    << " row_id " << row_id
                    << " run# " << lab.fm_.at(row_id).size()
                    << ". aborting.\n";
                return out;
            }
        }

    }
    mytimer.stop();

    std::cerr
        << "wrote to disk " << row_count
        << " rows and "
        << run_count << " runs. "
        << " time: " << (double)mytimer.elapsed_time_nano() / 1e9 << " s \n";
    return out;
}

/**
 * Append operator for CPDs. Used when building partial CPDs so we can join them
 * into a single one.
 *
 * Works by appending runs to the current instance and, if the order is not set,
 * copying the order.
 */
warthog::cpd::graph_oracle&
warthog::cpd::graph_oracle::operator+=(const warthog::cpd::graph_oracle &cpd)
{
    fm_.insert(fm_.end(), cpd.fm_.begin(), cpd.fm_.end());

    // Order is always read completely
    if (order_.size() == 0)
    {
        order_ = cpd.order_;
    }

    assert(order_ == cpd.order_);

    return *this;
}
