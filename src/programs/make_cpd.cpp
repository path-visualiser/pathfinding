/**
 * This file is used to create CPDs in an independent fashion.
 */
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <omp.h>

#include "cfg.h"
#include "graph_oracle.h"
#include "xy_graph.h"

using namespace std;

int
make_cpd(std::string xy_filename, int from, int to)
{
    warthog::graph::xy_graph g;
    warthog::cpd::graph_oracle cpd(&g);
    std::ifstream ifs(xy_filename);
    std::vector<uint32_t> order;

    if (!ifs.good())
    {
        std::cerr << "Cannot open file " << xy_filename << std::endl;
        return 1;
    }

    ifs >> g;
    ifs.close();

    uint32_t node_count = g.get_num_nodes();
    if (to < 0)
    {
        to = node_count;
    }

    warthog::timer t;
    t.start();

    std::cerr << "Computing node ordering" << std::endl;
    warthog::cpd::compute_dfs_preorder(&g, &order);

    std::cerr << "Computing Dijkstra labels" << std::endl;
#pragma omp parallel
    {
        int thread_count = omp_get_num_threads();
        int thread_id = omp_get_thread_num();

        // warthog ids are 0-indexed, so no need to do anything.
        warthog::sn_id_t
            node_begin = from + (node_count * thread_id) / thread_count;
        warthog::sn_id_t
            node_end = from + (node_count * (thread_id + 1)) / thread_count;
        warthog::sn_id_t source_id = node_begin;

        std::vector<warthog::cpd::fm_coll> s_row(node_count);
        // each thread has its own copy of Dijkstra and each
        // copy has a separate memory pool
        warthog::simple_graph_expansion_policy expander(&g);
        warthog::zero_heuristic h;
        warthog::pqueue_min queue;
        warthog::cpd::graph_oracle_listener listener;

        warthog::flexible_astar<
            warthog::zero_heuristic,
            warthog::simple_graph_expansion_policy,
            warthog::pqueue_min,
            warthog::cpd::graph_oracle_listener>
                dijk(&h, &expander, &queue, &listener);

        listener.oracle_ = &cpd;
        listener.source_id_ = &source_id;
        listener.s_row_ = &s_row;
        dijk.set_listener(&listener);

        while (source_id < node_end)
        {
            // No idea why we want to play with an int pointer...
            source_id++;
            warthog::cpd::compute_row(source_id, &cpd, &dijk, s_row);
        }
    }

    // convert the column order into a map: from vertex id to its ordered index
    warthog::helpers::value_index_swap_array(order);

    t.stop();
    std::cerr
        << "total preproc time (seconds): "
        << t.elapsed_time_sec() << "\n";

    return 0;
}

int
main(int argc, char *argv[])
{
    int from = 0;
    int to = -1;

    warthog::util::param valid_args[] =
    {
        {"from", required_argument, 0, 1},
        {"to", required_argument, 0, 1},
        {"input", required_argument, 0, 1},
        {0, 0, 0, 0}
    };

    warthog::util::cfg cfg;
    cfg.parse_args(argc, argv, valid_args);

    std::string s_from = cfg.get_param_value("from");
    std::string s_to = cfg.get_param_value("to");
    std::string fname = cfg.get_param_value("input");

    if (fname == "")
    {
        std::cerr << "Require argument --input missing." << std::endl;
        return 1;
    }

    if (s_from != "")
    {
        from = std::stoi(s_from);
    }

    if (s_to != "")
    {
        to = std::stoi(s_to);
    }

    return make_cpd(fname, from, to);
}
