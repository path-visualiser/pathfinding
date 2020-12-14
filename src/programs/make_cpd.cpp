/**
 * This file is used to create CPDs in an independent fashion.
 */
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <omp.h>

#include "bidirectional_graph_expansion_policy.h"
#include "cfg.h"
#include "graph_oracle.h"
#include "oracle_listener.h"
#include "log.h"
#include "xy_graph.h"

/**
 * Rebuild a CPD given a list of file containing its parts.
 *
 * The partial CPDs must be given in the order of the nodes.
 */
int
join_cpds(warthog::graph::xy_graph &g, std::string cpd_filename,
          std::vector<std::string> file_list, bool verbose)
{
    // Here, the type of the oracle does not matter.
    warthog::cpd::graph_oracle cpd(&g);
    cpd.compute_dfs_preorder();
    // convert the column order into a map: from vertex id to its ordered index
    cpd.value_index_swap_array();

    for (auto name: file_list)
    {
        warthog::cpd::graph_oracle part(&g);
        std::ifstream ifs(name);

        if (!ifs.good())
        {
            std::cerr << "Cannot open file " << name << std::endl;
            return EXIT_FAILURE;
        }

        ifs >> part;
        ifs.close();

        cpd.append_fm(part);
    }

    std::ofstream ofs(cpd_filename);

    if (!ofs.good())
    {
        std::cerr << "Could not open CPD file " << cpd_filename << std::endl;
        return EXIT_FAILURE;
    }

    info(verbose, "Writing results to", cpd_filename);
    ofs << cpd;
    ofs.close();

    return EXIT_SUCCESS;
}

template<warthog::cpd::symbol S>
int
make_cpd(warthog::graph::xy_graph &g, warthog::cpd::graph_oracle_base<S> &cpd,
         std::vector<warthog::cpd::oracle_listener*> &listeners,
         std::string cpd_filename, bool reverse, int from, int to,
         bool verbose=false)
{
    uint32_t node_count = g.get_num_nodes();
    if (to < 0 || node_count < (uint32_t) to)
    {
        to = node_count;
    }

    assert(to > 0);
    assert(from < to);
    assert((unsigned int)from < node_count);
    assert((unsigned int)to <= node_count);

    uint32_t node_end = to;
    unsigned char pct_done = 0;
    uint32_t nprocessed = 0;

    warthog::timer t;
    t.start();

    info(verbose, "Computing node ordering.");
    cpd.compute_dfs_preorder();

    info(verbose, "Computing Dijkstra labels.");
    std::cerr << "progress: [";
    for(uint32_t i = 0; i < 100; i++) { std::cerr <<" "; }
    std::cerr << "]\rprogress: [";

    #ifndef SINGLE_THREADED
    #pragma omp parallel
    #endif
    {
        int thread_count = omp_get_num_threads();
        int thread_id = omp_get_thread_num();

        // warthog ids are 0-indexed, so no need to do anything.
        warthog::sn_id_t source_id = from + thread_id;

        std::vector<warthog::cpd::fm_coll> s_row(node_count);
        // each thread has its own copy of Dijkstra and each
        // copy has a separate memory pool
        warthog::bidirectional_graph_expansion_policy expander(&g, reverse);
        warthog::zero_heuristic h;
        warthog::pqueue_min queue;
        warthog::flexible_astar<
            warthog::zero_heuristic,
            warthog::bidirectional_graph_expansion_policy,
            warthog::pqueue_min,
            warthog::cpd::oracle_listener>
            dijk(&h, &expander, &queue);

        listeners.at(thread_id)->set_run(&source_id, &s_row);
        dijk.set_listener(listeners.at(thread_id));

        while (source_id < node_end)
        {
            cpd.compute_row(source_id, &dijk, s_row);
            // We increment the source by the number of threads to *jump* to
            // that id.
            source_id += thread_count;
            #pragma omp critical
            {
                nprocessed++;

                if ((nprocessed * 100 / (to - from)) > pct_done)
                {
                    std::cerr << "=";
                    pct_done++;
                }
            }
        }
    }

    std::cerr << std::endl;
    // convert the column order into a map: from vertex id to its ordered index
    cpd.value_index_swap_array();

    t.stop();
    info(verbose, "total preproc time (seconds):", t.elapsed_time_sec());

    std::ofstream ofs(cpd_filename);

    if (!ofs.good())
    {
        std::cerr << "Could not open CPD file " << cpd_filename << std::endl;
        return EXIT_FAILURE;
    }

    info(verbose, "Writing results to", cpd_filename);
    ofs << cpd;
    ofs.close();

    for (auto l : listeners)
    {
        delete l;
    }

    return EXIT_SUCCESS;
}

int
main(int argc, char *argv[])
{
    int verbose = 0;
    warthog::util::param valid_args[] =
    {
        {"from", required_argument, 0, 1},
        {"to", required_argument, 0, 1},
        {"input", required_argument, 0, 1},
        {"output", required_argument, 0, 1},
        {"join", required_argument, 0, 1},
        {"type", required_argument, 0, 1},
        {"verbose", no_argument, &verbose, 1},
        {0, 0, 0, 0}
    };

    warthog::util::cfg cfg;
    cfg.parse_args(argc, argv, valid_args);

    std::string type = cfg.get_param_value("type");
    bool reverse;
    warthog::cpd::symbol cpd_type;

    if (type == "" || type == "fwd" || type == "forward")
    {
        cpd_type = warthog::cpd::FORWARD;
        reverse = false;
    }
    else if (type == "bwd" || type == "backward" || type == "rev" ||
             type == "reverse")
    {
        cpd_type = warthog::cpd::REVERSE;
        reverse = true;
    }
    else if (type == "bearing")
    {
        cpd_type = warthog::cpd::BEARING;
        reverse = true;
        // TODO Implement forward azimuth compression
    }
    else if (type == "table")
    {
        cpd_type = warthog::cpd::TABLE;
        reverse = true;
    }
    else
    {
        std::cerr << "Unknown CPD type '" << type << "'\n";
        return EXIT_FAILURE;
    }

    std::string s_from = cfg.get_param_value("from");
    std::string s_to = cfg.get_param_value("to");
    std::string xy_filename = cfg.get_param_value("input");
    std::string cpd_filename = cfg.get_param_value("output");

    if (xy_filename == "")
    {
        std::cerr << "Required argument --input [xy graph] missing."
                  << std::endl;
        return EXIT_FAILURE;
    }

    // We save the incoming edges in case we are building a reverse CPD
    warthog::graph::xy_graph g(0, "", reverse);
    std::ifstream ifs(xy_filename);

    if (!ifs.good())
    {
        std::cerr << "Cannot open file " << xy_filename << std::endl;
        return EXIT_FAILURE;
    }

    ifs >> g;
    ifs.close();

    int from = 0;
    int to = -1;

    if (s_from != "")
    {
        from = std::stoi(s_from);

        if (from < 0)
        {
            std::cerr << "Argument --from [node id] cannot be negative."
                      << std::endl;
            return EXIT_FAILURE;
        }
    }

    if (s_to != "")
    {
        to = std::stoi(s_to);
    }

    if (cpd_filename == "")
    {
        // Use default name
        cpd_filename = xy_filename;

        switch (cpd_type)
        {
            case warthog::cpd::REVERSE:
                cpd_filename += "-rev";
            case warthog::cpd::BEARING:
                cpd_filename += "-bearing";
                break;
            case warthog::cpd::TABLE:
                cpd_filename += "-table";
                break;
            default: // noop
                break;
        }

        cpd_filename += ".cpd";

        std::cerr << "No --output provided, defaulting to: " << cpd_filename
                  << std::endl;
    }

    if (cfg.get_num_values("join") > 0)
    {
        std::vector<std::string> names;
        std::string part;

        while (true)
        {
            part = cfg.get_param_value("join");

            if (part == "") { break; }
            names.push_back(part);
        }

        return join_cpds(g, cpd_filename, names, verbose);
    }
    else
    {
        #ifdef SINGLE_THREADED
        size_t nthreads = 1;
        #else
        size_t nthreads = omp_get_max_threads();
        #endif
        std::vector<warthog::cpd::oracle_listener*> listeners(nthreads);

        // We have to explicitly create and pass the different (sub-) types of
        // oracles and listeners or it messes with the template resolution.
        switch (cpd_type)
        {
            case warthog::cpd::REVERSE:
            {
                warthog::cpd::graph_oracle_base<warthog::cpd::REVERSE> cpd(&g);

                for (size_t t = 0; t < nthreads; t++)
                {
                    listeners.at(t) =
                        new warthog::cpd::reverse_oracle_listener(&cpd);
                }

                return make_cpd<warthog::cpd::REVERSE>(
                    g, cpd, listeners, cpd_filename, reverse, from, to, verbose);
            }

            case warthog::cpd::BEARING:
            {
                warthog::cpd::graph_oracle_base<warthog::cpd::BEARING> cpd(&g);

                for (size_t t = 0; t < nthreads; t++)
                {
                    listeners.at(t) =
                        new warthog::cpd::reverse_bearing_oracle_listener(&cpd);
                }

                return make_cpd<warthog::cpd::BEARING>(
                    g, cpd, listeners, cpd_filename, reverse, from, to, verbose);
            }

            case warthog::cpd::TABLE:
            {
                warthog::cpd::graph_oracle_base<warthog::cpd::TABLE> cpd(&g);

                for (size_t t = 0; t < nthreads; t++)
                {
                    listeners.at(t) =
                        new warthog::cpd::table_oracle_listener(&cpd);
                }

                return make_cpd<warthog::cpd::TABLE>(
                    g, cpd, listeners, cpd_filename, reverse, from, to, verbose);
            }

            // case warthog::cpd::FORWARD:
            default:
            {
                warthog::cpd::graph_oracle cpd(&g);

                for (size_t t = 0; t < nthreads; t++)
                {
                    listeners.at(t) =
                        new warthog::cpd::graph_oracle_listener(&cpd);
                }

                return make_cpd<warthog::cpd::FORWARD>(
                    g, cpd, listeners, cpd_filename, reverse, from, to, verbose);
            }
        }
    }
}
