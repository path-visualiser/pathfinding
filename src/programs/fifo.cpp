#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <getopt.h>
#include <csignal>
#include <iostream>
#include <fstream>
#include <vector>
#include <omp.h>

#define FIFO "/tmp/warthog.fifo"

#ifndef NDEBUG
#define VERBOSE true
#else
#define VERBOSE false
#endif

#include "cfg.h"
#include "cpd_graph_expansion_policy.h"
#include "cpd_heuristic.h"
#include "cpd_search.h"
#include "graph_oracle.h"
#include "json.hpp"
#include "log.h"
#include "solution.h"
#include "timer.h"
#include "xy_graph.h"

using namespace std;

//
// - Forwards
//
void reader();

//
// - Definitions
//
typedef warthog::sn_id_t t_query;

typedef std::tuple<unsigned int, // Nodes expanded
                   unsigned int, // Nodes inserted
                   unsigned int, // Nodes touched
                   unsigned int, // Nodes updated
                   unsigned int, // Surplus nodes
                   unsigned int, // Total length of paths
                   double,       // Difference between perturbed and opti
                   double,       // Difference between free flow and opti
                   double,       // Difference between free flow and perturbed
                   double,       // Time spent in A*
                   double>       // Time to do the search
    t_results;

typedef struct config
{
    // Default config
    config() : hscale(1.0), fscale(1.0), time(DBL_MAX), itrs(warthog::INF32),
               k_moves(warthog::INF32), prefix(0), threads(0), verbose(VERBOSE)
    {};

    double hscale;                // Modifier for heuristic's value
    double fscale;                // Quality tolerance
    double time;
    uint32_t itrs;
    uint32_t k_moves;
    uint32_t prefix;
    unsigned char threads;
    bool verbose;
    bool debug;
} config;

void
to_json(nlohmann::json& j, const config& c)
{
    j = {
        {"hscale", c.hscale}, {"fscale", c.fscale}, {"time", c.time},
        {"itrs", c.itrs}, {"k_moves", c.k_moves}, {"prefix", c.prefix},
        {"threads", c.threads}, {"verbose", c.verbose},
        {"debug", c.debug}
    };
}

void
from_json(const nlohmann::json& j, config &c)
{
    j.at("hscale").get_to(c.hscale);
    j.at("fscale").get_to(c.fscale);
    j.at("time").get_to(c.time);
    j.at("itrs").get_to(c.itrs);
    j.at("k_moves").get_to(c.k_moves);
    j.at("prefix").get_to(c.prefix);
    j.at("threads").get_to(c.threads);
    j.at("verbose").get_to(c.verbose);
    j.at("debug").get_to(c.debug);
}

std::ostream&
operator<<(std::ostream& os, config &c)
{
    nlohmann::json j = c;

    os << j;

    return os;
}

config&
operator>>(std::istream& is, config &c)
{
    nlohmann::json j;

    is >> j;
    c = j.get<config>(); // this will only read a complete object

    return c;
}

/**
 * Takes care of "default parameters" as we use a bunch of wildcards to
 * represent different unbounded values.
 *
 * TODO Should this be part of cpd search directly?
 */
void
sanitise_conf(config& conf)
{
    conf.fscale = std::max(1.0, conf.fscale);
    conf.hscale = std::max(1.0, conf.hscale);

    if (conf.itrs == 0)
    { conf.itrs = warthog::INF32; }

    if (conf.k_moves == 0)
    { conf.k_moves = warthog::INF32; }

    if (conf.time == 0)
    { conf.time = DBL_MAX; }

    // Do not touch conf.prefix

    // Enforce single threaded or use max threads
#ifdef SINGLE_THREADED
    conf.threads = 1;
#else
    if (conf.threads == 0)
    {
        conf.threads = omp_get_max_threads();
    }
#endif

}

// Global container, booh
vector<warthog::search*> algos;

//
// - Functions
//
void
signalHandler(int signum)
{
    warning(true, "Interrupt signal", signum, "received.");

    remove(FIFO);

    exit(signum);
}

/**
 * Helper function to run a fixed number of CPD extractions.
 */
void
extractions(warthog::sn_id_t start_id, warthog::sn_id_t target_id,
            uint32_t prefix, warthog::cpd::graph_oracle* cpd,
            warthog::graph::xy_graph* g, warthog::cost_t& base,
            warthog::cost_t &init)
{
    uint32_t p = 0;
    warthog::sn_id_t pid = start_id;

    while (p < prefix && pid != target_id)
    {
        uint32_t fm = cpd->get_move(pid, target_id);
        warthog::graph::node* n = g->get_node(pid);

        // Should not happen
        if (fm == warthog::cpd::CPD_FM_NONE)
        {
            error("Could not find path", start_id, "->", target_id);
            break;
        }

        assert(n->out_degree() > fm);
        warthog::graph::edge* e = n->outgoing_begin() + fm;
        pid = e->node_id_;
        base += warthog::cpd::label_to_wt(e->label_);
        init += e->wt_;
        p++;
    }
}

/**
 * The search function does a bunch of statistics out of the search. It takes a
 * configration object, an output pipe and a list of queries and processes them.
 */
void
cpd_search(config& conf,
           const std::string& fifo_out,
           const std::vector<t_query> &reqs)
{
    assert(reqs.size() % 2 == 0);
    size_t n_results = reqs.size() / 2;
    // Statistics
    unsigned int n_expanded = 0;
    unsigned int n_inserted = 0;
    unsigned int n_touched = 0;
    unsigned int n_updated = 0;
    unsigned int n_surplus = 0;
    unsigned int plen = 0;
    double t_astar = 0;
    double cost_diff = 0;
    double flow_diff = 0;
    double base_diff = 0;

    warthog::timer t;
    user(conf.verbose, "Preparing to process", n_results, "queries using",
         (int) conf.threads, "threads.");

    t.start();

#pragma omp parallel num_threads(conf.threads)                          \
    reduction(+ : t_astar, n_expanded, n_inserted, n_touched, n_updated, \
              n_surplus, cost_diff, flow_diff, base_diff)
    {
        // Parallel data
        const int thread_count = omp_get_num_threads();
        const int thread_id = omp_get_thread_num();

        warthog::timer t_thread;
        warthog::solution sol;

        warthog::cpd_search<
            warthog::cpd_heuristic,
            warthog::simple_graph_expansion_policy,
            warthog::pqueue_min>* alg = static_cast<warthog::cpd_search<
                warthog::cpd_heuristic,
                warthog::simple_graph_expansion_policy,
                warthog::pqueue_min>*>(algos.at(thread_id));

        // Setup algo's config; we assume sane inputs
        alg->get_heuristic()->set_hscale(conf.hscale);
        alg->set_max_time_cutoff(conf.time); // This needs to be in ns
        alg->set_max_expansions_cutoff(conf.itrs);
        alg->set_max_k_moves(conf.k_moves);
        alg->set_quality_cutoff(conf.fscale);

        // Instead of bothering with manual conversiong (think 'ceil()'), we use
        // the magic of "usual arithmetic" to achieve the right from/to values.
        size_t step = n_results * thread_id, from = step / thread_count,
            to = (step + n_results) / thread_count;

        t_thread.start();
        // Iterate over the *requests* then convert to ids ({o,d} pair)
        for (auto id = from; id < to; id += 1)
        {
            size_t i = id * 2;
            warthog::sn_id_t start_id = reqs.at(i);
            warthog::sn_id_t target_id = reqs.at(i + 1);
            warthog::cost_t base = 0;
            warthog::cost_t init = 0;

            if (conf.prefix > 0)
            {
                extractions(start_id, target_id, conf.prefix,
                            alg->get_heuristic()->get_oracle(),
                            alg->get_expander()->get_g(), base, init);
            }
            else
            {
                // Initial path
                alg->get_heuristic()->h(start_id, target_id, base, init);
                // Actual search
                warthog::problem_instance pi(start_id, target_id, conf.debug);
                alg->get_pathcost(pi, sol);
            }

            // Update stasts
            t_astar += sol.time_elapsed_nano_;
            n_expanded += sol.nodes_expanded_;
            n_inserted += sol.nodes_inserted_;
            n_touched += sol.nodes_touched_;
            n_updated += sol.nodes_updated_;
            n_surplus += sol.nodes_surplus_;
            plen += sol.path_.size();
            cost_diff += (init - sol.sum_of_edge_costs_) / init;
            flow_diff += (sol.sum_of_edge_costs_ - base) / sol.sum_of_edge_costs_;
            base_diff += (init - base) / init;
        }

        t_thread.stop();

#pragma omp critical
        trace(conf.verbose, "[", thread_id, "] Processed", to - from,
              "trips in", t_thread.elapsed_time_micro(), "us.");
    }

    t.stop();

    user(conf.verbose, "Processed", n_results, "in", t.elapsed_time_micro(),
         "us");

    std::ofstream out(fifo_out);

    debug(conf.verbose, "Spawned a writer on", fifo_out);
    out << n_expanded << "," << n_inserted << "," << n_touched << ","
        << n_updated << "," << n_surplus << "," << plen << ","
        << cost_diff << "," << flow_diff << "," << base_diff << ","
        << t_astar << "," << t.elapsed_time_nano()
        << std::endl;

    out.close();
}

/**
 * The reader thread reads the data passed to the pipe ('FIFO') in the following
 * order:
 *
 *  1. the configuration for the search;
 *
 *  2. the output pipe's name and the number of queries; and,
 *
 *  3. the queries as (o, d)-pairs.
 *
 * It then passes the data to the search function before calling itself again.
 */
void
reader()
{
    ifstream fd;
    config conf;
    string fifo_out;
    vector<t_query> lines;

    while (true)
    {
        fd.open(FIFO);
        debug(VERBOSE, "waiting for writers...");

        if (fd.good())
        {
            debug(VERBOSE, "Got a writer");
        }
        // else?

        // Start by reading config
        try
        {
            fd >> conf;
            sanitise_conf(conf);
        } // Ignore bad parsing and fall back on default conf
        catch (std::exception& e)
        {
            debug(conf.verbose, e.what());
        }

        trace(conf.verbose, conf);

        // Read output pipe and size of input
        size_t s = 0;
        fd >> fifo_out >> s;
        debug(conf.verbose, "Preparing to read", s, "items.");
        debug(conf.verbose, "Output to", fifo_out);

        warthog::sn_id_t o, d;
        size_t i = 0;

        lines.resize(s * 2);
        while (fd >> o >> d)
        {
            lines.at(i) = o;
            lines.at(i + 1) = d;
            i += 2;
        }
        fd.close();                 // TODO check if we need to keep this open

        trace(conf.verbose, "Read", int(lines.size() / 2), "queries.");
        assert(lines.size() == s * 2);

        DO_ON_DEBUG_IF(conf.debug)
        {
            for (size_t i = 0; i < lines.size(); i += 2)
            {
                // Not using `verbose` here, it's a lot of info...
                debug(conf.debug, lines.at(i), ",", lines.at(i + 1));
            }
        }

        if (lines.size() > 0)
        {
            cpd_search(conf, fifo_out, lines);
        }
    }
}
/**
 * The main takes care of loading the data and spawning the reader thread.
 */
int
main(int argc, char *argv[])
{
	// parse arguments
	warthog::util::param valid_args[] =
        {
            // {"help", no_argument, &print_help, 1},
            // {"checkopt",  no_argument, &checkopt, 1},
            // {"verbose",  no_argument, &verbose, 1},
            // {"noheader",  no_argument, &suppress_header, 1},
            {"input",  required_argument, 0, 1},
            // {"problem",  required_argument, 0, 1},
            {0,  0, 0, 0}
        };

	warthog::util::cfg cfg;
    bool label_as_lb = false;
    warthog::graph::xy_graph g;
    std::ifstream ifs;

	cfg.parse_args(argc, argv, "-f", valid_args);

    // TODO
    // if(argc == 1 || print_help)
    // {
	// 	help();
    //     exit(0);
    // }

    // We first load the xy_graph and its diff as we need them to be *read* in
    // reverse order.
    std::string xy_filename = cfg.get_param_value("input");
    if(xy_filename == "")
    {
        std::cerr << "parameter is missing: --input [xy-graph file]\n";
        return EXIT_FAILURE;
    }

    // Check if we have a second parameter in the --input
    std::string diff_filename = cfg.get_param_value("input");
    if (diff_filename == "")
    {
        diff_filename = xy_filename + ".diff";
    }

    ifs.open(diff_filename);
    if (!ifs.good())
    {
        ifs.close();
        ifs.open(xy_filename);
        ifs >> g;
        // No modification on the graph
        label_as_lb = false;
    }
    else
    {
        // Perturbed graph first
        ifs >> g;
        ifs.close();
        // Base graph second
        ifs.open(xy_filename);
        g.perturb(ifs);
        label_as_lb = true;
    }
    ifs.close();

    // read the cpd (create from scratch if one doesn't exist)
    warthog::cpd::graph_oracle oracle(&g);
    std::string cpd_filename = cfg.get_param_value("input");
    if(cpd_filename == "")
    {
        cpd_filename = xy_filename + ".cpd";
    }

    ifs.open(cpd_filename);
    if(ifs.is_open())
    {
        ifs >> oracle;
        ifs.close();
    }
    else
    {
        oracle.precompute();
        std::ofstream ofs(cpd_filename);
        ofs << oracle;
        ofs.close();
    }

#ifdef SINGLE_THREADED
    algos.resize(1);
#else
    algos.resize(omp_get_max_threads());
#endif

    for (auto& alg: algos)
    {
        warthog::simple_graph_expansion_policy* expander =
            new warthog::simple_graph_expansion_policy(&g);
        warthog::cpd_heuristic* h =
            new warthog::cpd_heuristic(&oracle, 1.0, label_as_lb);
        warthog::pqueue_min* open = new warthog::pqueue_min();

        alg = new warthog::cpd_search<
            warthog::cpd_heuristic,
            warthog::simple_graph_expansion_policy,
            warthog::pqueue_min>(h, expander, open);
    }

    user(VERBOSE, "Loaded", algos.size(), "search.");

    int status = mkfifo(FIFO, S_IFIFO | 0666);

    if (status < 0)
    {
        perror("mkfifo");
        return EXIT_FAILURE;
    }

    // Register signal handlers
    signal(SIGINT, signalHandler);
    signal(SIGTERM, signalHandler);
    reader();

    // We do not exit from here
    return EXIT_FAILURE;
}
