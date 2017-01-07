// programs/roadhog.cpp 
//
// Pulls together a variety of different algorithms for 
// routing on road graphs.
//
// @author: dharabor
// @created: 2016-11-24
//

#include "af_filter.h"
#include "bbaf_filter.h"
#include "bb_filter.h"
#include "bidirectional_search.h"
#include "cfg.h"
#include "ch_expansion_policy.h"
#include "chaf_expansion_policy.h"
#include "chafbb_expansion_policy.h"
#include "chase_search.h"
#include "constants.h"
#include "contraction.h"
#include "corner_point_graph.h"
#include "dcl_filter.h"
#include "dimacs_parser.h"
#include "down_distance_filter.h"
#include "dummy_filter.h"
#include "euclidean_heuristic.h"
#include "fch_af_expansion_policy.h"
#include "fch_bb_expansion_policy.h"
#include "fch_bbaf_expansion_policy.h"
#include "fch_dcl_expansion_policy.h"
#include "fch_dd_expansion_policy.h"
#include "fch_expansion_policy.h"
#include "fixed_graph_contraction.h"
#include "flexible_astar.h"
#include "graph_expansion_policy.h"
#include "lazy_graph_contraction.h"
#include "planar_graph.h"
#include "timer.h"
#include "zero_heuristic.h"

// stuff for searching corner point graphs
#include "ch_cpg_expansion_policy.h"
#include "fch_af_cpg_expansion_policy.h"
#include "fch_bb_cpg_expansion_policy.h"
#include "fch_bbaf_cpg_expansion_policy.h"
#include "fch_cpg_expansion_policy.h"

#include "getopt.h"

#include <fstream>
#include <functional>
#include <iomanip>
#include <memory>
#include <sstream>
#include <tr1/unordered_map>

// check computed solutions are optimal
int checkopt = 0;
// print debugging info during search
int verbose = 0;
// display program help on startup
int print_help = 0;

void
help()
{
	std::cerr << "Valid parameters:\n"
    << "\t--alg [ algorithm name (required) ]\n"
    << "\t--input [ algorithm-specific input files (omit to show options) ] \n"
    << "\t--problem [ ss or p2p problem file (required) ]\n"
	<< "\t--verbose (optional)\n"
    << "\nRecognised values for --alg:\n"
    << "\tastar, dijkstra, bi-astar, bi-dijkstra\n"
    << "\tch, chase, chaf, chaf-bb, ch-cpg\n"
    << "\tfch, fch-dd, fch-af, fch-bb, fch-bbaf, fch-dcl\n"
    << "\tfch-cpg, fch-bb-cpg\n"
    << "\nRecognised values for --input:\n "
    << "\ttoo many to list. missing input files will be listed at uuntime\n";
}

void
run_astar(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
        std::string alg_name, std::string gr, std::string co)
{
    warthog::graph::planar_graph g;
    g.load_dimacs(gr.c_str(), co.c_str());

    warthog::graph_expansion_policy expander(&g);

    warthog::euclidean_heuristic h(&g);
    warthog::flexible_astar<
        warthog::euclidean_heuristic, 
        warthog::graph_expansion_policy> alg(&h, &expander);
    alg.set_verbose(verbose);

    int i = 0;
    std::cout << "id\talg\texp\tgen\ttouch\tmicros\tplen\tmap\n";
    for(warthog::dimacs_parser::experiment_iterator it = parser.experiments_begin(); 
            it != parser.experiments_end(); it++)
    {
        warthog::dimacs_parser::experiment exp = (*it);
        double len = alg.get_length(warthog::problem_instance(
                    exp.source, (exp.p2p ? exp.target : warthog::INF)));

        std::cout << i++ <<"\t" << alg_name << "\t" 
        << alg.get_nodes_expanded() << "\t" 
        << alg.get_nodes_generated() << "\t"
        << alg.get_nodes_touched() << "\t"
        << alg.get_search_time()  << "\t"
        << len << "\t" 
        << gr << " " << parser.get_problemfile() << std::endl;
    }
}

void
run_bi_astar(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
        std::string alg_name, std::string gr, std::string co)
{
    warthog::graph::planar_graph g;
    g.load_dimacs(gr.c_str(), co.c_str(), false, true);
    warthog::graph_expansion_policy fexp(&g);

    warthog::graph::planar_graph backward_g;
    backward_g.load_dimacs(gr.c_str(), co.c_str(), true, true);
    warthog::graph_expansion_policy bexp(&g);

    warthog::euclidean_heuristic h(&g);
    warthog::bidirectional_search<warthog::euclidean_heuristic>
        alg(&fexp, &bexp, &h);
    alg.set_verbose(verbose);

    int i = 0;
    std::cout << "id\talg\texp\tgen\ttouch\tmicros\tplen\tmap\n";
    for(warthog::dimacs_parser::experiment_iterator it = parser.experiments_begin(); 
            it != parser.experiments_end(); it++)
    {
        warthog::dimacs_parser::experiment exp = (*it);
        double len = alg.get_length(warthog::problem_instance(
                    exp.source, (exp.p2p ? exp.target : warthog::INF)));

        std::cout << i++ <<"\t" << alg_name << "\t" 
        << alg.get_nodes_expanded() << "\t" 
        << alg.get_nodes_generated() << "\t"
        << alg.get_nodes_touched() << "\t"
        << alg.get_search_time()  << "\t"
        << len << "\t" 
        << gr << " " << parser.get_problemfile() << std::endl;
    }
}

void
run_bi_dijkstra(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
        std::string alg_name, std::string gr, std::string co)
{
    warthog::graph::planar_graph g;
    g.load_dimacs(gr.c_str(), co.c_str());
    warthog::graph_expansion_policy fexp(&g);

    warthog::graph::planar_graph backward_g;
    backward_g.load_dimacs(gr.c_str(), co.c_str(), true, true);
    warthog::graph_expansion_policy bexp(&backward_g);

    warthog::zero_heuristic h;
    warthog::bidirectional_search<warthog::zero_heuristic>
        alg(&fexp, &bexp, &h);
    alg.set_verbose(verbose);

    int i = 0;
    std::cout << "id\talg\texp\tgen\ttouch\tmicros\tplen\tmap\n";
    for(warthog::dimacs_parser::experiment_iterator it = parser.experiments_begin(); 
            it != parser.experiments_end(); it++)
    {
        warthog::dimacs_parser::experiment exp = (*it);
        double len = alg.get_length(warthog::problem_instance(
                    exp.source, (exp.p2p ? exp.target : warthog::INF)));

        std::cout << i++ <<"\t" << alg_name << "\t" 
        << alg.get_nodes_expanded() << "\t" 
        << alg.get_nodes_generated() << "\t"
        << alg.get_nodes_touched() << "\t"
        << alg.get_search_time()  << "\t"
        << len << "\t" 
        << gr << " " << parser.get_problemfile() << std::endl;
    }
}

void
run_dijkstra(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
        std::string alg_name, std::string gr, std::string co)
{
    warthog::graph::planar_graph g;
    g.load_dimacs(gr.c_str(), co.c_str());
    warthog::graph_expansion_policy expander(&g);

    warthog::zero_heuristic h;
    warthog::flexible_astar<warthog::zero_heuristic, 
        warthog::graph_expansion_policy> alg(&h, &expander);
    alg.set_verbose(verbose);

    int i = 0;
    for(warthog::dimacs_parser::experiment_iterator it = parser.experiments_begin(); 
            it != parser.experiments_end(); it++)
    {
        warthog::dimacs_parser::experiment exp = (*it);
        double len = alg.get_length(warthog::problem_instance(
                    exp.source, (exp.p2p ? exp.target : warthog::INF)));

        std::cout << i++ <<"\t" << alg_name << "\t" 
        << alg.get_nodes_expanded() << "\t" 
        << alg.get_nodes_generated() << "\t"
        << alg.get_nodes_touched() << "\t"
        << alg.get_search_time()  << "\t"
        << len << "\t" 
        << gr << " " << parser.get_problemfile() << std::endl;
    }
}

void
run_ch(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
        std::string alg_name, std::string gr, std::string co)
{
    std::string orderfile = cfg.get_param_value("input");
    if(orderfile == "")
    {
        std::cerr << "err; missing contraction order input file\n";
        return;
    }

    // load up the graph 
    warthog::graph::planar_graph g;
    g.load_dimacs(gr.c_str(), co.c_str(), false, true);

    // load up the node order
    std::vector<uint32_t> order;
    warthog::ch::load_node_order(orderfile.c_str(), order);
    warthog::ch::value_index_swap_dimacs(order);

    std::cerr << "preparing to search\n";
    warthog::zero_heuristic h;
    //warthog::euclidean_heuristic h(&g);
    warthog::ch_expansion_policy fexp(&g, &order);
    warthog::ch_expansion_policy bexp (&g, &order, true);
    warthog::bidirectional_search<warthog::zero_heuristic> alg(&fexp, &bexp, &h);
    //warthog::bidirectional_search<warthog::euclidean_heuristic> alg(&fexp, &bexp, &h);
    alg.set_verbose(verbose);

    std::cerr << "running experiments\n";
    int i = 0;
    std::cout << "id\talg\texp\tgen\ttouch\tmicros\tplen\tmap\n";
    for(warthog::dimacs_parser::experiment_iterator it = parser.experiments_begin(); 
            it != parser.experiments_end(); it++)
    {
        warthog::dimacs_parser::experiment exp = (*it);
        double len = alg.get_length(warthog::problem_instance(
                    exp.source, (exp.p2p ? exp.target : warthog::INF)));

        std::cout << i++ <<"\t" << alg_name << "\t" 
        << alg.get_nodes_expanded() << "\t" 
        << alg.get_nodes_generated() << "\t"
        << alg.get_nodes_touched() << "\t"
        << alg.get_search_time()  << "\t"
        << len << "\t" 
        << gr << " " << parser.get_problemfile() << std::endl;
    }
}

void
run_ch_astar(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
        std::string alg_name, std::string gr, std::string co)
{
    std::string orderfile = cfg.get_param_value("input");
    if(orderfile == "")
    {
        std::cerr << "err; missing contraction order input file\n";
        return;
    }

    // load up the graph 
    warthog::graph::planar_graph g;
    g.load_dimacs(gr.c_str(), co.c_str(), false, true);

    // load up the node order
    std::vector<uint32_t> order;
    warthog::ch::load_node_order(orderfile.c_str(), order);
    warthog::ch::value_index_swap_dimacs(order);

    std::cerr << "preparing to search\n";
    warthog::euclidean_heuristic h(&g);
    warthog::ch_expansion_policy fexp(&g, &order);
    warthog::ch_expansion_policy bexp (&g, &order, true);
    warthog::bidirectional_search<warthog::euclidean_heuristic> alg(&fexp, &bexp, &h);
    alg.set_verbose(verbose);

    std::cerr << "running experiments\n";
    int i = 0;
    std::cout << "id\talg\texp\tgen\ttouch\tmicros\tplen\tmap\n";
    for(warthog::dimacs_parser::experiment_iterator it = parser.experiments_begin(); 
            it != parser.experiments_end(); it++)
    {
        warthog::dimacs_parser::experiment exp = (*it);
        double len = alg.get_length(warthog::problem_instance(
                    exp.source, (exp.p2p ? exp.target : warthog::INF)));

        std::cout << i++ <<"\t" << alg_name << "\t" 
        << alg.get_nodes_expanded() << "\t" 
        << alg.get_nodes_generated() << "\t"
        << alg.get_nodes_touched() << "\t"
        << alg.get_search_time()  << "\t"
        << len << "\t" 
        << gr << " " << parser.get_problemfile() << std::endl;
    }
}

void
run_ch_cpg(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
        std::string alg_name, std::string gr, std::string co)
{
    std::string orderfile = cfg.get_param_value("input");
    std::string gridmapfile = cfg.get_param_value("input");
    if(!orderfile.compare("") || !gridmapfile.compare(""))
    {
        std::cerr << "err; insufficient input params for alg " << alg_name 
                  << ". required (in, order) "
                  << "--input [ gr file ] [co file]"
                  << " [contraction order file] [gridmap file]\n";
        return;
    }
    
    // load up the node order
    std::vector<uint32_t> order;
    if(!warthog::ch::load_node_order(orderfile.c_str(), order, true))
    {
        std::cerr << "err; could not load contraction order file\n";
        return;
    }

    // load up the grid
    std::shared_ptr<warthog::gridmap> map(
            new warthog::gridmap(gridmapfile.c_str()));
            
    // load up the graph 
    std::shared_ptr<warthog::graph::planar_graph> pg(
            new warthog::graph::planar_graph());
    pg->load_dimacs(gr.c_str(), co.c_str(), false, true);

    // wrapper
    warthog::graph::corner_point_graph cpg(map, pg);

    // expander
	warthog::zero_heuristic h;
	warthog::ch::ch_cpg_expansion_policy fexp(&cpg, &order);
	warthog::ch::ch_cpg_expansion_policy bexp(&cpg, &order, true);
    warthog::bidirectional_search<warthog::zero_heuristic> alg(&fexp, &bexp, &h);
	alg.set_verbose(verbose);

    std::cerr << "running experiments\n";
    int i = 0;
    std::cout << "id\talg\texp\tgen\ttouch\tmicros\tplen\tmap\n";
    for(warthog::dimacs_parser::experiment_iterator it = parser.experiments_begin(); 
            it != parser.experiments_end(); it++)
    {
        warthog::dimacs_parser::experiment exp = (*it);
        uint32_t startid = map.get()->to_padded_id(exp.source);
        uint32_t targetid = map.get()->to_padded_id(exp.target);
		double len = alg.get_length(warthog::problem_instance(
                    startid, targetid));

        std::cout << i++ <<"\t" << alg_name << "\t" 
        << alg.get_nodes_expanded() << "\t" 
        << alg.get_nodes_generated() << "\t"
        << alg.get_nodes_touched() << "\t"
        << alg.get_search_time()  << "\t"
        << len << "\t" 
        << gr << " " << parser.get_problemfile() << std::endl;
    }
}

void
run_chase(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
        std::string alg_name, std::string gr, std::string co)
{
    std::string orderfile = cfg.get_param_value("order");

    // load up the graph 
    warthog::graph::planar_graph g;
    g.load_dimacs(gr.c_str(), co.c_str(), false, true);

    // load up the node order
    std::vector<uint32_t> order;
    warthog::ch::load_node_order(orderfile.c_str(), order);
    warthog::ch::value_index_swap_dimacs(order);

    std::cerr << "preparing to search\n";
    warthog::euclidean_heuristic h(&g);
    warthog::ch_expansion_policy fexp(&g, &order);
    warthog::ch_expansion_policy bexp (&g, &order, true);
    warthog::chase_search<warthog::euclidean_heuristic> alg(&fexp, &bexp, &h);
    alg.set_verbose(verbose);

    std::cerr << "running experiments\n";
    int i = 0;
    std::cout << "id\talg\texp\tgen\ttouch\tmicros\tplen\tmap\n";
    for(warthog::dimacs_parser::experiment_iterator it = parser.experiments_begin(); 
            it != parser.experiments_end(); it++)
    {
        warthog::dimacs_parser::experiment exp = (*it);
        double len = alg.get_length(warthog::problem_instance(
                    exp.source, (exp.p2p ? exp.target : warthog::INF)));

        std::cout << i++ <<"\t" << alg_name << "\t" 
        << alg.get_nodes_expanded() << "\t" 
        << alg.get_nodes_generated() << "\t"
        << alg.get_nodes_touched() << "\t"
        << alg.get_search_time()  << "\t"
        << len << "\t" 
        << gr << " " << parser.get_problemfile() << std::endl;
    }
}

void
run_chaf(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
        std::string alg_name, std::string gr, std::string co)
{
    std::string orderfile = cfg.get_param_value("input");
    std::string arclabels_file = cfg.get_param_value("input");
    std::string partition_file = cfg.get_param_value("input");
    if(orderfile == "" || arclabels_file == "" || partition_file == "")
    {
        std::cerr << "err; insufficient input parameters for --alg "
                  << alg_name << ". required, in order:\n"
                  << " --input [gr file] [co file] "
                  << " [contraction order file] [arclabels file] "
                  << " [graph partition file]\n";
        return;
    }

    // load up the node order
    std::vector<uint32_t> order;
    bool sort_by_id = true;
    if(!warthog::ch::load_node_order(orderfile.c_str(), order, sort_by_id))
    {
        std::cerr << "err; could not load node order input file\n";
        return;
    }

    // load up the node partition info
    std::vector<uint32_t> part;
    if(!warthog::helpers::load_integer_labels_dimacs(
            partition_file.c_str(), part))
    {
        std::cerr << "err; could not load graph partition input file\n";
        return;
    }

    // load up the graph 
    warthog::graph::planar_graph g;
    if(!g.load_dimacs(gr.c_str(), co.c_str(), false, true))
    {
        std::cerr << "err; could not load gr or co input files (one or both)\n";
        return;
    }

    // load up the arc-flags
    warthog::af_filter filter(&g, &part);
    if(!filter.load_labels(arclabels_file.c_str()))
    {
        std::cerr << "err; could not load arcflags file\n";
        return;
    }

    std::cerr << "preparing to search\n";
    //warthog::zero_heuristic h;
    warthog::euclidean_heuristic h(&g);
    warthog::chaf_expansion_policy fexp(&g, &order, &filter);
    warthog::chaf_expansion_policy bexp (&g, &order, &filter, true);
    //warthog::bidirectional_search<warthog::zero_heuristic> alg(&fexp, &bexp, &h);
    warthog::bidirectional_search<warthog::euclidean_heuristic> alg(&fexp, &bexp, &h);
    alg.set_verbose(verbose);

    std::cerr << "running experiments\n";
    int i = 0;
    std::cout << "id\talg\texp\tgen\ttouch\tmicros\tplen\tmap\n";
    for(warthog::dimacs_parser::experiment_iterator it = parser.experiments_begin(); 
            it != parser.experiments_end(); it++)
    {
        warthog::dimacs_parser::experiment exp = (*it);
        double len = alg.get_length(warthog::problem_instance(
                    exp.source, (exp.p2p ? exp.target : warthog::INF)));

        std::cout << i++ <<"\t" << alg_name << "\t" 
        << alg.get_nodes_expanded() << "\t" 
        << alg.get_nodes_generated() << "\t"
        << alg.get_nodes_touched() << "\t"
        << alg.get_search_time()  << "\t"
        << len << "\t" 
        << g.get_filename() << " " << parser.get_problemfile() << std::endl;
    }
}

void
run_chaf_bb(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
        std::string alg_name, std::string gr, std::string co)
{
    std::string orderfile = cfg.get_param_value("input");
    std::string arclabels_file = cfg.get_param_value("input");
    std::string partition_file = cfg.get_param_value("input");
    if(orderfile == "" || arclabels_file == "" || partition_file == "")
    {
        std::cerr << "err; insufficient input parameters for --alg "
                  << alg_name << ". required, in order:\n"
                  << " --input [gr file] [co file] "
                  << " [contraction order file] [arclabels file] "
                  << " [graph partition file]\n";
        return;
    }

    // load up the node order
    std::vector<uint32_t> order;
    bool sort_by_id = true;
    if(!warthog::ch::load_node_order(orderfile.c_str(), order, sort_by_id))
    {
        std::cerr << "err; could not load node order input file\n";
        return;
    }

    // load up the node partition info
    std::vector<uint32_t> part;
    if(!warthog::helpers::load_integer_labels_dimacs(
            partition_file.c_str(), part))
    {
        std::cerr << "err; could not load graph partition input file\n";
        return;
    }

    // load up the graph 
    warthog::graph::planar_graph g;
    if(!g.load_dimacs(gr.c_str(), co.c_str(), false, true))
    {
        std::cerr << "err; could not load gr or co input files (one or both)\n";
        return;
    }
    // load up the arc-flags
    warthog::bbaf_filter filter(&g, &part);
    if(!filter.load_labels(arclabels_file.c_str()))
    {
        std::cerr << "err; could not load arcflags file\n";
        return;
    }

    std::cerr << "preparing to search\n";
    warthog::zero_heuristic h;
    //warthog::euclidean_heuristic h(&g);
    warthog::chafbb_expansion_policy fexp(&g, &order, &filter);
    warthog::chafbb_expansion_policy bexp (&g, &order, &filter, true);
    warthog::bidirectional_search<warthog::zero_heuristic> alg(&fexp, &bexp, &h);
    //warthog::bidirectional_search<warthog::euclidean_heuristic> alg(&fexp, &bexp, &h);
    alg.set_verbose(verbose);

    std::cerr << "running experiments\n";
    int i = 0;
    std::cout << "id\talg\texp\tgen\ttouch\tmicros\tplen\tmap\n";
    for(warthog::dimacs_parser::experiment_iterator it = parser.experiments_begin(); 
            it != parser.experiments_end(); it++)
    {
        warthog::dimacs_parser::experiment exp = (*it);
        double len = alg.get_length(warthog::problem_instance(
                    exp.source, (exp.p2p ? exp.target : warthog::INF)));

        std::cout << i++ <<"\t" << alg_name << "\t" 
        << alg.get_nodes_expanded() << "\t" 
        << alg.get_nodes_generated() << "\t"
        << alg.get_nodes_touched() << "\t"
        << alg.get_search_time()  << "\t"
        << len << "\t" 
        << g.get_filename() << " " << parser.get_problemfile() << std::endl;
    }
}

void
run_fch(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
        std::string alg_name, std::string gr, std::string co)
{
    std::string orderfile = cfg.get_param_value("input");
    if(orderfile == "")
    {
        std::cerr << "err; insufficient input parameters for --alg "
                  << alg_name << ". required, in order:\n"
                  << " --input [gr file] [co file] "
                  << " [contraction order file] \n";
        return;
    }

    // load up the graph 
    warthog::graph::planar_graph g;
    g.load_dimacs(gr.c_str(), co.c_str(), false, true);

    // load up the node order
    std::vector<uint32_t> order;
    warthog::ch::load_node_order(orderfile.c_str(), order);
    warthog::ch::value_index_swap_dimacs(order);

    std::cerr << "preparing to search\n";
    warthog::fch_expansion_policy fexp(&g, &order); 
    warthog::euclidean_heuristic h(&g);

    warthog::flexible_astar< warthog::euclidean_heuristic, 
        warthog::fch_expansion_policy> alg(&h, &fexp);
    alg.set_verbose(verbose);

    std::cerr << "running experiments\n";
    int i = 0;
    std::cout << "id\talg\texp\tgen\ttouch\tmicros\tplen\tmap\n";
    for(warthog::dimacs_parser::experiment_iterator it = parser.experiments_begin(); 
            it != parser.experiments_end(); it++)
    {
        warthog::dimacs_parser::experiment exp = (*it);
        double len = alg.get_length(warthog::problem_instance(
                    exp.source, (exp.p2p ? exp.target : warthog::INF)));

        std::cout << i++ <<"\t" << alg_name << "\t" 
        << alg.get_nodes_expanded() << "\t" 
        << alg.get_nodes_generated() << "\t"
        << alg.get_nodes_touched() << "\t"
        << alg.get_search_time()  << "\t"
        << len << "\t" 
        << g.get_filename() << " " << parser.get_problemfile() << std::endl;
    }
}

void
run_fch_dd(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
        std::string alg_name, std::string gr, std::string co)
{
    std::string orderfile = cfg.get_param_value("input");
    std::string arclabels_file = cfg.get_param_value("input");
    if(orderfile == "" || arclabels_file == "")
    {
        std::cerr << "err; insufficient input parameters for --alg "
                  << alg_name << ". required, in order:\n"
                  << " --input [gr file] [co file] "
                  << " [contraction order file] [arclabels file]\n";
        return;
    }

    // load up the node order
    std::vector<uint32_t> order;
    bool sort_by_id = true;
    if(!warthog::ch::load_node_order(orderfile.c_str(), order, sort_by_id))
    {
        std::cerr << "err; could not load node order input file\n";
        return;
    }

    // load up the graph 
    warthog::graph::planar_graph g;
    g.load_dimacs(gr.c_str(), co.c_str(), false, true);

    // load up the arc-flags
    warthog::down_distance_filter filter(&g);
    if(!filter.load_labels(arclabels_file.c_str()))
    {
        std::cerr << "err; could not load arcflags file\n";
        return;
    }

    warthog::euclidean_heuristic h(&g);
    warthog::fch_dd_expansion_policy fexp(&g, &order, &filter, &h);

    warthog::flexible_astar< warthog::euclidean_heuristic, 
       warthog::fch_dd_expansion_policy>
    alg(&h, &fexp);
    alg.set_verbose(verbose);

    std::cerr << "preparing to search\n";
    std::function<uint32_t(warthog::search_node*)> fn_get_apex = 
        [&order] (warthog::search_node* n) -> uint32_t
        {
            while(true)
            {
                warthog::search_node* p = n->get_parent();
                if(!p || order.at(p->get_id()) < order.at(n->get_id()))
                { break; }
                n = p;
            }
            return order.at(n->get_id());
        };

    std::function<uint32_t(uint32_t)> fn_redundant_expansions = 
        [&fexp, &order, &alg] (uint32_t apex) -> uint32_t
        {
            std::set<uint32_t> tmp;
            for(uint32_t i = 0; i < order.size(); i++)
            {
                warthog::search_node* n = fexp.get_ptr(i, alg.get_search_id());
                if(n && order.at(n->get_id()) > apex)
                { 
                    tmp.insert(n->get_id());
                }
            }

            for(uint32_t j = 0; j < order.size(); j++)
            {
                warthog::search_node* n = fexp.get_ptr(j, alg.get_search_id());
                warthog::search_node* m = n;
                while(m)
                {
                    if(tmp.find(m->get_id()) != tmp.end())
                    {
                        while(tmp.find(n->get_id()) == tmp.end())
                        {
                            tmp.insert(n->get_id());
                            n = n->get_parent();
                        }
                        break;
                    }
                    m = m->get_parent();
                }
            }
            return tmp.size();
        };


    std::cerr << "running experiments\n";
    int i = 0;
    std::cout << "id\talg\texp\tgen\ttouch\tmicros\tplen\tmap\n";
    for(warthog::dimacs_parser::experiment_iterator it = parser.experiments_begin(); 
            it != parser.experiments_end(); it++)
    {
        warthog::dimacs_parser::experiment exp = (*it);
        double len = alg.get_length(warthog::problem_instance(
                    exp.source, (exp.p2p ? exp.target : warthog::INF)));

        std::cout << i++ <<"\t" << alg_name << "\t" 
        << alg.get_nodes_expanded() << "\t" 
        << alg.get_nodes_generated() << "\t"
        << alg.get_nodes_touched() << "\t"
        << alg.get_search_time()  << "\t"
        << len << "\t" 
//            << fn_redundant_expansions(
//                    fn_get_apex(
//                        fexp.get_ptr(
//                            exp.target, 
//                            alg.get_searchid()))) << "\t"
        << g.get_filename() << " " << parser.get_problemfile() << std::endl;
    }
}

void
run_fch_af(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
        std::string alg_name, std::string gr, std::string co)
{
    std::string orderfile = cfg.get_param_value("input");
    std::string arclabels_file = cfg.get_param_value("input");
    std::string partition_file = cfg.get_param_value("input");
    if(orderfile == "" || arclabels_file == "" || partition_file == "")
    {
        std::cerr << "err; insufficient input parameters for --alg "
                  << alg_name << ". required, in order:\n"
                  << " --input [gr file] [co file] "
                  << " [contraction order file] [arclabels file] "
                  << " [graph partition file]\n";
        return;
    }

    // load up the node order
    std::vector<uint32_t> order;
    bool sort_by_id = true;
    if(!warthog::ch::load_node_order(orderfile.c_str(), order, sort_by_id))
    {
        std::cerr << "err; could not load node order input file\n";
        return;
    }

    // load up the node partition info
    std::vector<uint32_t> part;
    if(!warthog::helpers::load_integer_labels_dimacs(
            partition_file.c_str(), part))
    {
        std::cerr << "err; could not load graph partition input file\n";
        return;
    }

    // load up the graph 
    warthog::graph::planar_graph g;
    if(!g.load_dimacs(gr.c_str(), co.c_str(), false, true))
    {
        std::cerr << "err; could not load gr or co input files (one or both)\n";
        return;
    }

    // load up the arc-flags
    warthog::af_filter filter(&g, &part);
    if(!filter.load_labels(arclabels_file.c_str()))
    {
        std::cerr << "err; could not load arcflags file\n";
        return;
    }

    warthog::euclidean_heuristic h(&g);
    warthog::fch_af_expansion_policy fexp(&g, &order, &filter);

    warthog::flexible_astar< warthog::euclidean_heuristic, 
       warthog::fch_af_expansion_policy>
    alg(&h, &fexp);
    alg.set_verbose(verbose);

    std::cerr << "running experiments\n";
    int i = 0;
    std::cout << "id\talg\texp\tgen\ttouch\tmicros\tplen\tmap\n";
    for(warthog::dimacs_parser::experiment_iterator it = parser.experiments_begin(); 
            it != parser.experiments_end(); it++)
    {
        warthog::dimacs_parser::experiment exp = (*it);
        double len = alg.get_length(warthog::problem_instance(
                    exp.source, (exp.p2p ? exp.target : warthog::INF)));

        std::cout << i++ <<"\t" << alg_name << "\t" 
        << alg.get_nodes_expanded() << "\t" 
        << alg.get_nodes_generated() << "\t"
        << alg.get_nodes_touched() << "\t"
        << alg.get_search_time()  << "\t"
        << len << "\t" 
        << g.get_filename() << " " << parser.get_problemfile() << std::endl;
    }
}

void
run_fch_bb(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
        std::string alg_name, std::string gr, std::string co)
{
    std::string orderfile = cfg.get_param_value("input");
    std::string arclabels_file = cfg.get_param_value("input");
    if(orderfile == "" || arclabels_file == "")
    {
        std::cerr << "err; insufficient input parameters for --alg "
                  << alg_name << ". required, in order:\n"
                  << " --input [gr file] [co file] "
                  << " [contraction order file] [arclabels file]\n";
        return;
    }

    // load up the node order
    std::vector<uint32_t> order;
    bool sort_by_id = true;
    if(!warthog::ch::load_node_order(orderfile.c_str(), order, sort_by_id))
    {
        std::cerr << "err; could not load node order input file\n";
        return;
    }

    // load up the graph 
    warthog::graph::planar_graph g;
    if(!g.load_dimacs(gr.c_str(), co.c_str(), false, true))
    {
        std::cerr << "err; could not load gr or co input files (one or both)\n";
        return;
    }
    // load up the arc-flags
    warthog::bb_filter filter(&g);
    if(!filter.load_labels(arclabels_file.c_str()))
    {
        std::cerr << "err; could not load arcflags file\n";
        return;
    }

    warthog::euclidean_heuristic h(&g);
    warthog::fch_bb_expansion_policy fexp(&g, &order, &filter);

    warthog::flexible_astar< warthog::euclidean_heuristic, 
       warthog::fch_bb_expansion_policy>
    alg(&h, &fexp);
    alg.set_verbose(verbose);

    std::cerr << "running experiments\n";
    int i = 0;
    std::cout << "id\talg\texp\tgen\ttouch\tmicros\tplen\tmap\n";
    for(warthog::dimacs_parser::experiment_iterator it = parser.experiments_begin(); 
            it != parser.experiments_end(); it++)
    {
        warthog::dimacs_parser::experiment exp = (*it);
        double len = alg.get_length(warthog::problem_instance(
                    exp.source, (exp.p2p ? exp.target : warthog::INF)));

        std::cout << i++ <<"\t" << alg_name << "\t" 
        << alg.get_nodes_expanded() << "\t" 
        << alg.get_nodes_generated() << "\t"
        << alg.get_nodes_touched() << "\t"
        << alg.get_search_time()  << "\t"
        << len << "\t" 
        << g.get_filename() << " " << parser.get_problemfile() << std::endl;
    }
}

void
run_fch_bbaf(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
        std::string alg_name, std::string gr, std::string co)
{
    std::string orderfile = cfg.get_param_value("input");
    std::string arclabels_file = cfg.get_param_value("input");
    std::string partition_file = cfg.get_param_value("input");
    if(orderfile == "" || arclabels_file == "" || partition_file == "")
    {
        std::cerr << "err; insufficient input parameters for --alg "
                  << alg_name << ". required, in order:\n"
                  << " --input [gr file] [co file] "
                  << " [contraction order file] [arclabels file] "
                  << " [graph partition file]\n";
        return;
    }

    // load up the node order
    std::vector<uint32_t> order;
    bool sort_by_id = true;
    if(!warthog::ch::load_node_order(orderfile.c_str(), order, sort_by_id))
    {
        std::cerr << "err; could not load node order input file\n";
        return;
    }

    // load up the node partition info
    std::vector<uint32_t> part;
    if(!warthog::helpers::load_integer_labels_dimacs(
            partition_file.c_str(), part))
    {
        std::cerr << "err; could not load graph partition input file\n";
        return;
    }

    // load up the graph 
    warthog::graph::planar_graph g;
    if(!g.load_dimacs(gr.c_str(), co.c_str(), false, true))
    {
        std::cerr << "err; could not load gr or co input files (one or both)\n";
        return;
    }
    // load up the arc-flags
    warthog::bbaf_filter filter(&g, &part);
    if(!filter.load_labels(arclabels_file.c_str()))
    {
        std::cerr << "err; could not load arcflags file\n";
        return;
    }

    std::cerr << "preparing to search\n";
    warthog::euclidean_heuristic h(&g);
    warthog::fch_bbaf_expansion_policy fexp(&g, &order, &filter);

    warthog::flexible_astar< warthog::euclidean_heuristic, 
       warthog::fch_bbaf_expansion_policy>
    alg(&h, &fexp);
    alg.set_verbose(verbose);

    std::cerr << "running experiments\n";
    int i = 0;
    std::cout << "id\talg\texp\tgen\ttouch\tmicros\tplen\tmap\n";
    for(warthog::dimacs_parser::experiment_iterator it = parser.experiments_begin(); 
            it != parser.experiments_end(); it++)
    {
        warthog::dimacs_parser::experiment exp = (*it);
        double len = alg.get_length(warthog::problem_instance(
                    exp.source, (exp.p2p ? exp.target : warthog::INF)));

        std::cout << i++ <<"\t" << alg_name << "\t" 
        << alg.get_nodes_expanded() << "\t" 
        << alg.get_nodes_generated() << "\t"
        << alg.get_nodes_touched() << "\t"
        << alg.get_search_time()  << "\t"
        << len << "\t" 
        << g.get_filename() << " " << parser.get_problemfile() << std::endl;
    }
}

void
run_fch_bbaf_cpg(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
        std::string alg_name, std::string gr, std::string co)
{
    std::string orderfile = cfg.get_param_value("input");
    std::string arclabels_file = cfg.get_param_value("input");
    std::string partition_file = cfg.get_param_value("input");
    std::string gridmapfile = cfg.get_param_value("input");
    if( orderfile == "" || arclabels_file == "" || 
        partition_file == "" || gridmapfile == "" )
    {
        std::cerr << "err; insufficient input parameters for --alg "
                  << alg_name << ". required, in order:\n"
                  << " --input [gr file] [co file] "
                  << " [contraction order file] [arclabels file] "
                  << " [graph partition file] [gridmap file]\n";
        return;
    }
    
    // load up the node order
    std::vector<uint32_t> order;
    if(!warthog::ch::load_node_order(orderfile.c_str(), order, true))
    {
        std::cerr << "err; could not load contraction order file\n";
        return;
    }

    // we insert two extra elements in the event that we
    // need to insert the start or target. both have the lowest
    // possible rank in the hierarchy (0 and 1)
    // NB: along the way we need to increase all ranks by 2 
    for(uint32_t i = 0; i < order.size(); i++) order.at(i)+=2;
    order.push_back(0);
    order.push_back(1);

    // load up the node partition info
    std::vector<uint32_t> part;
    if(!warthog::helpers::load_integer_labels_dimacs(
            partition_file.c_str(), part))
    {
        std::cerr << "err; could not load graph partition input file\n";
        return;
    }

    // load up the grid
    std::shared_ptr<warthog::gridmap> map(
            new warthog::gridmap(gridmapfile.c_str()));
            
    // load up the graph 
    std::shared_ptr<warthog::graph::planar_graph> pg(
            new warthog::graph::planar_graph());
    pg->load_dimacs(gr.c_str(), co.c_str(), false, true);

    // graph wrapper
    warthog::graph::corner_point_graph cpg(map, pg);

    // load up the arc-flags
    warthog::bbaf_filter filter(pg.get(), &part);
    if(!filter.load_labels(arclabels_file.c_str()))
    {
        std::cerr << "err; could not load arcflags file\n";
        return;
    }

    warthog::euclidean_heuristic h(pg.get());
    h.set_hscale(warthog::ONE);

    warthog::fch_bbaf_cpg_expansion_policy fexp(&cpg, &order, &filter);

    warthog::flexible_astar< warthog::euclidean_heuristic, 
       warthog::fch_bbaf_cpg_expansion_policy>
    alg(&h, &fexp);
    alg.set_verbose(verbose);

    std::cerr << "running experiments\n";
    int i = 0;
    std::cout << "id\talg\texp\tgen\ttouch\tmicros\tplen\tmap\n";
    for(warthog::dimacs_parser::experiment_iterator it = parser.experiments_begin(); 
            it != parser.experiments_end(); it++)
    {
        warthog::dimacs_parser::experiment exp = (*it);
        uint32_t startid = map.get()->to_padded_id(exp.source);
        uint32_t targetid = map.get()->to_padded_id(exp.target);
		double len = alg.get_length(warthog::problem_instance(
                    startid, targetid));

        std::cout << i++ <<"\t" << alg_name << "\t" 
        << alg.get_nodes_expanded() << "\t" 
        << alg.get_nodes_generated() << "\t"
        << alg.get_nodes_touched() << "\t"
        << alg.get_search_time()  << "\t"
        << len << "\t" 
        << pg->get_filename() << " " << parser.get_problemfile() << std::endl;
    }
}

void
run_fch_dcl(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
        std::string alg_name, std::string gr, std::string co)
{
    std::string orderfile = cfg.get_param_value("input");
    std::string arclabels_file = cfg.get_param_value("input");
    if(orderfile == "" || arclabels_file == "")
    {
        std::cerr << "err; insufficient input parameters for --alg "
                  << alg_name << ". required, in order:\n"
                  << " --input [gr file] [co file] "
                  << " [contraction order file] [arclabels file]\n";
        return;
    }

    // load up the node order
    std::vector<uint32_t> order;
    bool sort_by_id = true;
    if(!warthog::ch::load_node_order(orderfile.c_str(), order, sort_by_id))
    {
        std::cerr << "err; could not load node order input file\n";
        return;
    }

    // load up the graph 
    warthog::graph::planar_graph g;
    g.load_dimacs(gr.c_str(), co.c_str(), false, true);

    // load up the arc-flags
    warthog::dcl_filter filter(&g);
    if(!filter.load_labels(arclabels_file.c_str()))
    {
        std::cerr << "err; could not load arcflags file\n";
        return;
    }

    warthog::euclidean_heuristic h(&g);
    warthog::fch_dcl_expansion_policy fexp(&g, &order, &filter);

    warthog::flexible_astar< warthog::euclidean_heuristic, 
       warthog::fch_dcl_expansion_policy>
    alg(&h, &fexp);
    alg.set_verbose(verbose);

    std::cerr << "preparing to search\n";
    std::function<uint32_t(warthog::search_node*)> fn_get_apex = 
        [&order] (warthog::search_node* n) -> uint32_t
        {
            while(true)
            {
                warthog::search_node* p = n->get_parent();
                if(!p || order.at(p->get_id()) < order.at(n->get_id()))
                { break; }
                n = p;
            }
            return order.at(n->get_id());
        };

    std::function<uint32_t(uint32_t)> fn_redundant_expansions = 
        [&fexp, &order, &alg] (uint32_t apex) -> uint32_t
        {
            std::set<uint32_t> tmp;
            for(uint32_t i = 0; i < order.size(); i++)
            {
                warthog::search_node* n = fexp.get_ptr(i, alg.get_search_id());
                if(n && order.at(n->get_id()) > apex)
                { 
                    tmp.insert(n->get_id());
                }
            }

            for(uint32_t j = 0; j < order.size(); j++)
            {
                warthog::search_node* n = fexp.get_ptr(j, alg.get_search_id());
                warthog::search_node* m = n;
                while(m)
                {
                    if(tmp.find(m->get_id()) != tmp.end())
                    {
                        while(tmp.find(n->get_id()) == tmp.end())
                        {
                            tmp.insert(n->get_id());
                            n = n->get_parent();
                        }
                        break;
                    }
                    m = m->get_parent();
                }
            }
            return tmp.size();
        };


    std::cerr << "running experiments\n";
    int i = 0;
    std::cout << "id\talg\texp\tgen\ttouch\tmicros\tplen\tmap\n";
    for(warthog::dimacs_parser::experiment_iterator it = parser.experiments_begin(); 
            it != parser.experiments_end(); it++)
    {
        warthog::dimacs_parser::experiment exp = (*it);
        double len = alg.get_length(warthog::problem_instance(
                    exp.source, (exp.p2p ? exp.target : warthog::INF)));

        std::cout << i++ <<"\t" << alg_name << "\t" 
        << alg.get_nodes_expanded() << "\t" 
        << alg.get_nodes_generated() << "\t"
        << alg.get_nodes_touched() << "\t"
        << alg.get_search_time()  << "\t"
        << len << "\t" 
//            << fn_redundant_expansions(
//                    fn_get_apex(
//                        fexp.get_ptr(
//                            exp.target, 
//                            alg.get_searchid()))) << "\t"
        << g.get_filename() << " " << parser.get_problemfile() << std::endl;
    }
}

void
run_fch_cpg(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
        std::string alg_name, std::string gr, std::string co)
{
    std::string orderfile = cfg.get_param_value("input");
    std::string gridmapfile = cfg.get_param_value("input");
    if(co == "" || gr == "" || orderfile == "" || gridmapfile == "")
    {
        std::cerr << "err; insufficient input params for alg " << alg_name 
                  << ". required (in, order) "
                  << "--input [ gr file ] [co file]"
                  << " [contraction order file] [gridmap file]\n";
        return;
    }
    
    // load up the node order
    std::vector<uint32_t> order;
    if(!warthog::ch::load_node_order(orderfile.c_str(), order, true))
    {
        std::cerr << "err; could not load contraction order file\n";
        return;
    }

    // load up the grid
    std::shared_ptr<warthog::gridmap> map(
            new warthog::gridmap(gridmapfile.c_str()));
            
    // load up the graph 
    std::shared_ptr<warthog::graph::planar_graph> pg(
            new warthog::graph::planar_graph());
    pg->load_dimacs(gr.c_str(), co.c_str(), false, true);

    // wrapper
    warthog::graph::corner_point_graph cpg(map, pg);

    // expander
    warthog::fch_cpg_expansion_policy fexp(&cpg, &order); 

    // heuristic 
    warthog::euclidean_heuristic h(pg.get());

    warthog::flexible_astar< warthog::euclidean_heuristic, 
        warthog::fch_cpg_expansion_policy> alg(&h, &fexp);
    alg.set_verbose(verbose);

    std::cerr << "running experiments\n";
    int i = 0;
    std::cout << "id\talg\texp\tgen\ttouch\tmicros\tplen\tmap\n";
    for(warthog::dimacs_parser::experiment_iterator it = parser.experiments_begin(); 
            it != parser.experiments_end(); it++)
    {
        warthog::dimacs_parser::experiment exp = (*it);

        uint32_t startid = map.get()->to_padded_id(exp.source);
        uint32_t targetid = map.get()->to_padded_id(exp.target);
		double len = alg.get_length(warthog::problem_instance(
                    startid, targetid));

        std::cout << i++ <<"\t" << alg_name << "\t" 
        << alg.get_nodes_expanded() << "\t" 
        << alg.get_nodes_generated() << "\t"
        << alg.get_nodes_touched() << "\t"
        << alg.get_search_time()  << "\t"
        << len << "\t" 
        << gr << " " << parser.get_problemfile() << std::endl;
    }
}

void
run_fch_af_cpg(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
        std::string alg_name, std::string gr, std::string co)
{
    std::string orderfile = cfg.get_param_value("input");
    std::string arclabels_file = cfg.get_param_value("input");
    std::string partition_file = cfg.get_param_value("input");
    std::string gridmapfile = cfg.get_param_value("input");
    if( orderfile == "" || arclabels_file == "" || 
        partition_file == "" || gridmapfile == "" )
    {
        std::cerr << "err; insufficient input parameters for --alg "
                  << alg_name << ". required, in order:\n"
                  << " --input [gr file] [co file] "
                  << " [contraction order file] [arclabels file] "
                  << " [graph partition file] [gridmap file]\n";
        return;
    }
    
    // load up the node order
    std::vector<uint32_t> order;
    if(!warthog::ch::load_node_order(orderfile.c_str(), order, true))
    {
        std::cerr << "err; could not load contraction order file\n";
        return;
    }

    // we insert two extra elements in the event that we
    // need to insert the start or target. both have the lowest
    // possible rank in the hierarchy (0 and 1)
    // NB: along the way we need to increase all ranks by 2 
    for(uint32_t i = 0; i < order.size(); i++) order.at(i)+=2;
    order.push_back(0);
    order.push_back(1);

    // load up the node partition info
    std::vector<uint32_t> part;
    if(!warthog::helpers::load_integer_labels_dimacs(
            partition_file.c_str(), part))
    {
        std::cerr << "err; could not load graph partition input file\n";
        return;
    }

    // load up the grid
    std::shared_ptr<warthog::gridmap> map(
            new warthog::gridmap(gridmapfile.c_str()));
            
    // load up the graph 
    std::shared_ptr<warthog::graph::planar_graph> pg(
            new warthog::graph::planar_graph());
    pg->load_dimacs(gr.c_str(), co.c_str(), false, true);

    // graph wrapper
    warthog::graph::corner_point_graph cpg(map, pg);

    // load up the arc-flags
    warthog::af_filter filter(pg.get(), &part);
    if(!filter.load_labels(arclabels_file.c_str()))
    {
        std::cerr << "err; could not load arcflags file\n";
        return;
    }

    // expander
    warthog::fch_af_cpg_expansion_policy fexp(&cpg, &order, &filter);

    // heuristic 
    warthog::euclidean_heuristic h(pg.get());
    h.set_hscale(warthog::ONE);

    warthog::flexible_astar< warthog::euclidean_heuristic, 
        warthog::fch_af_cpg_expansion_policy> alg(&h, &fexp);
    alg.set_verbose(verbose);

    std::cerr << "running experiments\n";
    int i = 0;
    std::cout << "id\talg\texp\tgen\ttouch\tmicros\tplen\tmap\n";
    for(warthog::dimacs_parser::experiment_iterator it = parser.experiments_begin(); 
            it != parser.experiments_end(); it++)
    {
        warthog::dimacs_parser::experiment exp = (*it);

        uint32_t startid = map.get()->to_padded_id(exp.source);
        uint32_t targetid = map.get()->to_padded_id(exp.target);
		double len = alg.get_length(warthog::problem_instance(
                    startid, targetid));

        std::cout << i++ <<"\t" << alg_name << "\t" 
        << alg.get_nodes_expanded() << "\t" 
        << alg.get_nodes_generated() << "\t"
        << alg.get_nodes_touched() << "\t"
        << alg.get_search_time()  << "\t"
        << len << "\t" 
        << gr << " " << parser.get_problemfile() << std::endl;
    }
}

void
run_fch_bb_cpg(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
        std::string alg_name, std::string gr, std::string co)
{
    std::string orderfile = cfg.get_param_value("input");
    std::string arclabels_file = cfg.get_param_value("input");
    std::string gridmapfile = cfg.get_param_value("input");
    if(orderfile == "" || arclabels_file == "" || gridmapfile == "")
    {
        std::cerr << "err; insufficient input parameters for --alg "
                  << alg_name << ". required, in order:\n"
                  << " --input [gr file] [co file] "
                  << " [contraction order file] [arclabels file]"
                  << " [gridmap file]\n";
        return;
    }

    // load up the node order
    std::vector<uint32_t> order;
    bool sort_by_id = true;
    if(!warthog::ch::load_node_order(orderfile.c_str(), order, sort_by_id))
    {
        std::cerr << "err; could not load node order input file\n";
        return;
    }

    // load up the grid
    std::shared_ptr<warthog::gridmap> map(
            new warthog::gridmap(gridmapfile.c_str()));

    // load up the graph 
    std::shared_ptr<warthog::graph::planar_graph> pg(
            new warthog::graph::planar_graph());
    pg->load_dimacs(gr.c_str(), co.c_str(), false, true);

    // wrapper
    warthog::graph::corner_point_graph cpg(map, pg);

    // load up the arc-flags
    warthog::bb_filter filter(pg.get());
    if(!filter.load_labels(arclabels_file.c_str()))
    {
        std::cerr << "err; could not load arcflags file\n";
        return;
    }

    // create a search algo
    warthog::euclidean_heuristic h(pg.get());
    warthog::fch_bb_cpg_expansion_policy fexp(&cpg, &order, &filter);
    warthog::flexible_astar< warthog::euclidean_heuristic, 
       warthog::fch_bb_cpg_expansion_policy>
    alg(&h, &fexp);
    alg.set_verbose(verbose);

    std::cerr << "running experiments\n";
    int i = 0;
    std::cout << "id\talg\texp\tgen\ttouch\tmicros\tplen\tmap\n";
    for(warthog::dimacs_parser::experiment_iterator it = parser.experiments_begin(); 
            it != parser.experiments_end(); it++)
    {
        warthog::dimacs_parser::experiment exp = (*it);
        uint32_t startid = map.get()->to_padded_id(exp.source);
        uint32_t targetid = map.get()->to_padded_id(exp.target);
		double len = alg.get_length(warthog::problem_instance(
                    startid, targetid));

        std::cout << i++ <<"\t" << alg_name << "\t" 
        << alg.get_nodes_expanded() << "\t" 
        << alg.get_nodes_generated() << "\t"
        << alg.get_nodes_touched() << "\t"
        << alg.get_search_time()  << "\t"
        << len << "\t" 
        << pg->get_filename() << " " << parser.get_problemfile() << std::endl;
    }
}

void
run_dimacs(warthog::util::cfg& cfg)
{
    std::string gr = cfg.get_param_value("input");
    std::string co = cfg.get_param_value("input");
    std::string problemfile = cfg.get_param_value("problem");
    std::string alg_name = cfg.get_param_value("alg");

    if((problemfile == ""))
    {
        std::cerr << "parameter is missing: --problem\n";
        return;
    }
    if((gr == "") || co == "")
    {
        std::cerr << "parameter is missing: --input [gr file] [co file]\n";
        return;
    }
    if((alg_name == ""))
    {
        std::cerr << "parameter is missing: --alg\n";
        return;
    }

    warthog::dimacs_parser parser;
    parser.load_instance(problemfile.c_str());
    if(parser.num_experiments() == 0)
    {
        std::cerr << "err; specified problem file contains no instances\n";
        return;
    }

    if(alg_name == "dijkstra")
    {
        run_dijkstra(cfg, parser, alg_name, gr, co);
    }
    else if(alg_name == "astar")
    {
        run_astar(cfg, parser, alg_name, gr, co);
    }
    else if(alg_name == "bi-dijkstra")
    {
        run_bi_dijkstra(cfg, parser, alg_name, gr, co);
    }
    else if(alg_name == "bi-astar")
    {
        run_bi_astar(cfg, parser, alg_name, gr, co);
    }
    else if(alg_name == "ch")
    {
        run_ch(cfg, parser, alg_name, gr, co);
    }
    else if(alg_name == "ch-astar")
    {
        run_ch_astar(cfg, parser, alg_name, gr, co);
    }
    else if(alg_name == "chase")
    {
        run_chase(cfg, parser, alg_name, gr, co);
    }
    else if(alg_name == "chaf")
    {
        run_chaf(cfg, parser, alg_name, gr, co);
    }
    else if(alg_name == "chaf-bb")
    {
        run_chaf_bb(cfg, parser, alg_name, gr, co);
    }
    else if(alg_name == "ch-cpg")
    {
        run_ch_cpg(cfg, parser, alg_name, gr, co);
    }
    else if(alg_name == "fch")
    {
        run_fch(cfg, parser, alg_name, gr, co);
    }
    else if(alg_name == "fch-dd")
    {
        run_fch_dd(cfg, parser, alg_name, gr, co);
    }
    else if(alg_name == "fch-af")
    {
        run_fch_af(cfg, parser, alg_name, gr, co);
    }
    else if(alg_name == "fch-bb")
    {
        run_fch_bb(cfg, parser, alg_name, gr, co);
    }
    else if(alg_name == "fch-bbaf")
    {
        run_fch_bbaf(cfg, parser, alg_name, gr, co);
    }
    else if(alg_name == "fch-dcl")
    {
        run_fch_dcl(cfg, parser, alg_name, gr, co);
    }
    else if(alg_name == "fch-cpg")
    {
        run_fch_cpg(cfg, parser, alg_name, gr, co);
    }
    else if(alg_name == "fch-af-cpg")
    {
        run_fch_af_cpg(cfg, parser, alg_name, gr, co);
    }
    else if(alg_name == "fch-bb-cpg")
    {
        run_fch_bb_cpg(cfg, parser, alg_name, gr, co);
    }
    else if(alg_name == "fch-bbaf-cpg")
    {
        run_fch_bbaf_cpg(cfg, parser, alg_name, gr, co);
    }
    else
    {
        std::cerr << "invalid search algorithm\n";
    }
}

int 
main(int argc, char** argv)
{
	// parse arguments
	warthog::util::param valid_args[] = 
	{
		{"alg",  required_argument, 0, 1},
		{"help", no_argument, &print_help, 1},
		{"checkopt",  no_argument, &checkopt, 1},
		{"verbose",  no_argument, &verbose, 1},
		{"input",  required_argument, 0, 1},
		{"problem",  required_argument, 0, 1},
	};

	warthog::util::cfg cfg;
	cfg.parse_args(argc, argv, "-f", valid_args);

    if(argc == 1 || print_help)
    {
		help();
        exit(0);
    }

    run_dimacs(cfg);
}

