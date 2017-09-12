// programs/roadhog.cpp 
//
// Pulls together a variety of different algorithms for 
// routing on road graphs.
//
// @author: dharabor
// @created: 2016-11-24
//

#include "af_filter.h"
#include "apex_filter.h"
#include "bbaf_filter.h"
#include "bb_filter.h"
#include "bidirectional_search.h"
#include "cfg.h"
#include "ch_expansion_policy.h"
#include "chaf_expansion_policy.h"
#include "chafbb_expansion_policy.h"
#include "chbb_expansion_policy.h"
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
#include "solution.h"
#include "timer.h"
#include "zero_heuristic.h"

// stuff for searching corner point graphs
#include "ch_cpg_expansion_policy.h"
#include "fch_af_cpg_expansion_policy.h"
#include "fch_bb_cpg_expansion_policy.h"
#include "fch_bbaf_cpg_expansion_policy.h"
#include "fch_cpg_expansion_policy.h"

#include "fch_jpg_expansion_policy.h"
#include "fch_bb_jpg_expansion_policy.h"
#include "fch_af_jpg_expansion_policy.h"

#include "getopt.h"

#include <fstream>
#include <functional>
#include <iomanip>
#include <memory>
#include <sstream>
#include <unordered_map>

// check computed solutions are optimal
int checkopt = 0;
// print debugging info during search
int verbose = 0;
// display program help on startup
int print_help = 0;

// suppress the header row when printing results? (default: no)
int suppress_header = 0;

uint32_t nruns = 1;

void
help()
{
	std::cerr << "Valid parameters:\n"
    << "\t--alg [ algorithm name (required) ]\n"
    << "\t--input [ algorithm-specific input files (omit to show options) ] \n"
    << "\t--problem [ ss or p2p problem file (required) ]\n"
	<< "\t--verbose (print debug info; omitting this param means no)\n"
	<< "\t--nruns [int (repeats per instance; default=" << nruns << ")]\n"
    << "\nRecognised values for --alg:\n"
    << "\tastar, dijkstra, bi-astar, bi-dijkstra\n"
    << "\tch, ch-astar, chase, chaf, chbb, chaf-bb, ch-cpg\n"
    << "\tfch, fchx, fch-dd, fch-af, fch-bb, fch-bbaf, fch-dcl\n"
    << "\tfch-cpg, fch-af-cpg, fch-bb-cpg, fch-bbaf-cpg\n"
    << "\tfch-jpg, fch-bb-jpg, fch-af-jpg\n"
    << "\nRecognised values for --input:\n "
    << "\ttoo many to list. missing input files will be listed at uuntime\n";
}

////////////////////////////////////////////////////////////////////////////
// ----- These two closures were useful once upon a time. Maybe again? --- 
////////////////////////////////////////////////////////////////////////////
//    std::function<uint32_t(warthog::search_node*)> fn_get_apex = 
//        [&order] (warthog::search_node* n) -> uint32_t
//        {
//            while(true)
//            {
//                warthog::search_node* p = n->get_parent();
//                if(!p || order.at(p->get_id()) < order.at(n->get_id()))
//                { break; }
//                n = p;
//            }
//            return order.at(n->get_id());
//        };
//
//    std::function<uint32_t(uint32_t)> fn_redundant_expansions = 
//        [&fexp, &order, &alg] (uint32_t apex) -> uint32_t
//        {
//            std::set<uint32_t> tmp;
//            for(uint32_t i = 0; i < order.size(); i++)
//            {
//                warthog::search_node* n = fexp.get_ptr(i, alg.get_search_id());
//                if(n && order.at(n->get_id()) > apex)
//                { 
//                    tmp.insert(n->get_id());
//                }
//            }
//
//            for(uint32_t j = 0; j < order.size(); j++)
//            {
//                warthog::search_node* n = fexp.get_ptr(j, alg.get_search_id());
//                warthog::search_node* m = n;
//                while(m)
//                {
//                    if(tmp.find(m->get_id()) != tmp.end())
//                    {
//                        while(tmp.find(n->get_id()) == tmp.end())
//                        {
//                            tmp.insert(n->get_id());
//                            n = n->get_parent();
//                        }
//                        break;
//                    }
//                    m = m->get_parent();
//                }
//            }
//            return tmp.size();
//        };

void
run_experiments( warthog::search* algo, std::string alg_name, 
        warthog::dimacs_parser& parser, std::ostream& out)
{
    std::cerr << "running experiments\n";
    std::cerr << "(averaging over " << nruns << " runs per instance)\n";

    if(!suppress_header)
    {
        std::cout 
            << "id\talg\texpanded\tinserted\tupdated\ttouched"
            << "\tmicros\tpcost\tplen\tmap\n";
    }
    uint32_t exp_id = 0;
    for(auto it = parser.experiments_begin(); 
            it != parser.experiments_end(); 
            it++)
    {
        warthog::dimacs_parser::experiment exp = (*it);
        warthog::solution sol;
        warthog::problem_instance pi(
                exp.source, (exp.p2p ? exp.target : warthog::INF), verbose);
        uint32_t expanded=0, inserted=0, updated=0, touched=0;
        double micros = 0;
        for(uint32_t i = 0; i < nruns; i++)
        {
            sol.reset();
            algo->get_path(pi, sol);

            expanded += sol.nodes_expanded_;
            inserted += sol.nodes_inserted_;
            touched += sol.nodes_touched_;
            updated += sol.nodes_updated_;
            micros += sol.time_elapsed_micro_;
        }

        out
            << exp_id++ <<"\t" 
            << alg_name << "\t" 
            << expanded / nruns << "\t" 
            << inserted / nruns << "\t"
            << updated / nruns << "\t"
            << touched / nruns << "\t"
            << micros / (double)nruns << "\t"
            << sol.sum_of_edge_costs_ << "\t" 
            << sol.path_.size() << "\t" 
            << parser.get_problemfile() 
            << std::endl;
    }
}

void
run_astar(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
        std::string alg_name, std::string gr, std::string co)
{
    warthog::graph::planar_graph g;
    if(!g.load_dimacs(gr.c_str(), co.c_str()))
    {
        std::cerr << "err; could not load gr or co input files " 
                  << "(one or both)\n";
        return;
    }

    warthog::graph_expansion_policy expander(&g);

    warthog::euclidean_heuristic h(&g);
    warthog::flexible_astar<
        warthog::euclidean_heuristic, 
        warthog::graph_expansion_policy> alg(&h, &expander);

    run_experiments(&alg, alg_name, parser, std::cout);
}

void
run_bi_astar(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
        std::string alg_name, std::string gr, std::string co)
{
    warthog::graph::planar_graph g;
    if(!g.load_dimacs(gr.c_str(), co.c_str(), false, true))
    {
        std::cerr << "err; could not load gr or co input files " 
                  << "(one or both)\n";
        return;
    }

    warthog::graph_expansion_policy fexp(&g);

    warthog::graph::planar_graph backward_g;
    if(!backward_g.load_dimacs(gr.c_str(), co.c_str(), true, true))
    {
        std::cerr << "err; could not load gr or co input files " 
                  << "(one or both)\n";
        return;
    }
    warthog::graph_expansion_policy bexp(&g);

    warthog::euclidean_heuristic h(&g);
    warthog::bidirectional_search<warthog::euclidean_heuristic>
        alg(&fexp, &bexp, &h);

    run_experiments(&alg, alg_name, parser, std::cout);
}

void
run_bi_dijkstra(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
        std::string alg_name, std::string gr, std::string co)
{
    warthog::graph::planar_graph g;
    if(!g.load_dimacs(gr.c_str(), co.c_str()))
    {
        std::cerr << "err; could not load gr or co input files " 
                  << "(one or both)\n";
        return;
    }
    warthog::graph_expansion_policy fexp(&g);

    warthog::graph::planar_graph backward_g;
    if(!backward_g.load_dimacs(gr.c_str(), co.c_str(), true, true))
    {
        std::cerr << "err; could not load gr or co input files " 
                  << "(one or both)\n";
        return;
    }
    warthog::graph_expansion_policy bexp(&backward_g);

    warthog::zero_heuristic h;
    warthog::bidirectional_search<warthog::zero_heuristic>
        alg(&fexp, &bexp, &h);

    run_experiments(&alg, alg_name, parser, std::cout);
}

void
run_dijkstra(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
        std::string alg_name, std::string gr, std::string co)
{
    warthog::graph::planar_graph g;

    if(!g.load_dimacs(gr.c_str(), co.c_str()))
    {
        std::cerr << "err; could not load gr or co input files " 
                  << "(one or both)\n";
        return;
    }

    warthog::graph_expansion_policy expander(&g);

    warthog::zero_heuristic h;
    warthog::flexible_astar<warthog::zero_heuristic, 
        warthog::graph_expansion_policy> alg(&h, &expander);

    run_experiments(&alg, alg_name, parser, std::cout);
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
    if(!g.load_dimacs(gr.c_str(), co.c_str(), false, true))
    {
        std::cerr << "err; could not load gr or co input files " 
                  << "(one or both)\n";
        return;
    }

    // load up the node order
    std::vector<uint32_t> order;
    if(!warthog::ch::load_node_order(orderfile.c_str(), order, true))
    {
        std::cerr << "err; could not load contraction order file\n";
        return;
    }

    std::cerr << "preparing to search\n";
    warthog::ch_expansion_policy fexp(&g, &order);
    warthog::ch_expansion_policy bexp (&g, &order, true);
    warthog::zero_heuristic h;
    warthog::bidirectional_search<warthog::zero_heuristic> alg(&fexp, &bexp, &h);

    run_experiments(&alg, alg_name, parser, std::cout);
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
    if(!g.load_dimacs(gr.c_str(), co.c_str(), false, true))
    {
        std::cerr 
            << "err; could not load gr or co input files (one or both)\n";
        return;
    }

    // load up the node order
    std::vector<uint32_t> order;
    if(!warthog::ch::load_node_order(orderfile.c_str(), order, true))
    {
        std::cerr << "err; could not load node order input file\n";
        return;
    }

    std::cerr << "preparing to search\n";
    warthog::euclidean_heuristic h(&g);
    warthog::ch_expansion_policy fexp(&g, &order);
    warthog::ch_expansion_policy bexp (&g, &order, true);
    warthog::bidirectional_search<warthog::euclidean_heuristic> alg(&fexp, &bexp, &h);

    run_experiments(&alg, alg_name, parser, std::cout);
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
    if(!pg->load_dimacs(gr.c_str(), co.c_str(), false, true))
    {
        std::cerr << "err; could not load gr or co input files " 
                  << "(one or both)\n";
        return;
    }

    // wrapper
    warthog::graph::corner_point_graph cpg(map, pg);

    // expander
	warthog::zero_heuristic h;
	warthog::ch::ch_cpg_expansion_policy fexp(&cpg, &order);
	warthog::ch::ch_cpg_expansion_policy bexp(&cpg, &order, true);
    warthog::bidirectional_search<warthog::zero_heuristic> alg(&fexp, &bexp, &h);

    run_experiments(&alg, alg_name, parser, std::cout);
}

//void
//run_chase(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
//        std::string alg_name, std::string gr, std::string co)
//{
//    std::string orderfile = cfg.get_param_value("input");
//
//    // load up the graph 
//    warthog::graph::planar_graph g;
//    if(!g.load_dimacs(gr.c_str(), co.c_str(), false, true))
//    {
//        std::cerr << "err; could not load gr or co input files " 
//                  << "(one or both)\n";
//        return;
//    }
//
//    // load up the node order
//    std::vector<uint32_t> order;
//    if(!warthog::ch::load_node_order(orderfile.c_str(), order, true))
//    {
//        std::cerr << "err; could not load contraction order file\n";
//        return;
//    }
//
//    std::cerr << "preparing to search\n";
//    warthog::zero_heuristic h;
//    warthog::ch_expansion_policy fexp(&g, &order);
//    warthog::ch_expansion_policy bexp (&g, &order, true);
//    warthog::chase_search<warthog::zero_heuristic> alg(&fexp, &bexp, &h);
//
//    run_experiments(&alg, alg_name, parser, std::cout);
//}

void
run_chase(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
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
    std::shared_ptr<warthog::label::af_labelling> afl(
        warthog::label::af_labelling::load(
            arclabels_file.c_str(), &g, &part));
    if(!afl.get())
    {
        std::cerr << "err; could not load arcflags file\n";
        return;
    }
    warthog::af_filter filter(afl.get());

    std::cerr << "preparing to search\n";
    warthog::zero_heuristic h;
    warthog::chaf_expansion_policy fexp(&g, &order, &filter);
    warthog::chaf_expansion_policy bexp (&g, &order, &filter, true);
    warthog::chase_search<
        warthog::chaf_expansion_policy, warthog::zero_heuristic> 
        alg(&fexp, &bexp, &h);

    run_experiments(&alg, alg_name, parser, std::cout);
}

void
run_chbb(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
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

    // load up the arc labels
    std::shared_ptr<warthog::label::bb_labelling> bbl
        (warthog::label::bb_labelling::load(arclabels_file.c_str(), &g));
    if(!bbl.get())
    {
        std::cerr << "err; could not load arcflags file\n";
        return;
    }
    warthog::bb_filter filter(bbl.get());

    std::cerr << "preparing to search\n";
    warthog::chbb_expansion_policy fexp(&g, &order, &filter);
    warthog::chbb_expansion_policy bexp (&g, &order, &filter, true);
    warthog::zero_heuristic h;
    warthog::bidirectional_search<warthog::zero_heuristic> 
        alg(&fexp, &bexp, &h);

    //warthog::euclidean_heuristic h(&g);
    //warthog::bidirectional_search<warthog::euclidean_heuristic> 
    //   alg(&fexp, &bexp, &h);

    run_experiments(&alg, alg_name, parser, std::cout);
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
    std::shared_ptr<warthog::label::af_labelling> afl(
        warthog::label::af_labelling::load(
            arclabels_file.c_str(), &g, &part));
    if(!afl.get())
    {
        std::cerr << "err; could not load arcflags file\n";
        return;
    }
    warthog::af_filter filter(afl.get());

    std::cerr << "preparing to search\n";
    warthog::chaf_expansion_policy fexp(&g, &order, &filter);
    warthog::chaf_expansion_policy bexp (&g, &order, &filter, true);
    warthog::zero_heuristic h;
    warthog::bidirectional_search<warthog::zero_heuristic> 
        alg(&fexp, &bexp, &h);

    //warthog::euclidean_heuristic h(&g);
    //warthog::bidirectional_search<warthog::euclidean_heuristic> 
    //   alg(&fexp, &bexp, &h);

    run_experiments(&alg, alg_name, parser, std::cout);
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
    std::shared_ptr<warthog::label::bbaf_labelling> lab(
            warthog::label::bbaf_labelling::load(
                    arclabels_file.c_str(), &g, &part));
    if(!lab.get())
    {
        std::cerr << "err; could not load arcflags file\n";
        return;
    }

    std::cerr << "preparing to search\n";
    warthog::zero_heuristic h;
    //warthog::euclidean_heuristic h(&g);
    warthog::chafbb_expansion_policy fexp(&g, &order, lab.get());
    warthog::chafbb_expansion_policy bexp (&g, &order, lab.get(), true);
    warthog::bidirectional_search<warthog::zero_heuristic> 
        alg(&fexp, &bexp, &h);
    //warthog::bidirectional_search<warthog::euclidean_heuristic> alg(&fexp, &bexp, &h);

    run_experiments(&alg, alg_name, parser, std::cout);
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
    if(!g.load_dimacs(gr.c_str(), co.c_str(), false, true))
    {
        std::cerr 
            << "err; could not load gr or co input files (one or both)\n";
        return;
    }

    // load up the node order
    std::vector<uint32_t> order;
    if(!warthog::ch::load_node_order(orderfile.c_str(), order, true))
    {
        std::cerr << "err; could not load node order input file\n";
        return;
    }

    std::cerr << "preparing to search\n";
    warthog::fch_expansion_policy fexp(&g, &order); 
    warthog::euclidean_heuristic h(&g);

    warthog::flexible_astar< warthog::euclidean_heuristic, 
        warthog::fch_expansion_policy> alg(&h, &fexp);
    
    // extra metric; how many nodes do we expand above the apex?
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

    run_experiments(&alg, alg_name, parser, std::cout);
}

void
run_fch_apex_experiment(warthog::util::cfg& cfg, 
        warthog::dimacs_parser& parser, std::string alg_name, 
        std::string gr, std::string co)
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
    warthog::apex_filter filter(&order);

    std::string apex_filter_type = cfg.get_param_value("alg");
    if(apex_filter_type == "no_paa")
    {
        // don't prune up nodes after the apex is reached
        // (weaker filter)
        filter.prune_after_apex_ = false;
        alg_name = alg_name + "_" + apex_filter_type;
    }
    if(apex_filter_type == "no_pba")
    {
        // don't prune down nodes before the apex is reached
        // (weaker filter)
        filter.prune_before_apex_ = false;
        alg_name = alg_name + "_" + apex_filter_type;
    }

    // reference algo
    warthog::flexible_astar<warthog::euclidean_heuristic, 
        warthog::fch_expansion_policy>
            fch(&h, &fexp);

    // fchx
    warthog::flexible_astar<warthog::euclidean_heuristic, 
        warthog::fch_expansion_policy, warthog::apex_filter> 
            fchx(&h, &fexp, &filter);

    std::cerr << "running fch apex experiment\n";
    if(!suppress_header)
    {
        std::cout 
            << "id\talg\texpanded\tinserted\tupdated\ttouched"
            << "\tmicros\tpcost\tplen\tmap\n";
    }
    uint32_t exp_id = 0;
    for(auto it = parser.experiments_begin(); 
            it != parser.experiments_end(); 
            it++)
    {
        warthog::dimacs_parser::experiment exp = (*it);
		warthog::solution sol;
        warthog::problem_instance pi(
                exp.source, (exp.p2p ? exp.target : warthog::INF), verbose);
        // run once
        fch.get_path(pi, sol);

        // identify the apex of the path
        if(sol.path_.size() > 0)
        {
            filter.set_apex(sol.path_.at(0));
            for(uint32_t i = 1; i < sol.path_.size(); i++)
            {
                if(order.at(filter.get_apex()) < order.at(sol.path_.at(i)))
                {
                    filter.set_apex(sol.path_.at(i));
                }
            }
        }

        // run the experiment again, using the apex to prune
        sol.reset();
        fchx.get_path(pi, sol);

        std::cout
            << exp_id++ <<"\t" 
            << alg_name << "\t" 
            << sol.nodes_expanded_ << "\t" 
            << sol.nodes_inserted_ << "\t"
            << sol.nodes_updated_ << "\t"
            << sol.nodes_touched_ << "\t"
            << sol.time_elapsed_micro_ << "\t"
            << sol.sum_of_edge_costs_ << "\t" 
            << sol.path_.size() << "\t" 
            << parser.get_problemfile() 
            << std::endl;
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

    run_experiments(&alg, alg_name, parser, std::cout);
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
    std::shared_ptr<warthog::label::af_labelling> afl
        (warthog::label::af_labelling::load
        (arclabels_file.c_str(), &g, &part));
    if(!afl.get())
    {
        std::cerr << "err; could not load arcflags file\n";
        return;
    }

    warthog::af_filter filter(afl.get());
    warthog::euclidean_heuristic h(&g);
    warthog::fch_af_expansion_policy fexp(&g, &order, &filter);

    warthog::flexible_astar< warthog::euclidean_heuristic, 
       warthog::fch_af_expansion_policy>
    alg(&h, &fexp);

    run_experiments(&alg, alg_name, parser, std::cout);
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
        std::cerr 
            << "err; could not load gr or co input files (one or both)\n";
        return;
    }
    // load up the arc labels
    std::shared_ptr<warthog::label::bb_labelling> bbl
        (warthog::label::bb_labelling::load(arclabels_file.c_str(), &g));
    if(!bbl.get())
    {
        std::cerr << "err; could not load arcflags file\n";
        return;
    }
    warthog::bb_filter filter(bbl.get());

    warthog::euclidean_heuristic h(&g);
    warthog::fch_bb_expansion_policy fexp(&g, &order, &filter);

    warthog::flexible_astar< warthog::euclidean_heuristic, 
       warthog::fch_bb_expansion_policy>
    alg(&h, &fexp);

    run_experiments(&alg, alg_name, parser, std::cout);

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
    std::shared_ptr<warthog::label::bbaf_labelling> lab(
            warthog::label::bbaf_labelling::load(
                arclabels_file.c_str(), &g, &part));

    if(!lab.get())
    {
        std::cerr << "err; could not load arcflags file\n";
        return;
    }

    std::cerr << "preparing to search\n";
    warthog::euclidean_heuristic h(&g);
    warthog::fch_bbaf_expansion_policy fexp(&g, &order, lab.get());

    warthog::flexible_astar< warthog::euclidean_heuristic, 
       warthog::fch_bbaf_expansion_policy>
    alg(&h, &fexp);

    run_experiments(&alg, alg_name, parser, std::cout);
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
    std::shared_ptr<warthog::label::bbaf_labelling> lab(
            warthog::label::bbaf_labelling::load(
                arclabels_file.c_str(), pg.get(), &part));
    if(!lab.get())
    {
        std::cerr << "err; could not load arcflags file\n";
        return;
    }

    // search algo
    warthog::euclidean_heuristic h(pg.get());
    warthog::fch_bbaf_cpg_expansion_policy fexp(&cpg, &order, lab.get());
    warthog::flexible_astar< warthog::euclidean_heuristic, 
       warthog::fch_bbaf_cpg_expansion_policy>
    alg(&h, &fexp);

    run_experiments(&alg, alg_name, parser, std::cout);
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

    run_experiments(&alg, alg_name, parser, std::cout);
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

    run_experiments(&alg, alg_name, parser, std::cout);
}

void
run_fch_jpg(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
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
    warthog::fch_jpg_expansion_policy fexp(&cpg, &order); 

    // heuristic 
    warthog::euclidean_heuristic h(pg.get());

    warthog::flexible_astar< warthog::euclidean_heuristic, 
        warthog::fch_jpg_expansion_policy> alg(&h, &fexp);

    run_experiments(&alg, alg_name, parser, std::cout);
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
    std::shared_ptr<warthog::label::af_labelling> afl(
        warthog::label::af_labelling::load(
            arclabels_file.c_str(), pg.get(), &part));
    if(!afl.get())
    {
        std::cerr << "err; could not load arcflags file\n";
        return;
    }
    warthog::af_filter filter(afl.get());

    // search algo
    warthog::euclidean_heuristic h(pg.get());
    warthog::fch_af_cpg_expansion_policy fexp(&cpg, &order, afl.get());
    warthog::flexible_astar< warthog::euclidean_heuristic, 
        warthog::fch_af_cpg_expansion_policy> alg(&h, &fexp);

    run_experiments(&alg, alg_name, parser, std::cout);
}

void
run_fch_af_jpg(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
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
    std::shared_ptr<warthog::label::af_labelling> afl(
            warthog::label::af_labelling::load(
                arclabels_file.c_str(), pg.get(), &part));
    if(!afl.get())
    {
        std::cerr << "err; could not load arcflags file\n";
        return;
    }

    // search algo
    warthog::euclidean_heuristic h(pg.get());
    warthog::fch_af_jpg_expansion_policy fexp(&cpg, &order, afl.get());
    warthog::flexible_astar< warthog::euclidean_heuristic, 
        warthog::fch_af_jpg_expansion_policy> alg(&h, &fexp);

    run_experiments(&alg, alg_name, parser, std::cout);
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
    std::shared_ptr<warthog::label::bb_labelling> bbl(
            warthog::label::bb_labelling::load(
                arclabels_file.c_str(), pg.get()));
    if(!bbl.get())
    {
        std::cerr << "err; could not load arcflags file\n";
        return;
    }

    // create a search algo
    warthog::euclidean_heuristic h(pg.get());
    warthog::fch_bb_cpg_expansion_policy fexp(&cpg, &order, bbl.get());
    warthog::flexible_astar< warthog::euclidean_heuristic, 
       warthog::fch_bb_cpg_expansion_policy>
    alg(&h, &fexp);

    run_experiments(&alg, alg_name, parser, std::cout);
}

void
run_fch_bb_jpg(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
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
    warthog::graph::corner_point_graph jpg(map, pg);

    // load up the arc-flags
    std::shared_ptr<warthog::label::bb_labelling> bbl
        (warthog::label::bb_labelling::load(
            arclabels_file.c_str(), pg.get()));
    if(!bbl.get())
    {
        std::cerr << "err; could not load arcflags file\n";
        return;
    }

    // create a search algo
    warthog::euclidean_heuristic h(pg.get());
    warthog::fch_bb_jpg_expansion_policy fexp(&jpg, &order, bbl.get());
    warthog::flexible_astar< warthog::euclidean_heuristic, 
       warthog::fch_bb_jpg_expansion_policy>
    alg(&h, &fexp);

    run_experiments(&alg, alg_name, parser, std::cout);
}

void
run_dimacs(warthog::util::cfg& cfg)
{
    std::string gr = cfg.get_param_value("input");
    std::string co = cfg.get_param_value("input");
    std::string problemfile = cfg.get_param_value("problem");
    std::string alg_name = cfg.get_param_value("alg");
    std::string par_nruns = cfg.get_param_value("nruns");

    if(par_nruns != "")
    {
       char* end;
       nruns = strtol(par_nruns.c_str(), &end, 10);
    }


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
    else if(alg_name == "chase")
    {
        run_chase(cfg, parser, alg_name, gr, co);
    }
    else if(alg_name == "ch-astar")
    {
        run_ch_astar(cfg, parser, alg_name, gr, co);
    }
    else if(alg_name == "chbb")
    {
        run_chbb(cfg, parser, alg_name, gr, co);
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
    else if(alg_name == "fchx")
    {
        run_fch_apex_experiment(cfg, parser, alg_name, gr, co);
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
    else if(alg_name == "fch-jpg")
    {
        run_fch_jpg(cfg, parser, alg_name, gr, co);
    }
    else if(alg_name == "fch-bb-jpg")
    {
        run_fch_bb_jpg(cfg, parser, alg_name, gr, co);
    }
    else if(alg_name == "fch-af-jpg")
    {
        run_fch_af_jpg(cfg, parser, alg_name, gr, co);
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
		{"nruns",  required_argument, 0, 1},
		{"help", no_argument, &print_help, 1},
		{"checkopt",  no_argument, &checkopt, 1},
		{"verbose",  no_argument, &verbose, 1},
		{"noheader",  no_argument, &suppress_header, 1},
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



