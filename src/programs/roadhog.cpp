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
#include "bch_search.h"
#include "bch_expansion_policy.h"
#include "bch_bb_expansion_policy.h"
#include "bidirectional_graph_expansion_policy.h"
#include "bidirectional_search.h"
#include "cfg.h"
#include "constants.h"
#include "contraction.h"
#include "cpd_heuristic.h"
#include "depth_first_search.h"
#include "dimacs_parser.h"
#include "euclidean_heuristic.h"
#include "fch_bb_expansion_policy.h"
#include "fch_dfs_expansion_policy.h"
#include "fch_expansion_policy.h"
#include "firstmove_labelling.h"
#include "fixed_graph_contraction.h"
#include "flexible_astar.h"
#include "graph_expansion_policy.h"
#include "graph_oracle.h"
#include "graph_oracle_expansion_policy.h"
#include "lazy_graph_contraction.h"
#include "xy_graph.h"
#include "solution.h"
#include "timer.h"
#include "zero_heuristic.h"

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

long nruns = 1;

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
    << "\tbch, bch-astar, bch-af, bch-bb, bch-bbaf\n"
    << "\tfch, fch-af, fch-bb, fch-bbaf, fch-dfs\n"
    << "\tdfs, cpd, cpd-search\n"
    << "\nRecognised values for --input:\n "
    << "\ttoo many to list. missing input files will be listed at runtime\n";
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
            << "\tnanos\tpcost\tplen\tmap\n";
    }
    uint32_t exp_id = 0;
    for(auto it = parser.experiments_begin(); 
            it != parser.experiments_end(); 
            it++)
    {
        warthog::dimacs_parser::experiment exp = (*it);
        warthog::solution sol;
        uint32_t start_id = exp.source;
        uint32_t target_id = exp.p2p ? exp.target : warthog::INF32;
        warthog::problem_instance pi(start_id, target_id, verbose);
        uint32_t expanded=0, inserted=0, updated=0, touched=0;
        double nano_time = DBL_MAX;
        for(uint32_t i = 0; i < nruns; i++)
        {
            sol.reset();
            algo->get_distance(pi, sol);

            expanded += sol.nodes_expanded_;
            inserted += sol.nodes_inserted_;
            touched += sol.nodes_touched_;
            updated += sol.nodes_updated_;
            nano_time = nano_time < sol.time_elapsed_nano_ 
                            ?  nano_time : sol.time_elapsed_nano_;
        }

        out
            << exp_id++ <<"\t" 
            << alg_name << "\t" 
            << expanded / nruns << "\t" 
            << inserted / nruns << "\t"
            << updated / nruns << "\t"
            << touched / nruns << "\t"
            << nano_time << "\t" /// (double)nruns << "\t"
            << sol.sum_of_edge_costs_ << "\t" 
            << (int32_t)((sol.path_.size() == 0) ? -1 : (int32_t)(sol.path_.size()-1)) << "\t" 
            << parser.get_problemfile() 
            << std::endl;
    }
}

void
run_experiments( std::function<void(warthog::problem_instance*, warthog::solution*)>& algo_fn, 
        std::string alg_name, warthog::dimacs_parser& parser, std::ostream& out)
{
    std::cerr << "running experiments\n";
    std::cerr << "(averaging over " << nruns << " runs per instance)\n";

    if(!suppress_header)
    {
        std::cout 
            << "id\talg\texpanded\tinserted\tupdated\ttouched"
            << "\tnanos\tpcost\tplen\tmap\n";
    }
    uint32_t exp_id = 0;
    for(auto it = parser.experiments_begin(); 
            it != parser.experiments_end(); 
            it++)
    {
        warthog::dimacs_parser::experiment exp = (*it);
        warthog::solution sol;
        uint32_t start_id = exp.source;
        uint32_t target_id = exp.p2p ? exp.target : warthog::INF32;
        warthog::problem_instance pi(start_id, target_id, verbose);
        uint32_t expanded=0, inserted=0, updated=0, touched=0;
        double nano_time = DBL_MAX;
        for(uint32_t i = 0; i < nruns; i++)
        {
            sol.reset();
            algo_fn(&pi, &sol);

            expanded += sol.nodes_expanded_;
            inserted += sol.nodes_inserted_;
            touched += sol.nodes_touched_;
            updated += sol.nodes_updated_;
            nano_time = nano_time < sol.time_elapsed_nano_ 
                            ?  nano_time : sol.time_elapsed_nano_;
        }

        out
            << exp_id++ <<"\t" 
            << alg_name << "\t" 
            << expanded / nruns << "\t" 
            << inserted / nruns << "\t"
            << updated / nruns << "\t"
            << touched / nruns << "\t"
            << nano_time << "\t" /// (double)nruns << "\t"
            << sol.sum_of_edge_costs_ << "\t" 
            << (int32_t)((sol.path_.size() == 0) ? -1 : (int32_t)(sol.path_.size()-1)) << "\t" 
            << parser.get_problemfile() 
            << std::endl;
    }
}

void
run_astar(warthog::util::cfg& cfg, 
    warthog::dimacs_parser& parser, std::string alg_name)
{
    std::string xy_filename = cfg.get_param_value("input");
    if(xy_filename == "")
    {
        std::cerr << "parameter is missing: --input [xy-graph file]\n";
        return;
    }

    warthog::graph::xy_graph g;
    std::ifstream ifs(xy_filename);
    warthog::graph::read_xy(ifs, g);

    warthog::simple_graph_expansion_policy expander(&g);
    warthog::euclidean_heuristic h(&g);
    warthog::pqueue_min open;

    warthog::flexible_astar<
        warthog::euclidean_heuristic, 
        warthog::simple_graph_expansion_policy, 
        warthog::pqueue_min> 
            alg(&h, &expander, &open);

    run_experiments(&alg, alg_name, parser, std::cout);
}

void
run_dijkstra(warthog::util::cfg& cfg, 
    warthog::dimacs_parser& parser, std::string alg_name )
{
    // load up the graph
    std::string xy_filename = cfg.get_param_value("input");
    if(xy_filename == "")
    {
        std::cerr << "err; require --input [xy-graph file]\n";
        return;
    }
    warthog::graph::xy_graph g;
    std::ifstream ifs(xy_filename);
    warthog::graph::read_xy(ifs, g);

    warthog::simple_graph_expansion_policy expander(&g);
    warthog::zero_heuristic h;
    warthog::pqueue<warthog::cmp_less_search_node_f_only, warthog::min_q> open;

    warthog::flexible_astar<
        warthog::zero_heuristic, 
        warthog::simple_graph_expansion_policy, 
        warthog::pqueue<warthog::cmp_less_search_node_f_only, warthog::min_q>> 
            alg(&h, &expander, &open);

    run_experiments(&alg, alg_name, parser, std::cout);
}

void
run_bi_astar( warthog::util::cfg& cfg, 
        warthog::dimacs_parser& parser, std::string alg_name)
{
    std::string xy_filename = cfg.get_param_value("input");
    if(xy_filename == "")
    {
        std::cerr << "parameter is missing: --input [xy-graph file]\n";
        return;
    }

    warthog::graph::xy_graph g;
    std::ifstream ifs(xy_filename);
    warthog::graph::read_xy(ifs, g, true);

    warthog::bidirectional_graph_expansion_policy fexp(&g, false);
    warthog::bidirectional_graph_expansion_policy bexp(&g, true);

    warthog::euclidean_heuristic h(&g);
    warthog::bidirectional_search<
        warthog::euclidean_heuristic,
        warthog::bidirectional_graph_expansion_policy> 
            alg(&fexp, &bexp, &h);

    run_experiments(&alg, alg_name, parser, std::cout);
}

void
run_bi_dijkstra( warthog::util::cfg& cfg,
     warthog::dimacs_parser& parser, std::string alg_name )
{
    std::string xy_filename = cfg.get_param_value("input");
    if(xy_filename == "")
    {
        std::cerr << "parameter is missing: --input [xy-graph file]\n";
        return;
    }

    std::ifstream ifs(xy_filename);
    warthog::graph::xy_graph g;
    warthog::graph::read_xy(ifs, g, true);
    ifs.close();

    warthog::bidirectional_graph_expansion_policy fexp(&g, false);
    warthog::bidirectional_graph_expansion_policy bexp(&g, true);

    warthog::zero_heuristic h;
    warthog::bidirectional_search<
        warthog::zero_heuristic,
        warthog::bidirectional_graph_expansion_policy>
            alg(&fexp, &bexp, &h);

    run_experiments(&alg, alg_name, parser, std::cout);
}


void
run_bch(warthog::util::cfg& cfg, 
        warthog::dimacs_parser& parser, std::string alg_name)
{
    std::string chd_file = cfg.get_param_value("input");
    if(chd_file == "")
    {
        std::cerr << "err; missing chd input file\n";
        return;
    }

    std::shared_ptr<warthog::ch::ch_data> 
        chd(warthog::ch::load_ch_data(chd_file.c_str(), warthog::ch::UP_ONLY));

    std::cerr << "preparing to search\n";
    warthog::bch_expansion_policy fexp(chd.get()->g_);
    warthog::bch_expansion_policy bexp (chd.get()->g_, true);
    warthog::zero_heuristic h;
    warthog::bch_search<
        warthog::zero_heuristic, 
        warthog::bch_expansion_policy> 
            alg(&fexp, &bexp, &h);

    run_experiments(&alg, alg_name, parser, std::cout);
}

void
run_bch_backwards_only(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
        std::string alg_name)
{
    std::string chd_file = cfg.get_param_value("input");
    if(chd_file == "")
    {
        std::cerr << "err; missing chd input file\n";
        return;
    }

    std::shared_ptr<warthog::ch::ch_data> 
        chd(warthog::ch::load_ch_data(chd_file.c_str(), warthog::ch::UP_ONLY));

    std::cerr << "preparing to search\n";
    warthog::bch_expansion_policy bexp (chd.get()->g_, true);
    warthog::zero_heuristic h;
    warthog::pqueue_min open;

    warthog::flexible_astar<
        warthog::zero_heuristic, 
        warthog::bch_expansion_policy,
        warthog::pqueue_min> 
            alg(&h, &bexp, &open);

    std::cerr << "running experiments\n";
    std::cerr << "(averaging over " << nruns << " runs per instance)\n";

    if(!suppress_header)
    {
        std::cout 
            << "id\talg\texpanded\tinserted\tupdated\ttouched"
            << "\tnanos\tpcost\tplen\tmap\n";
    }
    uint32_t exp_id = 0;
    for(auto it = parser.experiments_begin(); 
            it != parser.experiments_end(); 
            it++)
    {
        warthog::dimacs_parser::experiment exp = (*it);
        warthog::solution sol;
        warthog::problem_instance pi(exp.source, warthog::INF32, verbose);
        uint32_t expanded=0, inserted=0, updated=0, touched=0;
        double nano_time = DBL_MAX;
        for(uint32_t i = 0; i < nruns; i++)
        {
            sol.reset();
            alg.get_path(pi, sol);

            expanded += sol.nodes_expanded_;
            inserted += sol.nodes_inserted_;
            touched += sol.nodes_touched_;
            updated += sol.nodes_updated_;
            nano_time = nano_time < sol.time_elapsed_nano_ 
                            ?  nano_time : sol.time_elapsed_nano_;
        }

        std::cout
            << exp_id++ <<"\t" 
            << alg_name << "\t" 
            << expanded / nruns << "\t" 
            << inserted / nruns << "\t"
            << updated / nruns << "\t"
            << touched / nruns << "\t"
            << nano_time << "\t" /// (double)nruns << "\t"
            << sol.sum_of_edge_costs_ << "\t" 
            << sol.path_.size() << "\t" 
            << parser.get_problemfile() 
            << std::endl;
    }
}

void
run_bch_astar(warthog::util::cfg& cfg, 
              warthog::dimacs_parser& parser, std::string alg_name)
{
    std::string chd_file = cfg.get_param_value("input");
    if(chd_file == "")
    {
        std::cerr << "err; required --input [chd file]\n";
        return;
    }

    warthog::ch::ch_data& chd = 
        *(warthog::ch::load_ch_data(chd_file.c_str(), warthog::ch::UP_ONLY));

    std::cerr << "preparing to search\n";
    warthog::euclidean_heuristic h(chd.g_);
    warthog::bch_expansion_policy fexp(chd.g_);
    warthog::bch_expansion_policy bexp (chd.g_, true);
    warthog::bch_search<
        warthog::euclidean_heuristic,
        warthog::bch_expansion_policy>
            alg(&fexp, &bexp, &h);

    run_experiments(&alg, alg_name, parser, std::cout);
}

void
run_bch_bb(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
        std::string alg_name)
{
    // load up the contraction hierarchy
    std::string chd_file = cfg.get_param_value("input");
    std::string arclabels_file = cfg.get_param_value("input");
    if(chd_file == "" || arclabels_file == "")
    {
        std::cerr << "err; require --input [chd file] [arclabels file]\n";
        return;
    }
    warthog::ch::ch_data& chd = 
        *(warthog::ch::load_ch_data(chd_file.c_str(), warthog::ch::UP_ONLY));

    // load up the edge labels
    warthog::label::bb_labelling *fwd_lab_ptr=0, *bwd_lab_ptr=0;
    bool result = warthog::label::bb_labelling::load_bch_labels(
            arclabels_file.c_str(), chd.g_, chd.level_,
            fwd_lab_ptr, bwd_lab_ptr);
    if(!result)
    {
        std::cerr << "err; could not load arc labels\n";
        return;
    }

    std::shared_ptr<warthog::label::bb_labelling> fwd_lab(fwd_lab_ptr);
    std::shared_ptr<warthog::label::bb_labelling> bwd_lab(bwd_lab_ptr);

    warthog::bb_filter fwd_filter(fwd_lab.get());
    warthog::bb_filter bwd_filter(bwd_lab.get());

    std::cerr << "preparing to search\n";
    warthog::bch_bb_expansion_policy fexp(chd.g_, &fwd_filter);
    warthog::bch_bb_expansion_policy bexp (chd.g_, &bwd_filter, true);
    warthog::zero_heuristic h;
    warthog::bch_search<
        warthog::zero_heuristic,
        warthog::bch_bb_expansion_policy> 
            alg(&fexp, &bexp, &h);

    run_experiments(&alg, alg_name, parser, std::cout);
}

void
run_fch(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
        std::string alg_name)
{
    std::string chd_file = cfg.get_param_value("input");
    if(chd_file == "")
    {
        std::cerr << "err; missing chd input file\n";
        return;
    }

    // load up the graph 
    warthog::ch::ch_data& chd = 
        *(warthog::ch::load_ch_data(
            chd_file.c_str(), warthog::ch::UP_DOWN));

    std::cerr << "preparing to search\n";
    warthog::fch_expansion_policy fexp(chd.g_, chd.level_); 
    warthog::euclidean_heuristic h(chd.g_);
    warthog::pqueue_min open;

    warthog::flexible_astar< 
        warthog::euclidean_heuristic, 
        warthog::fch_expansion_policy,
        warthog::pqueue_min> 
            alg(&h, &fexp, &open);
    
    // extra metric; how many nodes do we expand above the apex?
    std::function<uint32_t(warthog::search_node*)> fn_get_apex = 
    [&chd, &fexp] (warthog::search_node* n) -> uint32_t
    {
        while(true)
        {
            warthog::search_node* p = fexp.generate(n->get_parent());
            if(!p || chd.level_->at(p->get_id()) < chd.level_->at(n->get_id()))
            { break; }
            n = p;
        }
        return chd.level_->at(n->get_id());
    };

    run_experiments(&alg, alg_name, parser, std::cout);
}

void
run_fch_bb_dfs(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
        std::string alg_name)
{
    std::string alg_params = cfg.get_param_value("alg");
    std::string chd_file = cfg.get_param_value("input");
    if(chd_file == "")
    {
        std::cerr << "err; require --input [chd file]\n";
        return;
    }

    // load up the graph 
    std::shared_ptr<warthog::ch::ch_data> 
        chd(warthog::ch::load_ch_data(chd_file.c_str(), warthog::ch::UP_DOWN));

    std::cerr << "preparing to search\n";
    // define the workload
    double cutoff = 1;
    if(alg_params != "")
    {
        int32_t pct_dijkstra = std::stoi(alg_params.c_str());
        if(!(pct_dijkstra >= 0 && pct_dijkstra <= 100))
        {
            std::cerr << "dijkstra percentage must be in range 0-100\n";
            return;
        }
        cutoff = pct_dijkstra > 0 ? (1 - ((double)pct_dijkstra)/100) : 0;
        alg_name += "-dijk-";
        alg_name += std::to_string(pct_dijkstra);
    }

    warthog::util::workload_manager workload(chd.get()->g_->get_num_nodes());
    for(size_t i = 0; i < chd.get()->g_->get_num_nodes(); i++)
    {
        if(chd.get()->level_->at(i) >= (uint32_t)(chd.get()->level_->size()*cutoff))
        { workload.set_flag((uint32_t)i, true); }
    }

    // load up the labelling
    std::string arclab_file =  chd_file + "." + alg_name + "." + "label";
    warthog::label::dfs_labelling* lab = warthog::label::dfs_labelling::load(
            arclab_file.c_str(), chd.get()->g_, chd.get()->level_);

    if(lab == 0)
    {
        lab = warthog::label::dfs_labelling::compute(
                chd.get()->g_, chd.get()->level_, &workload);
        std::cerr << "precompute finished. saving result to " 
            << arclab_file << "...";
        warthog::label::dfs_labelling::save(arclab_file.c_str(), *lab);
        std::cerr << "done.\n";
    }

    warthog::fch_dfs_expansion_policy fexp(chd.get()->g_, chd.get()->level_, lab);
    warthog::euclidean_heuristic h(chd.get()->g_);
    warthog::pqueue_min open;

    warthog::flexible_astar<
        warthog::euclidean_heuristic, 
        warthog::fch_dfs_expansion_policy,
        warthog::pqueue_min> 
            alg(&h, &fexp, &open);
    
    // extra metric; how many nodes do we expand above the apex?
    std::function<uint32_t(warthog::search_node*)> fn_get_apex = 
    [chd, &fexp] (warthog::search_node* n) -> uint32_t
    {
        while(true)
        {
            warthog::search_node* p = fexp.generate(n->get_parent());
            if(!p || chd.get()->level_->at(p->get_id()) < chd.get()->level_->at(n->get_id()))
            { break; }
            n = p;
        }
        return chd.get()->level_->at(n->get_id());
    };

    run_experiments(&alg, alg_name, parser, std::cout);
    delete lab;
}

void
run_fch_bb(warthog::util::cfg& cfg, warthog::dimacs_parser& parser, 
        std::string alg_name)
{
    std::string chd_file = cfg.get_param_value("input");
    std::string arclabels_file = cfg.get_param_value("input");
    if(chd_file == "" || arclabels_file == "")
    {
        std::cerr << "err; require --input [chd file] [arclabels file]\n";
        return;
    }

    // load up the graph 
    std::shared_ptr<warthog::ch::ch_data> 
        chd(warthog::ch::load_ch_data(chd_file.c_str(), warthog::ch::UP_DOWN));

    // load up the arc labels
    std::shared_ptr<warthog::label::bb_labelling> bbl
        (warthog::label::bb_labelling::load(arclabels_file.c_str(), chd.get()->g_));
    if(!bbl.get())
    {
        std::cerr << "err; could not load arcflags file\n";
        return;
    }
    warthog::bb_filter filter(bbl.get());

    warthog::euclidean_heuristic h(chd.get()->g_);
    warthog::fch_bb_expansion_policy fexp(chd.get()->g_, chd.get()->level_, &filter);
    warthog::pqueue_min open;

    warthog::flexible_astar< 
        warthog::euclidean_heuristic, 
        warthog::fch_bb_expansion_policy,
        warthog::pqueue_min>
            alg(&h, &fexp, &open);

    run_experiments(&alg, alg_name, parser, std::cout);

}

void
run_cpd_search(warthog::util::cfg& cfg, 
    warthog::dimacs_parser& parser, std::string alg_name)
{
    std::string xy_filename = cfg.get_param_value("input");
    if(xy_filename == "")
    {
        std::cerr << "parameter is missing: --input [xy-graph file]\n";
        return;
    }

    std::string cpd_filename = cfg.get_param_value("cpd");
    if(cpd_filename == "")
    {
        cpd_filename = xy_filename + ".cpd";
    }

    warthog::graph::xy_graph g;

    std::ifstream ifs(xy_filename);
    warthog::graph::read_xy(ifs, g);
    ifs.close();

    warthog::cpd::graph_oracle oracle(&g);
    ifs.open(cpd_filename);
    if(ifs.is_open())
    {
        ifs >> oracle;
    }
    else
    {
        oracle.precompute();
        std::ofstream ofs(cpd_filename);
        ofs << oracle;
        std::cerr << "writing " << cpd_filename << std::endl;
    }

//    warthog::simple_graph_expansion_policy expander(&g);
//    warthog::euclidean_heuristic h(&g);
//    warthog::pqueue_min open;
//
//    warthog::flexible_astar<
//        warthog::euclidean_heuristic, 
//        warthog::simple_graph_expansion_policy, 
//        warthog::pqueue_min> 
//            alg(&h, &expander, &open);
//
//    run_experiments(&alg, alg_name, parser, std::cout);
}

void
run_cpd(warthog::util::cfg& cfg, 
    warthog::dimacs_parser& parser, std::string alg_name)
{
    std::string xy_filename = cfg.get_param_value("input");
    if(xy_filename == "")
    {
        std::cerr << "parameter is missing: --input [xy-graph file]\n";
        return;
    }


    warthog::graph::xy_graph g;
    std::ifstream ifs(xy_filename);
    warthog::graph::read_xy(ifs, g);
    ifs.close();


    warthog::graph::node* n = g.get_node(433);
    for(uint32_t i = 0; i < n->out_degree(); i++)
    {
        warthog::graph::edge* e = n->outgoing_begin() + i;
        e->print(std::cerr);
        std::cerr << std::endl;
    }

    warthog::cpd::graph_oracle oracle(&g);

    std::string cpd_filename = cfg.get_param_value("cpd");
    if(cpd_filename == "")
    {
        cpd_filename = xy_filename + ".cpd";
    }

    ifs.open(cpd_filename);
    if(ifs.is_open())
    {
        ifs >> oracle;
    }
    else
    {
        std::cerr << "precomputing... " <<std::endl;
        oracle.precompute();
        std::ofstream ofs(cpd_filename);
        ofs << oracle;
        std::cerr << "writing " << cpd_filename << std::endl;
    }

    warthog::graph_oracle_expansion_policy expander(&oracle);
    warthog::zero_heuristic h;
    warthog::pqueue_min open;

    std::function<void(warthog::problem_instance*, warthog::solution*)> 
            cpd_extract = 
    [&oracle] (warthog::problem_instance* pi, warthog::solution* sol) -> void
    {
        sol->sum_of_edge_costs_ = 
            oracle.get_distance(pi->start_id_, pi->target_id_);
    };

    run_experiments(cpd_extract, alg_name, parser, std::cout);

}

void
run_dfs(warthog::util::cfg& cfg, 
    warthog::dimacs_parser& parser, std::string alg_name)
{
    std::string xy_filename = cfg.get_param_value("input");
    if(xy_filename == "")
    {
        std::cerr << "parameter is missing: --input [xy-graph file]\n";
        return;
    }

    warthog::graph::xy_graph g;

    std::ifstream ifs(xy_filename);
    warthog::graph::read_xy(ifs, g);
    ifs.close();

    warthog::simple_graph_expansion_policy expander(&g);
    warthog::zero_heuristic h;
    warthog::pqueue_min open;

    warthog::depth_first_search<
        warthog::zero_heuristic, 
        warthog::simple_graph_expansion_policy, 
        warthog::pqueue_min> 
            alg(&h, &expander, &open);

    run_experiments(&alg, alg_name, parser, std::cout);
}

void
run_dimacs(warthog::util::cfg& cfg)
{
    std::string alg_name = cfg.get_param_value("alg");
    std::string par_nruns = cfg.get_param_value("nruns");
    std::string problemfile = cfg.get_param_value("problem");

    if((alg_name == ""))
    {
        std::cerr << "parameter is missing: --alg\n";
        return;
    }

    if((problemfile == ""))
    {
        std::cerr << "parameter is missing: --problem\n";
        return;
    }

    if(par_nruns != "")
    {
       char* end;
       nruns = strtol(par_nruns.c_str(), &end, 10);
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
        run_dijkstra(cfg, parser, alg_name);
    }
    else if(alg_name == "astar")
    {
        run_astar(cfg, parser, alg_name);
    }
    else if(alg_name == "bi-dijkstra")
    {
        run_bi_dijkstra(cfg, parser, alg_name);
    }
    else if(alg_name == "bi-astar")
    {
        run_bi_astar(cfg, parser, alg_name);
    }
    else if(alg_name == "bch")
    {
        run_bch(cfg, parser, alg_name);
    }
    else if(alg_name == "bchb")
    {
        run_bch_backwards_only(cfg, parser, alg_name);
    }
    else if(alg_name == "bch-astar")
    {
        run_bch_astar(cfg, parser, alg_name);
    }
    else if(alg_name == "bch-bb")
    {
        run_bch_bb(cfg, parser, alg_name);
    }
    else if(alg_name == "fch")
    {
        run_fch(cfg, parser, alg_name);
    }
    else if(alg_name == "fch-bb")
    {
        run_fch_bb(cfg, parser, alg_name);
    }
    else if(alg_name == "fch-bb-dfs")
    {
        run_fch_bb_dfs(cfg, parser, alg_name);
    }
    else if(alg_name == "cpd-search")
    {
        run_cpd_search(cfg, parser, alg_name);
    }
    else if(alg_name == "cpd")
    {
        run_cpd(cfg, parser, alg_name);
    }
    else if(alg_name == "dfs")
    {
        run_dfs(cfg, parser, alg_name);
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



