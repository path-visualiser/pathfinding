// gridhog.cpp
//
// Pulls together a variety of different algorithms 
// for pathfinding on grid graphs.
//
// @author: dharabor
// @created: 2016-11-23
//

#include "cbs.h"
#include "cbs_ll_heuristic.h"
#include "cfg.h"
#include "constants.h"
#include "flexible_astar.h"
#include "gridmap.h"
#include "gridmap_expansion_policy.h"
#include "cbs_ll_expansion_policy.h"
#include "jps_expansion_policy.h"
#include "jps2_expansion_policy.h"
#include "jpsplus_expansion_policy.h"
#include "jps2plus_expansion_policy.h"
#include "ll_expansion_policy.h"
#include "manhattan_heuristic.h"
#include "octile_heuristic.h"
#include "scenario_manager.h"
#include "timer.h"
#include "labelled_gridmap.h"
#include "vl_gridmap_expansion_policy.h"
#include "zero_heuristic.h"

#include "getopt.h"

#include <fstream>
#include <functional>
#include <iomanip>
#include <sstream>
#include <unordered_map>
#include <memory>

#include "time_constraints.h"

// check computed solutions are optimal
int checkopt = 0;
// print debugging info during search
int verbose = 0;
// display program help on startup
int print_help = 0;

void
help()
{
	std::cerr << "valid parameters:\n"
	<< "\t--alg []\n"
	<< "\t--scen [scenario filename]\n"
	<< "\t--gen [map filename] \n"
	<< "\t--checkopt (optional)\n"
	<< "\t--verbose (optional)\n"
    << "\nRecognised values for --alg:\n"
    << "\tcbs_ll, cbs_ll_w, dijkstra, astar, astar_wgm, 4c_astar, sssp\n"
    << "\tjps, jps2, jps+, jps2+, jps\n"
    << "\tcpg, jpg\n";
}

bool
check_optimality(warthog::solution& sol, warthog::experiment* exp)
{
	uint32_t precision = 1;
	double epsilon = (1 / (int)pow(10, precision)) / 2;

	double delta = fabs(sol.sum_of_edge_costs_ - exp->distance());
	if( fabs(delta - epsilon) > epsilon)
	{
		std::stringstream strpathlen;
		strpathlen << std::fixed << std::setprecision(exp->precision());
		strpathlen << sol.sum_of_edge_costs_;

		std::stringstream stroptlen;
		stroptlen << std::fixed << std::setprecision(exp->precision());
		stroptlen << exp->distance();

		std::cerr << std::setprecision(exp->precision());
		std::cerr << "optimality check failed!" << std::endl;
		std::cerr << std::endl;
		std::cerr << "optimal path length: "<<stroptlen.str()
			<<" computed length: ";
		std::cerr << strpathlen.str()<<std::endl;
		std::cerr << "precision: " << precision << " epsilon: "<<epsilon<<std::endl;
		std::cerr<< "delta: "<< delta << std::endl;
		exit(1);
	}
    return true;
}

void
run_experiments(warthog::search* algo, std::string alg_name,
        warthog::scenario_manager& scenmgr, bool verbose, bool checkopt,
        std::ostream& out)
{
	std::cout 
        << "id\talg\texpanded\tinserted\tupdated\ttouched"
        << "\tnanos\tpcost\tplen\tmap\n";
	for(unsigned int i=0; i < scenmgr.num_experiments(); i++)
	{
		warthog::experiment* exp = scenmgr.get_experiment(i);

		uint32_t startid = exp->starty() * exp->mapwidth() + exp->startx();
		uint32_t goalid = exp->goaly() * exp->mapwidth() + exp->goalx();
        warthog::problem_instance pi(startid, goalid, verbose);
        warthog::solution sol;

        algo->get_path(pi, sol);

		out
            << i<<"\t" 
            << alg_name << "\t" 
            << sol.nodes_expanded_ << "\t" 
            << sol.nodes_inserted_ << "\t"
            << sol.nodes_updated_ << "\t"
            << sol.nodes_touched_ << "\t"
            << sol.time_elapsed_nano_ << "\t"
            << sol.sum_of_edge_costs_ << "\t" 
            << (sol.path_.size()-1) << "\t" 
            << scenmgr.last_file_loaded() 
            << std::endl;

        if(checkopt) { check_optimality(sol, exp); }
	}
}


void
run_jpsplus(warthog::scenario_manager& scenmgr, std::string alg_name)
{
    warthog::gridmap map(scenmgr.get_experiment(0)->map().c_str());
	warthog::jpsplus_expansion_policy expander(&map);
	warthog::octile_heuristic heuristic(map.width(), map.height());
    warthog::pqueue_min open;

	warthog::flexible_astar<
		warthog::octile_heuristic,
	   	warthog::jpsplus_expansion_policy,
        warthog::pqueue_min> 
            astar(&heuristic, &expander, &open);

    run_experiments(&astar, alg_name, scenmgr, 
            verbose, checkopt, std::cout);

	std::cerr << "done. total memory: "<< astar.mem() + scenmgr.mem() << "\n";
}

void
run_jps2plus(warthog::scenario_manager& scenmgr, std::string alg_name)
{
    warthog::gridmap map(scenmgr.get_experiment(0)->map().c_str());
	warthog::jps2plus_expansion_policy expander(&map);
	warthog::octile_heuristic heuristic(map.width(), map.height());
    warthog::pqueue_min open;

	warthog::flexible_astar<
		warthog::octile_heuristic,
	   	warthog::jps2plus_expansion_policy,
        warthog::pqueue_min> astar(&heuristic, &expander, &open);

    std::function<void(warthog::search_node*)> relax_fn = 
            [&](warthog::search_node* n)
            {
                expander.update_parent_direction(n);
            };
    astar.apply_on_relax(relax_fn);

	std::cerr << "done. total memory: "<< astar.mem() + scenmgr.mem() << "\n";
}

void
run_jps2(warthog::scenario_manager& scenmgr, std::string alg_name)
{
    warthog::gridmap map(scenmgr.get_experiment(0)->map().c_str());
	warthog::jps2_expansion_policy expander(&map);
	warthog::octile_heuristic heuristic(map.width(), map.height());
    warthog::pqueue_min open;

	warthog::flexible_astar<
		warthog::octile_heuristic,
	   	warthog::jps2_expansion_policy,
        warthog::pqueue_min> 
            astar(&heuristic, &expander, &open);

    std::function<void(warthog::search_node*)> relax_fn = 
            [&](warthog::search_node* n)
            {
                expander.update_parent_direction(n);
            };
    astar.apply_on_relax(relax_fn);

    run_experiments(&astar, alg_name, scenmgr, 
            verbose, checkopt, std::cout);
	std::cerr << "done. total memory: "<< astar.mem() + scenmgr.mem() << "\n";
}

void
run_jps(warthog::scenario_manager& scenmgr, std::string alg_name)
{
    warthog::gridmap map(scenmgr.get_experiment(0)->map().c_str());
	warthog::jps_expansion_policy expander(&map);
	warthog::octile_heuristic heuristic(map.width(), map.height());
    warthog::pqueue_min open;

	warthog::flexible_astar<
		warthog::octile_heuristic,
	   	warthog::jps_expansion_policy,
        warthog::pqueue_min> 
            astar(&heuristic, &expander, &open);

    run_experiments(&astar, alg_name, scenmgr, 
            verbose, checkopt, std::cout);
	std::cerr << "done. total memory: "<< astar.mem() + scenmgr.mem() << "\n";
}

void
run_astar(warthog::scenario_manager& scenmgr, std::string alg_name)
{
    warthog::gridmap map(scenmgr.get_experiment(0)->map().c_str());
	warthog::gridmap_expansion_policy expander(&map);
	warthog::octile_heuristic heuristic(map.width(), map.height());
    warthog::pqueue_min open;

	warthog::flexible_astar<
		warthog::octile_heuristic,
	   	warthog::gridmap_expansion_policy, 
        warthog::pqueue_min> 
            astar(&heuristic, &expander, &open);

    run_experiments(&astar, alg_name, scenmgr, 
            verbose, checkopt, std::cout);
	std::cerr << "done. total memory: "<< astar.mem() + scenmgr.mem() << "\n";
}

void
run_4c_astar(warthog::scenario_manager& scenmgr, std::string alg_name)
{
    warthog::gridmap map(scenmgr.get_experiment(0)->map().c_str());
	warthog::gridmap_expansion_policy expander(&map, true);
	warthog::manhattan_heuristic heuristic(map.width(), map.height());
    warthog::pqueue_min open;

	warthog::flexible_astar<
		warthog::manhattan_heuristic,
	   	warthog::gridmap_expansion_policy, 
        warthog::pqueue_min> 
            astar(&heuristic, &expander, &open);

    run_experiments(&astar, alg_name, scenmgr, 
            verbose, checkopt, std::cout);
	std::cerr << "done. total memory: "<< astar.mem() + scenmgr.mem() << "\n";
}

void
run_cbs_ll(warthog::scenario_manager& scenmgr, std::string alg_name)
{
    warthog::gridmap gm(scenmgr.get_experiment(0)->map().c_str());
	warthog::cbs_ll_heuristic heuristic;
	warthog::cbs_ll_expansion_policy expander(&gm, &heuristic);

    warthog::reservation_table restab(gm.width()*gm.height());
    warthog::cbs::cmp_cbs_ll_lessthan lessthan(&restab);
    warthog::cbs::pqueue_cbs_ll open(&lessthan);

	warthog::flexible_astar<
		warthog::cbs_ll_heuristic,
	   	warthog::cbs_ll_expansion_policy,
        warthog::cbs::pqueue_cbs_ll>
            astar(&heuristic, &expander, &open);

    // precompute heuristic values for each target location
    std::cerr << "precomputing heuristic values...";
    std::vector<uint32_t> target_locations;
	for(unsigned int i=0; i < scenmgr.num_experiments(); i++)
    {
		warthog::experiment* exp = scenmgr.get_experiment(i);
		uint32_t goalid = exp->goaly() * exp->mapwidth() + exp->goalx();
        target_locations.push_back(goalid);
    }
    heuristic.compute_h_values(target_locations, &gm);
    std::cerr << "done\n";

    run_experiments(&astar, alg_name, scenmgr, 
            verbose, checkopt, std::cout);
	std::cerr << "done. total memory: "<< astar.mem() + scenmgr.mem() << "\n";
}

// cbs low-level with variable edge costs
// (each action still takes one timestep, regardless of cost)
void
run_cbs_ll_w(warthog::scenario_manager& scenmgr, std::string alg_name)
{
    warthog::gridmap gm(scenmgr.get_experiment(0)->map().c_str());
	warthog::cbs_ll_heuristic heuristic;
	warthog::ll_expansion_policy expander(&gm, &heuristic);

    warthog::reservation_table restab(gm.width()*gm.height());
    warthog::cbs::cmp_cbs_ll_lessthan lessthan(&restab);
    warthog::cbs::pqueue_cbs_ll open(&lessthan);

	warthog::flexible_astar<
		warthog::cbs_ll_heuristic,
	   	warthog::ll_expansion_policy,
        warthog::cbs::pqueue_cbs_ll>
            astar(&heuristic, &expander, &open);

    // precompute heuristic values for each target location
    std::cerr << "precomputing heuristic values...";
    std::vector<uint32_t> target_locations;
	for(unsigned int i=0; i < scenmgr.num_experiments(); i++)
    {
		warthog::experiment* exp = scenmgr.get_experiment(i);
		uint32_t goalid = exp->goaly() * exp->mapwidth() + exp->goalx();
        target_locations.push_back(goalid);
    }
    heuristic.compute_h_values(target_locations, &gm);
    std::cerr << "done\n";

    run_experiments(&astar, alg_name, scenmgr, 
            verbose, checkopt, std::cout);
	std::cerr << "done. total memory: "<< astar.mem() + scenmgr.mem() << "\n";

}

void
run_dijkstra(warthog::scenario_manager& scenmgr, std::string alg_name)
{
    warthog::gridmap map(scenmgr.get_experiment(0)->map().c_str());
	warthog::gridmap_expansion_policy expander(&map);
	warthog::zero_heuristic heuristic;
    warthog::pqueue_min open;

	warthog::flexible_astar<
		warthog::zero_heuristic,
	   	warthog::gridmap_expansion_policy,
        warthog::pqueue_min> 
            astar(&heuristic, &expander, &open);

    run_experiments(&astar, alg_name, scenmgr, 
            verbose, checkopt, std::cout);
	std::cerr << "done. total memory: "<< astar.mem() + scenmgr.mem() << "\n";
}

void
run_wgm_astar(warthog::scenario_manager& scenmgr, std::string alg_name)
{
    warthog::vl_gridmap map(scenmgr.get_experiment(0)->map().c_str());
	warthog::vl_gridmap_expansion_policy expander(&map);
	warthog::octile_heuristic heuristic(map.width(), map.height());
    warthog::pqueue_min open;
    
    // cheapest terrain (movingai benchmarks) has ascii value '.'; we scale
    // all heuristic values accordingly (otherwise the heuristic doesn't 
    // impact f-values much and search starts to behave like dijkstra)
    heuristic.set_hscale('.');

	warthog::flexible_astar<
		warthog::octile_heuristic,
	   	warthog::vl_gridmap_expansion_policy,
        warthog::pqueue_min> 
            astar(&heuristic, &expander, &open);

    run_experiments(&astar, alg_name, scenmgr, 
            verbose, checkopt, std::cout);
	std::cerr << "done. total memory: "<< astar.mem() + scenmgr.mem() << "\n";
}

void
run_wgm_sssp(warthog::scenario_manager& scenmgr, std::string alg_name)
{
    warthog::vl_gridmap map(scenmgr.get_experiment(0)->map().c_str());
	warthog::vl_gridmap_expansion_policy expander(&map);
	warthog::zero_heuristic heuristic;
    warthog::pqueue_min open;

	warthog::flexible_astar<
		warthog::zero_heuristic,
	   	warthog::vl_gridmap_expansion_policy,
        warthog::pqueue_min> 
            astar(&heuristic, &expander, &open);

    run_experiments(&astar, alg_name, scenmgr, 
            verbose, checkopt, std::cout);
	std::cerr << "done. total memory: "<< astar.mem() + scenmgr.mem() << "\n";
}

void
run_sssp(warthog::scenario_manager& scenmgr, std::string alg_name)
{
    warthog::gridmap map(scenmgr.get_experiment(0)->map().c_str());
	warthog::gridmap_expansion_policy expander(&map);
	warthog::zero_heuristic heuristic;
    warthog::pqueue_min open;

	warthog::flexible_astar<
		warthog::zero_heuristic,
	   	warthog::gridmap_expansion_policy, 
        warthog::pqueue_min> 
            astar(&heuristic, &expander, &open);

    run_experiments(&astar, alg_name, scenmgr, 
            verbose, checkopt, std::cout);
	std::cerr << "done. total memory: "<< astar.mem() + scenmgr.mem() << "\n";
}

int 
main(int argc, char** argv)
{
	// parse arguments
	warthog::util::param valid_args[] = 
	{
		{"scen",  required_argument, 0, 0},
		{"alg",  required_argument, 0, 1},
		{"gen", required_argument, 0, 3},
		{"help", no_argument, &print_help, 1},
		{"checkopt",  no_argument, &checkopt, 1},
		{"verbose",  no_argument, &verbose, 1},
		{"format",  required_argument, 0, 1},
	};

	warthog::util::cfg cfg;
	cfg.parse_args(argc, argv, "-f", valid_args);

    if(argc == 1 || print_help)
    {
		help();
        exit(0);
    }

    std::string sfile = cfg.get_param_value("scen");
    std::string alg = cfg.get_param_value("alg");
    std::string gen = cfg.get_param_value("gen");

	if(gen != "")
	{
		warthog::scenario_manager sm;
		warthog::gridmap gm(gen.c_str());
		sm.generate_experiments(&gm, 1000) ;
		sm.write_scenario(std::cout);
        exit(0);
	}

	if(alg == "" || sfile == "")
	{
        std::cerr << "Err. Must specify a scenario file and search algorithm. Try --help for options.\n";
		exit(0);
	}

	warthog::scenario_manager scenmgr;
	scenmgr.load_scenario(sfile.c_str());

    if(alg == "jps+")
    {
        run_jpsplus(scenmgr, alg);
    }

    else if(alg == "jps2")
    {
        run_jps2(scenmgr, alg);
    }

    else if(alg == "jps2+")
    {
        run_jps2plus(scenmgr, alg);
    }

    else if(alg == "jps")
    {
        run_jps(scenmgr, alg);
    }

    else if(alg == "dijkstra")
    {
        run_dijkstra(scenmgr, alg); 
    }

    else if(alg == "astar")
    {
        run_astar(scenmgr, alg); 
    }
    else if(alg == "4c_astar")
    {
        run_4c_astar(scenmgr, alg); 
    }

    else if(alg == "cbs_ll")
    {
        run_cbs_ll(scenmgr, alg); 
    }
    else if(alg == "cbs_ll_w")
    {
        run_cbs_ll_w(scenmgr, alg); 
    }

    else if(alg == "astar_wgm")
    {
        run_wgm_astar(scenmgr, alg); 
    }

    else if(alg == "sssp")
    {
        run_sssp(scenmgr, alg);
    }

    else if(alg == "sssp")
    {
        run_wgm_sssp(scenmgr, alg); 
    }
    else
    {
        std::cerr << "err; invalid search algorithm: " << alg << "\n";
    }
}


