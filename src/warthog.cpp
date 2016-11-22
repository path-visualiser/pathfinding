// warthog.cpp
//
// @author: dharabor
// @created: August 2012
//

#include "arcflags_filter.h"
#include "bb_af_filter.h"
#include "bbox_filter.h"
#include "bidirectional_search.h"
#include "cfg.h"
#include "chaf_expansion_policy.h"
#include "chafbb_expansion_policy.h"
#include "chase_search.h"
#include "ch_expansion_policy.h"
#include "cpg_expansion_policy.h"
#include "constants.h"
#include "contraction.h"
#include "corner_graph.h"
#include "corner_graph_expansion_policy.h"
#include "corner_point_locator.h"
#include "dimacs_parser.h"
#include "dcl_filter.h"
#include "down_distance_filter.h"
#include "dummy_filter.h"
#include "euclidean_heuristic.h"
#include "fixed_graph_contraction.h"
#include "flexible_astar.h"
#include "fwd_ch_af_expansion_policy.h"
#include "fwd_ch_bbaf_expansion_policy.h"
#include "fwd_ch_expansion_policy.h"
#include "fwd_ch_dcl_expansion_policy.h"
#include "fwd_ch_dd_expansion_policy.h"
#include "fwd_ch_bb_expansion_policy.h"
#include "graph_expansion_policy.h"
#include "gridmap.h"
#include "gridmap_expansion_policy.h"
#include "jpg_expansion_policy.h"
#include "jps_expansion_policy.h"
#include "jps_expansion_policy_wgm.h"
#include "jps2_expansion_policy.h"
#include "jpsplus_expansion_policy.h"
#include "jps2plus_expansion_policy.h"
#include "lazy_graph_contraction.h"
#include "octile_heuristic.h"
#include "planar_graph.h"
#include "scenario_manager.h"
#include "timer.h"
#include "weighted_gridmap.h"
#include "wgridmap_expansion_policy.h"
#include "zero_heuristic.h"

#include "getopt.h"

#include <fstream>
#include <functional>
#include <iomanip>
#include <sstream>
#include <tr1/unordered_map>
#include <memory>

// check computed solutions are optimal
int checkopt = 0;
// print debugging info during search
int verbose = 0;
// display program help on startup
int print_help = 0;

// 0: map is a unweighted grid
// 1: map is a weighted grid
// 2: map is a directed graph in DIMACS format
int format;

void
help()
{
	std::cerr << "valid parameters:\n"
    << "--format [ grid | wgm ]; followed by\n"
	<< "\t--alg [astar | cpg | jpg | jps | jps2 | jps+ | jps2+ | jps | jps2_sssp | sssp ]\n"
	<< "\t--map [map filename] (only with dimacs format)\n"
	<< "\t--scen [scenario filename]\n"
	<< "\t--gen [map filename] (only with grid format)\n"
	<< "\t--checkopt (optional)\n"
    << "\n OR \n\n"
    << "--format [dimacs]; followed by\n"
    << "\t--input [gr file] [co file] [(optional) other input files]\n"
    << "\t--problem [ss or p2p problem file]\n"
    << "\t--alg [astar | dijkstra | bi-astar | bi-dijkstra | ch | chf | chf-dd ]\n"
	<< "\t--order [order-of-contraction file] (only with --alg ch)\n"
	<< "\t--verbose (optional)\n";
}

void
check_optimality(double len, warthog::experiment* exp)
{
	if(!checkopt)
	{
		return;
	}

	uint32_t precision = 1;
	double epsilon = (1 / (int)pow(10, precision)) / 2;

	int delta = abs(len - exp->distance());
	if( abs(delta - epsilon) > epsilon)
	{
		std::stringstream strpathlen;
		strpathlen << std::fixed << std::setprecision(exp->precision());
		strpathlen << len;

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
}

void
run_jpsplus(warthog::scenario_manager& scenmgr)
{
    warthog::gridmap map(scenmgr.get_experiment(0)->map().c_str());
	warthog::jpsplus_expansion_policy expander(&map);
	warthog::octile_heuristic heuristic(map.width(), map.height());

	warthog::flexible_astar<
		warthog::octile_heuristic,
	   	warthog::jpsplus_expansion_policy> astar(&heuristic, &expander);
	astar.set_verbose(verbose);

	std::cout << "id\talg\texpd\tgend\ttouched\ttime\tcost\tsfile\n";
	for(unsigned int i=0; i < scenmgr.num_experiments(); i++)
	{
		warthog::experiment* exp = scenmgr.get_experiment(i);

		int startid = exp->starty() * exp->mapwidth() + exp->startx();
		int goalid = exp->goaly() * exp->mapwidth() + exp->goalx();
		double len = astar.get_length(
				map.to_padded_id(startid),
			   	map.to_padded_id(goalid));
		if(len == warthog::INF)
		{
			len = 0;
		}

		std::cout << i<<"\t" << "jps+" << "\t" 
		<< astar.get_nodes_expanded() << "\t" 
		<< astar.get_nodes_generated() << "\t"
		<< astar.get_nodes_touched() << "\t"
		<< astar.get_search_time()  << "\t"
		<< len << "\t" 
		<< scenmgr.last_file_loaded() << std::endl;

		check_optimality(len, exp);
	}
	std::cerr << "done. total memory: "<< astar.mem() + scenmgr.mem() << "\n";
}

void
run_jpg(warthog::scenario_manager& scenmgr)
{

    std::shared_ptr<warthog::gridmap> map(
            new warthog::gridmap(
                scenmgr.get_experiment(0)->map().c_str()));

    std::shared_ptr<warthog::graph::corner_graph> cpg(
            new warthog::graph::corner_graph(map));

    warthog::jps::corner_point_locator cpl(map.get());
	warthog::jps::jpg_expansion_policy expander(cpg.get());

	warthog::octile_heuristic heuristic(map->width(), map->height());
	warthog::flexible_astar<
		warthog::octile_heuristic,
	   	warthog::jps::jpg_expansion_policy> 
            astar(&heuristic, &expander);
	astar.set_verbose(verbose);

	std::cout << "id\talg\texpd\tgend\ttouched\ttime\tcost\tsfile\n";
    warthog::timer mytimer;
	for(unsigned int i=0; i < scenmgr.num_experiments(); i++)
	{
		warthog::experiment* exp = scenmgr.get_experiment(i);
		uint32_t startid = map->to_padded_id( 
                exp->starty() * exp->mapwidth() + exp->startx());
		uint32_t goalid = map->to_padded_id(
                exp->goaly() * exp->mapwidth() + exp->goalx());

        mytimer.start();
        cpg->insert(startid, goalid);
        startid = cpg->get_inserted_start_id();
        goalid = cpg->get_inserted_target_id();
		double len = astar.get_length(startid, goalid);
        mytimer.stop();

		if(len == warthog::INF)
		{
			len = 0;
		}
        len /= warthog::SCALE_FACTOR_DBL_TO_INT;

		std::cout << i<<"\t" << "jpg" << "\t" 
		<< astar.get_nodes_expanded() << "\t" 
		<< astar.get_nodes_generated() << "\t"
		<< astar.get_nodes_touched() << "\t"
		<< mytimer.elapsed_time_micro() << "\t"
		<< len << "\t" 
		<< scenmgr.last_file_loaded() << std::endl;

		check_optimality(len, exp);
	}
	std::cerr << "done. total memory: "<< astar.mem() + scenmgr.mem() << "\n";
}

void
run_cpg(warthog::scenario_manager& scenmgr)
{
    std::shared_ptr<warthog::gridmap> map(
            new warthog::gridmap(scenmgr.get_experiment(0)->map().c_str()));
    std::shared_ptr<warthog::graph::corner_graph> cpg(
            new warthog::graph::corner_graph(map));
	warthog::jps::corner_graph_expansion_policy expander(cpg.get());

	warthog::octile_heuristic heuristic(map->width(), map->height());
	warthog::flexible_astar<
		warthog::octile_heuristic,
	   	warthog::jps::corner_graph_expansion_policy> 
            astar(&heuristic, &expander);
	astar.set_verbose(verbose);

	std::cout << "id\talg\texpd\tgend\ttouched\ttime\tcost\tsfile\n";
    warthog::timer mytimer;
	for(unsigned int i=0; i < scenmgr.num_experiments(); i++)
	{
		warthog::experiment* exp = scenmgr.get_experiment(i);
		uint32_t startid = map->to_padded_id( 
                exp->starty() * exp->mapwidth() + exp->startx());
		uint32_t goalid = map->to_padded_id(
                exp->goaly() * exp->mapwidth() + exp->goalx());

        mytimer.start();
        cpg->insert(startid, goalid);
        startid = cpg->get_inserted_start_id();
        goalid = cpg->get_inserted_target_id();
		double len = astar.get_length(startid, goalid);
        mytimer.stop();

		if(len == warthog::INF)
		{
			len = 0;
		}
        len /= warthog::SCALE_FACTOR_DBL_TO_INT;

		std::cout << i<<"\t" << "cpg" << "\t" 
		<< astar.get_nodes_expanded() << "\t" 
		<< astar.get_nodes_generated() << "\t"
		<< astar.get_nodes_touched() << "\t"
		<< mytimer.elapsed_time_micro() << "\t"
		<< len << "\t" 
		<< scenmgr.last_file_loaded() << std::endl;

		check_optimality(len, exp);
	}
	std::cerr << "done. total memory: "<< astar.mem() + scenmgr.mem() << "\n";
}

void
run_ch_cpg(warthog::scenario_manager& scenmgr, warthog::util::cfg& cfg)
{
    std::string grfile = cfg.get_param_value("gr");
    std::string cofile = cfg.get_param_value("co");
    std::string orderfile = cfg.get_param_value("order");
    if(!grfile.compare("") || !cofile.compare("") || !orderfile.compare(""))
    {
        std::cerr << "insufficient params; algo ch_cpg requires "
                  << "--gr [file] --co [file] --order [file]\n";
    }
    
    // load up the node order
    std::vector<uint32_t> order;
    warthog::ch::load_node_order(orderfile.c_str(), order);
    warthog::ch::value_index_swap_dimacs(order);

    // we insert two extra elements in the event that we
    // need to insert the start or target. both have the lowest
    // possible rank in the hierarchy (0 and 1)
    // NB: along the way we need to increase all ranks by 2 
    for(uint32_t i = 0; i < order.size(); i++) order.at(i)+=2;
    order.push_back(0);
    order.push_back(1);
    
    // load up the grid
    std::shared_ptr<warthog::gridmap> map(
            new warthog::gridmap(
                scenmgr.get_experiment(0)->map().c_str()));
            
    // load up the graph 
    std::shared_ptr<warthog::graph::planar_graph> pg(
            new warthog::graph::planar_graph());
    pg->load_dimacs(grfile.c_str(), cofile.c_str(), false, true);

    // wrapper
    warthog::graph::corner_graph cpg(map, pg);

    // expander
	warthog::zero_heuristic h;
	warthog::ch::cpg_expansion_policy fexp(&cpg, &order);
	warthog::ch::cpg_expansion_policy bexp(&cpg, &order, true);
    warthog::bidirectional_search<warthog::zero_heuristic> alg(&fexp, &bexp, &h);
	alg.set_verbose(verbose);

	std::cout << "id\talg\texpd\tgend\ttouched\ttime\tcost\tsfile\n";
    warthog::timer mytimer;
	for(unsigned int i=0; i < scenmgr.num_experiments(); i++)
	{
		warthog::experiment* exp = scenmgr.get_experiment(i);
		uint32_t startid = map->to_padded_id( 
                exp->starty() * exp->mapwidth() + exp->startx());
		uint32_t goalid = map->to_padded_id(
                exp->goaly() * exp->mapwidth() + exp->goalx());

        mytimer.start();
        cpg.insert(startid, goalid);
        startid = cpg.get_inserted_start_id();
        goalid = cpg.get_inserted_target_id();
		double len = alg.get_length(startid, goalid);
        mytimer.stop();

		if(len == warthog::INF)
		{
			len = 0;
		}
        len /= warthog::SCALE_FACTOR_DBL_TO_INT;

		std::cout << i<<"\t" << "ch_cpg" << "\t" 
		<< alg.get_nodes_expanded() << "\t" 
		<< alg.get_nodes_generated() << "\t"
		<< alg.get_nodes_touched() << "\t"
		<< mytimer.elapsed_time_micro() << "\t"
		<< len << "\t" 
		<< scenmgr.last_file_loaded() << std::endl;

		check_optimality(len, exp);
	}
	std::cerr << "done. total memory: "<< alg.mem() + scenmgr.mem() << "\n";
}

void
run_jps2plus(warthog::scenario_manager& scenmgr)
{
    warthog::gridmap map(scenmgr.get_experiment(0)->map().c_str());
	warthog::jps2plus_expansion_policy expander(&map);
	warthog::octile_heuristic heuristic(map.width(), map.height());

	warthog::flexible_astar<
		warthog::octile_heuristic,
	   	warthog::jps2plus_expansion_policy> astar(&heuristic, &expander);
	astar.set_verbose(verbose);

    std::function<void(warthog::search_node*)> relax_fn = 
            [&](warthog::search_node* n)
            {
                expander.update_parent_direction(n);
            };
    astar.apply_on_relax(relax_fn);

	std::cout << "id\talg\texpd\tgend\ttouched\ttime\tcost\tsfile\n";
	for(unsigned int i=0; i < scenmgr.num_experiments(); i++)
	{
		warthog::experiment* exp = scenmgr.get_experiment(i);

		int startid = exp->starty() * exp->mapwidth() + exp->startx();
		int goalid = exp->goaly() * exp->mapwidth() + exp->goalx();
		double len = astar.get_length(
				map.to_padded_id(startid),
			   	map.to_padded_id(goalid));
		if(len == warthog::INF)
		{
			len = 0;
		}

		std::cout << i<<"\t" << "jps2+" << "\t" 
		<< astar.get_nodes_expanded() << "\t" 
		<< astar.get_nodes_generated() << "\t"
		<< astar.get_nodes_touched() << "\t"
		<< astar.get_search_time()  << "\t"
		<< len << "\t" 
		<< scenmgr.last_file_loaded() << std::endl;

		check_optimality(len, exp);
	}
	std::cerr << "done. total memory: "<< astar.mem() + scenmgr.mem() << "\n";
}

void
run_jps2(warthog::scenario_manager& scenmgr)
{
    warthog::gridmap map(scenmgr.get_experiment(0)->map().c_str());
	warthog::jps2_expansion_policy expander(&map);
	warthog::octile_heuristic heuristic(map.width(), map.height());

	warthog::flexible_astar<
		warthog::octile_heuristic,
	   	warthog::jps2_expansion_policy> astar(&heuristic, &expander);
	astar.set_verbose(verbose);

    std::function<void(warthog::search_node*)> relax_fn = 
            [&](warthog::search_node* n)
            {
                expander.update_parent_direction(n);
            };
    astar.apply_on_relax(relax_fn);
            
	std::cout << "id\talg\texpd\tgend\ttouched\ttime\tcost\tsfile\n";
	for(unsigned int i=0; i < scenmgr.num_experiments(); i++)
	{
		warthog::experiment* exp = scenmgr.get_experiment(i);

		int startid = exp->starty() * exp->mapwidth() + exp->startx();
		int goalid = exp->goaly() * exp->mapwidth() + exp->goalx();
		double len = astar.get_length(
				map.to_padded_id(startid),
			   	map.to_padded_id(goalid));
		if(len == warthog::INF)
		{
			len = 0;
		}

		std::cout << i<<"\t" << "jps2" << "\t" 
		<< astar.get_nodes_expanded() << "\t" 
		<< astar.get_nodes_generated() << "\t"
		<< astar.get_nodes_touched() << "\t"
		<< astar.get_search_time()  << "\t"
		<< len << "\t" 
		<< scenmgr.last_file_loaded() << std::endl;

		check_optimality(len, exp);
	}
	std::cerr << "done. total memory: "<< astar.mem() + scenmgr.mem() << "\n";
}

void
run_jps(warthog::scenario_manager& scenmgr)
{
    warthog::gridmap map(scenmgr.get_experiment(0)->map().c_str());
	warthog::jps_expansion_policy expander(&map);
	warthog::octile_heuristic heuristic(map.width(), map.height());

	warthog::flexible_astar<
		warthog::octile_heuristic,
	   	warthog::jps_expansion_policy> astar(&heuristic, &expander);
	astar.set_verbose(verbose);

	std::cout << "id\talg\texpd\tgend\ttouched\ttime\tcost\tsfile\n";
	for(unsigned int i=0; i < scenmgr.num_experiments(); i++)
	{
		warthog::experiment* exp = scenmgr.get_experiment(i);

		int startid = exp->starty() * exp->mapwidth() + exp->startx();
		int goalid = exp->goaly() * exp->mapwidth() + exp->goalx();
		double len = astar.get_length(
				map.to_padded_id(startid),
			   	map.to_padded_id(goalid));
		if(len == warthog::INF)
		{
			len = 0;
		}

		std::cout << i<<"\t" << "jps" << "\t" 
		<< astar.get_nodes_expanded() << "\t" 
		<< astar.get_nodes_generated() << "\t"
		<< astar.get_nodes_touched() << "\t"
		<< astar.get_search_time()  << "\t"
		<< len << "\t" 
		<< scenmgr.last_file_loaded() << std::endl;

		check_optimality(len, exp);
	}
	std::cerr << "done. total memory: "<< astar.mem() + scenmgr.mem() << "\n";
}

void
run_jps2_sssp(warthog::scenario_manager& scenmgr)
{
    warthog::gridmap map(scenmgr.get_experiment(0)->map().c_str());
	warthog::jps2_expansion_policy expander(&map);
	//warthog::octile_heuristic heuristic(map.width(), map.height());
	warthog::zero_heuristic heuristic;

	warthog::flexible_astar<
		//warthog::octile_heuristic,
        warthog::zero_heuristic,
	   	warthog::jps2_expansion_policy> astar(&heuristic, &expander);
	astar.set_verbose(verbose);

	std::cout << "id\talg\texpd\tgend\ttouched\ttime\tcost\tsfile\n";
	for(unsigned int i=0; i < scenmgr.num_experiments(); i++)
	{
		warthog::experiment* exp = scenmgr.get_experiment(i);

		int startid = exp->starty() * exp->mapwidth() + exp->startx();
		int goalid = exp->mapwidth() * exp->mapheight();
		double len = astar.get_length(
				map.to_padded_id(startid),
			   	map.to_padded_id(goalid));
		if(len == warthog::INF)
		{
			len = 0;
		}

		std::cout << i<<"\t" << "jps2_sssp" << "\t" 
		<< astar.get_nodes_expanded() << "\t" 
		<< astar.get_nodes_generated() << "\t"
		<< astar.get_nodes_touched() << "\t"
		<< astar.get_search_time()  << "\t"
		<< len << "\t" 
		<< scenmgr.last_file_loaded() << std::endl;

		check_optimality(len, exp);
	}
	std::cerr << "done. total memory: "<< astar.mem() + scenmgr.mem() << "\n";
}

void
run_astar(warthog::scenario_manager& scenmgr)
{
    warthog::gridmap map(scenmgr.get_experiment(0)->map().c_str());
	warthog::gridmap_expansion_policy expander(&map);
	warthog::octile_heuristic heuristic(map.width(), map.height());

	warthog::flexible_astar<
		warthog::octile_heuristic,
	   	warthog::gridmap_expansion_policy> astar(&heuristic, &expander);
	astar.set_verbose(verbose);


	std::cout << "id\talg\texpd\tgend\ttouched\ttime\tcost\tsfile\n";
	for(unsigned int i=0; i < scenmgr.num_experiments(); i++)
	{
		warthog::experiment* exp = scenmgr.get_experiment(i);

		int startid = exp->starty() * exp->mapwidth() + exp->startx();
		int goalid = exp->goaly() * exp->mapwidth() + exp->goalx();
		double len = astar.get_length(
				map.to_padded_id(startid), 
				map.to_padded_id(goalid));
		if(len == warthog::INF)
		{
			len = 0;
		}

		std::cout << i<<"\t" << "astar" << "\t" 
		<< astar.get_nodes_expanded() << "\t" 
		<< astar.get_nodes_generated() << "\t"
		<< astar.get_nodes_touched() << "\t"
		<< astar.get_search_time()  << "\t"
		<< len << "\t" 
		<< scenmgr.last_file_loaded() << std::endl;

		check_optimality(len, exp);
	}
	std::cerr << "done. total memory: "<< astar.mem() + scenmgr.mem() << "\n";
}

void
run_wgm_astar(warthog::scenario_manager& scenmgr)
{
    warthog::weighted_gridmap map(scenmgr.get_experiment(0)->map().c_str());
	warthog::wgridmap_expansion_policy expander(&map);
	warthog::octile_heuristic heuristic(map.width(), map.height());
    
    // cheapest terrain (movingai benchmarks) has ascii value '.'; we scale
    // all heuristic values accordingly (otherwise the heuristic doesn't 
    // impact f-values much and search starts to behave like dijkstra)
    heuristic.set_hscale('.');

	warthog::flexible_astar<
		warthog::octile_heuristic,
	   	warthog::wgridmap_expansion_policy> astar(&heuristic, &expander);
	astar.set_verbose(verbose);

	std::cout << "id\talg\texpd\tgend\ttouched\ttime\tcost\tsfile\n";
	for(unsigned int i=0; i < scenmgr.num_experiments(); i++)
	{
		warthog::experiment* exp = scenmgr.get_experiment(i);

		int startid = exp->starty() * exp->mapwidth() + exp->startx();
		int goalid = exp->goaly() * exp->mapwidth() + exp->goalx();
		double len = astar.get_length(
				map.to_padded_id(startid), 
				map.to_padded_id(goalid));
		if(len == warthog::INF)
		{
			len = 0;
		}

		std::cout << i<<"\t" << "astar_wgm" << "\t" 
		<< astar.get_nodes_expanded() << "\t" 
		<< astar.get_nodes_generated() << "\t"
		<< astar.get_nodes_touched() << "\t"
		<< astar.get_search_time()  << "\t"
		<< len << "\t" 
		<< scenmgr.last_file_loaded() << std::endl;

		check_optimality(len, exp);
	}
	std::cerr << "done. total memory: "<< astar.mem() + scenmgr.mem() << "\n";
}

void
run_wgm_sssp(warthog::scenario_manager& scenmgr)
{
    warthog::weighted_gridmap map(scenmgr.get_experiment(0)->map().c_str());
	warthog::wgridmap_expansion_policy expander(&map);
	warthog::zero_heuristic heuristic;

	warthog::flexible_astar<
		warthog::zero_heuristic,
	   	warthog::wgridmap_expansion_policy> astar(&heuristic, &expander);
	astar.set_verbose(verbose);

	std::cout << "id\talg\texpd\tgend\ttouched\ttime\tsfile\n";
	for(unsigned int i=0; i < scenmgr.num_experiments(); i++)
	{
		warthog::experiment* exp = scenmgr.get_experiment(i);

		int startid = exp->starty() * exp->mapwidth() + exp->startx();
		astar.get_length(map.to_padded_id(startid), warthog::INF);

		std::cout << i<<"\t" << "sssp_wgm" << "\t" 
		<< astar.get_nodes_expanded() << "\t" 
		<< astar.get_nodes_generated() << "\t"
		<< astar.get_nodes_touched() << "\t"
		<< astar.get_search_time()  << "\t"
		<< scenmgr.last_file_loaded() << std::endl;
	}
	std::cerr << "done. total memory: "<< astar.mem() + scenmgr.mem() << "\n";
}

void
run_jps_wgm(warthog::scenario_manager& scenmgr)
{
    warthog::weighted_gridmap map(scenmgr.get_experiment(0)->map().c_str());
	warthog::jps_expansion_policy_wgm expander(&map);
	warthog::octile_heuristic heuristic(map.width(), map.height());
    // cheapest terrain (movingai benchmarks) has ascii value '.'; we scale
    // all heuristic values accordingly (otherwise the heuristic doesn't 
    // impact f-values much and search starts to behave like dijkstra)
    heuristic.set_hscale('.');  

	warthog::flexible_astar<
		warthog::octile_heuristic,
	   	warthog::jps_expansion_policy_wgm> astar(&heuristic, &expander);
	astar.set_verbose(verbose);

	std::cout << "id\talg\texpd\tgend\ttouched\ttime\tcost\tsfile\n";
	for(unsigned int i=0; i < scenmgr.num_experiments(); i++)
	{
		warthog::experiment* exp = scenmgr.get_experiment(i);

		int startid = exp->starty() * exp->mapwidth() + exp->startx();
		int goalid = exp->goaly() * exp->mapwidth() + exp->goalx();
		double len = astar.get_length(
				map.to_padded_id(startid),
			   	map.to_padded_id(goalid));
		if(len == warthog::INF)
		{
			len = 0;
		}

		std::cout << i<<"\t" << "jps_wgm" << "\t" 
		<< astar.get_nodes_expanded() << "\t" 
		<< astar.get_nodes_generated() << "\t"
		<< astar.get_nodes_touched() << "\t"
		<< astar.get_search_time()  << "\t"
		<< len << "\t" 
		<< scenmgr.last_file_loaded() << std::endl;

		check_optimality(len, exp);
	}
	std::cerr << "done. total memory: "<< astar.mem() + scenmgr.mem() << "\n";
}

uint32_t
get_apex(std::vector<uint32_t>& order, warthog::search_node* goal)
{
    std::function<uint32_t(warthog::search_node*)> fn_get_apex = 
        [&order] (warthog::search_node* n) -> uint32_t
        {
            if(n == 0) { return warthog::INF; } 

            while(true)
            {
                warthog::search_node* p = n->get_parent();
                if(!p || order.at(p->get_id()) < order.at(n->get_id()))
                { break; }
                n = p;
            }
            return order.at(n->get_id());
        };

    return fn_get_apex(goal);
}

// extra junk not used currently
//
//        std::function<uint32_t(uint32_t)> fn_redundant_expansions = 
//            [&fexp, &order, &alg] (uint32_t apex) -> uint32_t
//            {
//                std::set<uint32_t> tmp;
//                for(uint32_t i = 0; i < order.size(); i++)
//                {
//                    warthog::search_node* n = fexp.get_ptr(i, alg.get_searchid());
//                    if(n && order.at(n->get_id()) > apex)
//                    { 
//                        tmp.insert(n->get_id());
//                    }
//                }
//
//                for(uint32_t j = 0; j < order.size(); j++)
//                {
//                    warthog::search_node* n = fexp.get_ptr(j, alg.get_searchid());
//                    warthog::search_node* m = n;
//                    while(m)
//                    {
//                        if(tmp.find(m->get_id()) != tmp.end())
//                        {
//                            while(tmp.find(n->get_id()) == tmp.end())
//                            {
//                                tmp.insert(n->get_id());
//                                n = n->get_parent();
//                            }
//                            break;
//                        }
//                        m = m->get_parent();
//                    }
//                }
//                return tmp.size();
//            };

void
run_dimacs(warthog::util::cfg& cfg)
{

    std::string grfile = cfg.get_param_value("input");
    std::string cofile = cfg.get_param_value("input");
    std::string problemfile = cfg.get_param_value("problem");
    std::string alg_name = cfg.get_param_value("alg");


    if((problemfile == ""))
    {
        std::cerr << "parameter is missing: --problem\n";
        return;
    }
    if((grfile == ""))
    {
        std::cerr << "parameter is missing: --gr\n";
        return;
    }
    if((cofile == ""))
    {
        std::cerr << "parameter is missing: --co\n";
        return;
    }
    if((alg_name == ""))
    {
        std::cerr << "parameter is missing: --alg\n";
        return;
    }

    warthog::dimacs_parser parser;
    parser.load_instance(problemfile.c_str());

    if(alg_name == "dijkstra")
    {
        warthog::graph::planar_graph g;
        g.load_dimacs(grfile.c_str(), cofile.c_str());
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
            double len = alg.get_length(exp.source, (exp.p2p ? exp.target : warthog::INF));

            std::cout << i++ <<"\t" << alg_name << "\t" 
            << alg.get_nodes_expanded() << "\t" 
            << alg.get_nodes_generated() << "\t"
            << alg.get_nodes_touched() << "\t"
            << alg.get_search_time()  << "\t"
            << len << "\t" 
            << grfile << " " << problemfile << std::endl;
        }
    }
    else if(alg_name == "astar")
    {
        warthog::graph::planar_graph g;
        g.load_dimacs(grfile.c_str(), cofile.c_str());

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
            double len = alg.get_length(exp.source, (exp.p2p ? exp.target : warthog::INF));

            std::cout << i++ <<"\t" << alg_name << "\t" 
            << alg.get_nodes_expanded() << "\t" 
            << alg.get_nodes_generated() << "\t"
            << alg.get_nodes_touched() << "\t"
            << alg.get_search_time()  << "\t"
            << len << "\t" 
            << grfile << " " << problemfile << std::endl;
        }
    }
    else if(alg_name == "bi-dijkstra")
    {
        warthog::graph::planar_graph g;
        g.load_dimacs(grfile.c_str(), cofile.c_str());
        warthog::graph_expansion_policy fexp(&g);

        warthog::graph::planar_graph backward_g;
        backward_g.load_dimacs(grfile.c_str(), cofile.c_str(), true, true);
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
            double len = alg.get_length(exp.source, (exp.p2p ? exp.target : warthog::INF));

            std::cout << i++ <<"\t" << alg_name << "\t" 
            << alg.get_nodes_expanded() << "\t" 
            << alg.get_nodes_generated() << "\t"
            << alg.get_nodes_touched() << "\t"
            << alg.get_search_time()  << "\t"
            << len << "\t" 
            << grfile << " " << problemfile << std::endl;
        }
    }
    else if(alg_name == "bi-astar")
    {
        warthog::graph::planar_graph g;
        g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true);
        warthog::graph_expansion_policy fexp(&g);

        warthog::graph::planar_graph backward_g;
        backward_g.load_dimacs(grfile.c_str(), cofile.c_str(), true, true);
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
            double len = alg.get_length(exp.source, (exp.p2p ? exp.target : warthog::INF));

            std::cout << i++ <<"\t" << alg_name << "\t" 
            << alg.get_nodes_expanded() << "\t" 
            << alg.get_nodes_generated() << "\t"
            << alg.get_nodes_touched() << "\t"
            << alg.get_search_time()  << "\t"
            << len << "\t" 
            << grfile << " " << problemfile << std::endl;
        }
    }
    else if(alg_name == "ch")
    {
        std::string orderfile = cfg.get_param_value("order");

        // load up the graph 
        warthog::graph::planar_graph g;
        g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true);

        // load up the node order
        std::vector<uint32_t> order;
        warthog::ch::load_node_order(orderfile.c_str(), order);
        warthog::ch::value_index_swap_dimacs(order);

        std::cerr << "preparing to search\n";
        //warthog::zero_heuristic h;
        warthog::euclidean_heuristic h(&g);
        warthog::ch_expansion_policy fexp(&g, &order);
        warthog::ch_expansion_policy bexp (&g, &order, true);
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
            double len = alg.get_length(exp.source, (exp.p2p ? exp.target : warthog::INF));

            std::cout << i++ <<"\t" << alg_name << "\t" 
            << alg.get_nodes_expanded() << "\t" 
            << alg.get_nodes_generated() << "\t"
            << alg.get_nodes_touched() << "\t"
            << alg.get_search_time()  << "\t"
            << len << "\t" 
            << grfile << " " << problemfile << std::endl;
        }
    }
    else if(alg_name == "chase")
    {
        std::string orderfile = cfg.get_param_value("order");

        // load up the graph 
        warthog::graph::planar_graph g;
        g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true);

        // load up the node order
        std::vector<uint32_t> order;
        warthog::ch::load_node_order(orderfile.c_str(), order);
        warthog::ch::value_index_swap_dimacs(order);

        std::cerr << "preparing to search\n";
        //warthog::zero_heuristic h;
        warthog::euclidean_heuristic h(&g);
        warthog::ch_expansion_policy fexp(&g, &order);
        warthog::ch_expansion_policy bexp (&g, &order, true);
        //warthog::bidirectional_search<warthog::zero_heuristic> alg(&fexp, &bexp, &h);
        warthog::chase_search<warthog::euclidean_heuristic> alg(&fexp, &bexp, &h);
        alg.set_verbose(verbose);

        std::cerr << "running experiments\n";
        int i = 0;
        std::cout << "id\talg\texp\tgen\ttouch\tmicros\tplen\tmap\n";
        for(warthog::dimacs_parser::experiment_iterator it = parser.experiments_begin(); 
                it != parser.experiments_end(); it++)
        {
            warthog::dimacs_parser::experiment exp = (*it);
            double len = alg.get_length(exp.source, (exp.p2p ? exp.target : warthog::INF));

            std::cout << i++ <<"\t" << alg_name << "\t" 
            << alg.get_nodes_expanded() << "\t" 
            << alg.get_nodes_generated() << "\t"
            << alg.get_nodes_touched() << "\t"
            << alg.get_search_time()  << "\t"
            << len << "\t" 
            << grfile << " " << problemfile << std::endl;
        }
    }
    else if(alg_name == "chaf")
    {
        std::string orderfile = cfg.get_param_value("order");

        std::string arclabels_file = cfg.get_param_value("filter");
        std::cerr << "filter param " << arclabels_file << std::endl;

        std::string partition_file = cfg.get_param_value("filter");
        std::cerr << "filter param " << partition_file << std::endl;

        // load up the graph 
        warthog::graph::planar_graph g;
        g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true);

        // load up the node order
        std::vector<uint32_t> order;
        warthog::ch::load_node_order(orderfile.c_str(), order);
        warthog::ch::value_index_swap_dimacs(order);

        if(order.size() == 0)
        {
            std::cerr << "err; could not load node order file\n";
            return;
        }

        // load up the node partition info
        std::vector<uint32_t> part;
        warthog::helpers::load_integer_labels_dimacs(
                partition_file.c_str(), part);

        if(part.size() == 0)
        {
            std::cerr << "err; could not load node partition file\n";
            return;
        }
        
        // load up the arc-flags
        warthog::arcflags_filter filter(&g, &order, &part, arclabels_file.c_str());

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
            double len = alg.get_length(exp.source, (exp.p2p ? exp.target : warthog::INF));

            std::cout << i++ <<"\t" << alg_name << "\t" 
            << alg.get_nodes_expanded() << "\t" 
            << alg.get_nodes_generated() << "\t"
            << alg.get_nodes_touched() << "\t"
            << alg.get_search_time()  << "\t"
            << len << "\t" 
            << grfile << " " << problemfile << std::endl;
        }
    }
    else if(alg_name == "chaf-bb")
    {
        std::string orderfile = cfg.get_param_value("order");

        std::string arclabels_file = cfg.get_param_value("filter");
        std::cerr << "filter param " << arclabels_file << std::endl;

        std::string partition_file = cfg.get_param_value("filter");
        std::cerr << "filter param " << partition_file << std::endl;

        // load up the graph 
        warthog::graph::planar_graph g;
        g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true);

        // load up the node order
        std::vector<uint32_t> order;
        warthog::ch::load_node_order(orderfile.c_str(), order);
        warthog::ch::value_index_swap_dimacs(order);
        
        if(order.size() == 0)
        {
            std::cerr << "err; could not load node order file\n";
            return;
        }

        // load up the node partition info
        std::vector<uint32_t> part;
        warthog::helpers::load_integer_labels_dimacs(
                partition_file.c_str(), part);

        if(part.size() == 0)
        {
            std::cerr << "err; could not load node partition file\n";
            return;
        }
        
        // load up the arc labels
        warthog::bb_af_filter filter(&g, &order, &part, arclabels_file.c_str());

        std::cerr << "preparing to search\n";
        //warthog::zero_heuristic h;
        warthog::euclidean_heuristic h(&g);
        warthog::chafbb_expansion_policy fexp(&g, &order, &filter);
        warthog::chafbb_expansion_policy bexp (&g, &order, &filter, true);
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
            double len = alg.get_length(exp.source, (exp.p2p ? exp.target : warthog::INF));

            std::cout << i++ <<"\t" << alg_name << "\t" 
            << alg.get_nodes_expanded() << "\t" 
            << alg.get_nodes_generated() << "\t"
            << alg.get_nodes_touched() << "\t"
            << alg.get_search_time()  << "\t"
            << len << "\t" 
            << grfile << " " << problemfile << std::endl;
        }
    }
    else if(alg_name == "chf")
    {
        std::string orderfile = cfg.get_param_value("order");

        // load up the graph 
        warthog::graph::planar_graph g;
        g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true);

        // load up the node order
        std::vector<uint32_t> order;
        warthog::ch::load_node_order(orderfile.c_str(), order);
        warthog::ch::value_index_swap_dimacs(order);

        if(order.size() == 0)
        {
            std::cerr << "err; could not load node order file\n";
            return;
        }

        std::cerr << "preparing to search\n";
        warthog::fwd_ch_expansion_policy fexp(&g, &order); 
        warthog::euclidean_heuristic h(&g);


        warthog::flexible_astar< warthog::euclidean_heuristic, 
            warthog::fwd_ch_expansion_policy>
                alg(&h, &fexp);
        alg.set_verbose(verbose);

        std::cerr << "running experiments\n";
        int i = 0;
        std::cout << "id\talg\texp\tgen\ttouch\tmicros\tplen\tmap\n";
        for(warthog::dimacs_parser::experiment_iterator it = parser.experiments_begin(); 
                it != parser.experiments_end(); it++)
        {
            warthog::dimacs_parser::experiment exp = (*it);
            double len = alg.get_length(exp.source, (exp.p2p ? exp.target : warthog::INF));

            std::cout << i++ <<"\t" << alg_name << "\t" 
            << alg.get_nodes_expanded() << "\t" 
            << alg.get_nodes_generated() << "\t"
            << alg.get_nodes_touched() << "\t"
            << alg.get_search_time()  << "\t"
            << len << "\t" 
            << grfile << " " << problemfile << std::endl;
        }
    }
    else if(alg_name == "chf-dd")
    {
        std::string orderfile = cfg.get_param_value("order");
        std::string arclabels_file = cfg.get_param_value("filter");
        std::cerr << "filter param " << arclabels_file << std::endl;

        // load up the graph 
        warthog::graph::planar_graph g;
        g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true);

        // load up the node order
        std::vector<uint32_t> order;
        warthog::ch::load_node_order(orderfile.c_str(), order);
        warthog::ch::value_index_swap_dimacs(order);

        if(order.size() == 0)
        {
            std::cerr << "err; could not load node order file\n";
            return;
        }

        std::cerr << "preparing to search\n";
        warthog::down_distance_filter filter(arclabels_file.c_str(), &g, &order);
        warthog::euclidean_heuristic h(&g);
        warthog::fwd_ch_dd_expansion_policy fexp(&g, &order, &filter, &h);

        warthog::flexible_astar< warthog::euclidean_heuristic, 
           warthog::fwd_ch_dd_expansion_policy>
        alg(&h, &fexp);
        alg.set_verbose(verbose);

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
                    warthog::search_node* n = fexp.get_ptr(i, alg.get_searchid());
                    if(n && order.at(n->get_id()) > apex)
                    { 
                        tmp.insert(n->get_id());
                    }
                }

                for(uint32_t j = 0; j < order.size(); j++)
                {
                    warthog::search_node* n = fexp.get_ptr(j, alg.get_searchid());
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
            double len = alg.get_length(exp.source, (exp.p2p ? exp.target : warthog::INF));

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
            << grfile << " " << problemfile << std::endl;
        }
    }
    else if(alg_name == "chf-af")
    {
        std::string orderfile = cfg.get_param_value("order");

        // arclabels file
        std::string arclabels_file = cfg.get_param_value("filter");
        std::cerr << "filter param " << arclabels_file << std::endl;

        // partition file
        std::string partition_file = cfg.get_param_value("filter");
        std::cerr << "filter param " << partition_file << std::endl;

        std::cerr << "preparing to search\n";

        // load up the graph 
        warthog::graph::planar_graph g;
        g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true);

        // load up the node order
        std::vector<uint32_t> order;
        warthog::ch::load_node_order(orderfile.c_str(), order, true);

        if(order.size() == 0)
        {
            std::cerr << "err; could not load node order file\n";
            return;
        }

        // load up the node partition info
        std::vector<uint32_t> part;
        warthog::helpers::load_integer_labels_dimacs(
                partition_file.c_str(), part);

        if(part.size() == 0)
        {
            std::cerr << "err; could not load node partition file\n";
            return;
        }

        // load up the arc labels
        warthog::arcflags_filter filter(&g, &order, &part, arclabels_file.c_str());

        warthog::euclidean_heuristic h(&g);
        warthog::fwd_ch_af_expansion_policy fexp(&g, &order, &filter);

        warthog::flexible_astar< warthog::euclidean_heuristic, 
           warthog::fwd_ch_af_expansion_policy>
        alg(&h, &fexp);
        alg.set_verbose(verbose);

        std::cerr << "running experiments\n";
        int i = 0;
        std::cout << "id\talg\texp\tgen\ttouch\tmicros\tplen\tmap\n";
        for(warthog::dimacs_parser::experiment_iterator it = parser.experiments_begin(); 
                it != parser.experiments_end(); it++)
        {
            warthog::dimacs_parser::experiment exp = (*it);

            //fexp.set_apex(warthog::INF);
            double len = alg.get_length(exp.source, (exp.p2p ? exp.target : warthog::INF));
            std::cout << i <<"\t" << alg_name << "\t" 
            << alg.get_nodes_expanded() << "\t" 
            << alg.get_nodes_generated() << "\t"
            << alg.get_nodes_touched() << "\t"
            << alg.get_search_time()  << "\t"
            << len << "\t";
            std::cout << problemfile << std::endl;

            // run the experiment again but now suppose we have an 
            // oracle that tells us the apex of the path. this search
            // never expands down nodes before the apex and never expands
            // any up nodes after the apex
//            std::cout << "apex " << get_apex(order, fexp.get_ptr(exp.target, alg.get_searchid())) << std::endl;
//            fexp.set_apex(get_apex(order, fexp.get_ptr(exp.target, alg.get_searchid())));
//            len = alg.get_length(exp.source, (exp.p2p ? exp.target : warthog::INF));
//            std::cout << i <<"\t" << alg_name << "\t" 
//            << alg.get_nodes_expanded() << "\t" 
//            << alg.get_nodes_generated() << "\t"
//            << alg.get_nodes_touched() << "\t"
//            << alg.get_search_time()  << "\t"
//            << len << "\t";
//            std::cout << problemfile << std::endl;

            i++;
        }
    }
    else if(alg_name == "chf-bb")
    {
        std::string orderfile = cfg.get_param_value("order");
        std::string arclabels_file = cfg.get_param_value("filter");
        std::cerr << "filter param " << arclabels_file << std::endl;

        // load up the graph 
        warthog::graph::planar_graph g;
        g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true);

        // load up the node order
        std::vector<uint32_t> order;
        warthog::ch::load_node_order(orderfile.c_str(), order);
        warthog::ch::value_index_swap_dimacs(order);

        if(order.size() == 0)
        {
            std::cerr << "err; could not load node order file\n";
            return;
        }

        std::cerr << "preparing to search\n";
        warthog::bbox_filter filter(arclabels_file.c_str(), &g, &order);
        warthog::euclidean_heuristic h(&g);
        warthog::fwd_ch_bb_expansion_policy fexp(&g, &order, &filter);

        warthog::flexible_astar< warthog::euclidean_heuristic, 
           warthog::fwd_ch_bb_expansion_policy>
        alg(&h, &fexp);
        alg.set_verbose(verbose);

        std::cerr << "running experiments\n";
        int i = 0;
        std::cout << "id\talg\texp\tgen\ttouch\tmicros\tplen\tmap\n";
        for(warthog::dimacs_parser::experiment_iterator it = parser.experiments_begin(); 
                it != parser.experiments_end(); it++)
        {
            warthog::dimacs_parser::experiment exp = (*it);

            fexp.set_apex(warthog::INF);
            double len = alg.get_length(exp.source, (exp.p2p ? exp.target : warthog::INF));
            std::cout << i <<"\t" << alg_name << "\t" 
            << alg.get_nodes_expanded() << "\t" 
            << alg.get_nodes_generated() << "\t"
            << alg.get_nodes_touched() << "\t"
            << alg.get_search_time()  << "\t"
            << len << "\t";
            std::cout << problemfile << std::endl;

            // run the experiment again but now suppose we have an 
            // oracle that tells us the apex of the path. this search
            // never expands down nodes before the apex and never expands
            // any up nodes after the apex
//            fexp.set_apex(get_apex(order, fexp.get_ptr(exp.target, alg.get_searchid())));
//            len = alg.get_length(exp.source, (exp.p2p ? exp.target : warthog::INF));
//            std::cout << i <<"\t" << alg_name << "\t" 
//            << alg.get_nodes_expanded() << "\t" 
//            << alg.get_nodes_generated() << "\t"
//            << alg.get_nodes_touched() << "\t"
//            << alg.get_search_time()  << "\t"
//            << len << "\t";
//            std::cout << problemfile << std::endl;
//
            i++;
        }
    }
    else if(alg_name == "chf-bbaf")
    {
        std::string orderfile = cfg.get_param_value("order");

        // arclabels file
        std::string arclabels_file = cfg.get_param_value("filter");
        std::cerr << "filter param " << arclabels_file << std::endl;

        // partition file
        std::string partition_file = cfg.get_param_value("filter");
        std::cerr << "filter param " << partition_file << std::endl;

        std::cerr << "preparing to search\n";

        // load up the graph 
        warthog::graph::planar_graph g;
        g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true);

        // load up the node order
        std::vector<uint32_t> order;
        warthog::ch::load_node_order(orderfile.c_str(), order, true);

        if(order.size() == 0)
        {
            std::cerr << "err; could not load node order file\n";
            return;
        }

        // load up the node partition info
        std::vector<uint32_t> part;
        warthog::helpers::load_integer_labels_dimacs(
                partition_file.c_str(), part);

        if(part.size() == 0)
        {
            std::cerr << "err; could not load node partition file\n";
            return;
        }

        // load up the arc labels
        warthog::bb_af_filter filter(&g, &order, &part, arclabels_file.c_str());

        warthog::euclidean_heuristic h(&g);
        warthog::fwd_ch_bbaf_expansion_policy fexp(&g, &order, &filter);

        warthog::flexible_astar< warthog::euclidean_heuristic, 
           warthog::fwd_ch_bbaf_expansion_policy>
        alg(&h, &fexp);
        alg.set_verbose(verbose);

        std::cerr << "running experiments\n";
        int i = 0;
        std::cout << "id\talg\texp\tgen\ttouch\tmicros\tplen\tmap\n";
        for(warthog::dimacs_parser::experiment_iterator it = parser.experiments_begin(); 
                it != parser.experiments_end(); it++)
        {
            warthog::dimacs_parser::experiment exp = (*it);

            //fexp.set_apex(warthog::INF);
            double len = alg.get_length(exp.source, (exp.p2p ? exp.target : warthog::INF));
            std::cout << i <<"\t" << alg_name << "\t" 
            << alg.get_nodes_expanded() << "\t" 
            << alg.get_nodes_generated() << "\t"
            << alg.get_nodes_touched() << "\t"
            << alg.get_search_time()  << "\t"
            << len << "\t";
            std::cout << problemfile << std::endl;
            //std::cout << "apex " << get_apex(order, fexp.get_ptr(exp.target, alg.get_searchid())) << std::endl;

            // run the experiment again but now suppose we have an 
            // oracle that tells us the apex of the path. this search
            // never expands down nodes before the apex and never expands
            // any up nodes after the apex
//            fexp.set_apex(get_apex(order, fexp.get_ptr(exp.target, alg.get_searchid())));
//            len = alg.get_length(exp.source, (exp.p2p ? exp.target : warthog::INF));
//            std::cout << i <<"\t" << alg_name << "\t" 
//            << alg.get_nodes_expanded() << "\t" 
//            << alg.get_nodes_generated() << "\t"
//            << alg.get_nodes_touched() << "\t"
//            << alg.get_search_time()  << "\t"
//            << len << "\t";
//            std::cout << problemfile << std::endl;

            i++;
        }
    }
    else if(alg_name == "chf-dcl")
    {
        std::string orderfile = cfg.get_param_value("order");
        std::string arclabels_file = cfg.get_param_value("filter");
        std::cerr << "filter param " << arclabels_file << std::endl;

        // load up the graph 
        warthog::graph::planar_graph g;
        g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true);

        // load up the node order
        std::vector<uint32_t> order;
        warthog::ch::load_node_order(orderfile.c_str(), order);
        warthog::ch::value_index_swap_dimacs(order);

        if(order.size() == 0)
        {
            std::cerr << "err; could not load node order file\n";
            return;
        }

        std::cerr << "preparing to search\n";
        warthog::dcl_filter filter(arclabels_file.c_str(), &g, &order);
        warthog::euclidean_heuristic h(&g);
        warthog::fwd_ch_dcl_expansion_policy fexp(&g, &order, &filter);

        warthog::flexible_astar< warthog::euclidean_heuristic, 
           warthog::fwd_ch_dcl_expansion_policy>
        alg(&h, &fexp);
        alg.set_verbose(verbose);

        std::cerr << "running experiments\n";
        std::cout << "id\talg\texp\tgen\ttouch\tmicros\tplen\tmap\n";
        for(warthog::dimacs_parser::experiment_iterator it = parser.experiments_begin(); 
                it != parser.experiments_end(); it++)
        {
            warthog::dimacs_parser::experiment exp = (*it);
            fexp.set_apex(warthog::INF);

            int i = (it - parser.experiments_begin());
            double len = alg.get_length(exp.source, (exp.p2p ? exp.target : warthog::INF));
            std::cout << i <<"\t" << alg_name << "\t" 
            << alg.get_nodes_expanded() << "\t" 
            << alg.get_nodes_generated() << "\t"
            << alg.get_nodes_touched() << "\t"
            << alg.get_search_time()  << "\t"
            << len << "\t";
            std::cout << problemfile << std::endl;
            std::cout << "apex " << get_apex(order, fexp.get_ptr(exp.target, alg.get_searchid())) << std::endl;

            // run the experiment again but now suppose we have an 
            // oracle that tells us the apex of the path. this search
            // never expands down nodes before the apex and never expands
            // any up nodes after the apex
        //    fexp.set_apex(get_apex(order, fexp.get_ptr(exp.target, alg.get_searchid())));
        //    len = alg.get_length(exp.source, (exp.p2p ? exp.target : warthog::INF));
        //    std::cout << i <<"\t" << alg_name << "\t" 
        //    << alg.get_nodes_expanded() << "\t" 
        //    << alg.get_nodes_generated() << "\t"
        //    << alg.get_nodes_touched() << "\t"
        //    << alg.get_search_time()  << "\t"
        //    << len << "\t";
        //    std::cout << problemfile << std::endl;

        }
    }
    else
    {
        std::cerr << "invalid search algorithm\n";
    }
}

void
run_grid(warthog::util::cfg& cfg)
{
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

	// run experiments
	if(alg == "" || sfile == "")
	{
        std::cerr << "Err. Must specify a scenario file and search algorithm. Try --help for options.\n";
		exit(0);
	}

	warthog::scenario_manager scenmgr;
	scenmgr.load_scenario(sfile.c_str());

    int wgm = cfg.get_param_value("format") == "wgm";
    std::cerr << "wgm: " << (wgm ? "true" : "false") << std::endl;
    if(alg == "jps+")
    {
        run_jpsplus(scenmgr);
    }

    else if(alg == "jps2")
    {
        run_jps2(scenmgr);
    }

    else if(alg == "jps2+")
    {
        run_jps2plus(scenmgr);
    }

    else if(alg == "jps")
    {
        if(wgm)
        {
            run_jps_wgm(scenmgr);
        }
        else
        {
            run_jps(scenmgr);
        }
    }

    else if(alg == "astar")
    {
        if(wgm) 
        { 
            run_wgm_astar(scenmgr); 
        }
        else 
        { 
            run_astar(scenmgr); 
        }
    }

    else if(alg == "sssp")
    {
        if(wgm) 
        { 
            run_wgm_sssp(scenmgr); 
        }
        else 
        { 
            //run_astar(scenmgr); 
        }
    }
   
    else if(alg == "jps2_sssp")
    {
        run_jps2_sssp(scenmgr);
    }
    
    else if(alg == "jpg")
    {
        run_jpg(scenmgr);
    }

    else if(alg == "cpg")
    {
        run_cpg(scenmgr);
    }

    else if(alg == "ch_cpg")
    {
        run_ch_cpg(scenmgr, cfg);
    }
    else
    {
        std::cerr << "err; invalid search algorithm: " << alg << "\n";
    }
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
		{"gr",  required_argument, 0, 1},
		{"co",  required_argument, 0, 1},
		{"problem",  required_argument, 0, 1},
		{"format",  required_argument, 0, 1},
		{"order",  required_argument, 0, 1},
		{"part",  required_argument, 0, 1},
        {"filter", required_argument, 0, 1}
	};

	warthog::util::cfg cfg;
	cfg.parse_args(argc, argv, "-f", valid_args);

    if(argc == 1 || print_help)
    {
		help();
        exit(0);
    }


    std::string format = cfg.get_param_value("format");
    if(format == "dimacs")
    {
        run_dimacs(cfg);
    }
    if((format == "grid") || (format == "wgm"))
    {
        run_grid(cfg);
    }
}


