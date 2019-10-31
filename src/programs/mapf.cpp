// programs/mapf.cpp
// 
// Pulls together a variefy of algorithms for 
// solving MAPF instances 
//
// @author: dharabor
// @created: 2019-10-29
//

#include "cbs.h"
#include "cbs_ll_heuristic.h"
#include "cfg.h"
#include "constants.h"
#include "flexible_astar.h"
#include "gridmap.h"
#include "gridmap_expansion_policy.h"
#include "cbs_ll_expansion_policy.h"
#include "manhattan_heuristic.h"
#include "scenario_manager.h"
#include "timer.h"
#include "sipp_expansion_policy.h"

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
	<< "\t--alg [algorithm name]\n"
	<< "\t--scen [scenario filename]\n"
	<< "\t--gen [map filename] \n"
	<< "\t--verbose (optional)\n"
    << "\nRecognised values for --alg:\n"
    << "\tcbs_ll, sipp\n";
}

// run SIPP as a prioritised planning algorithm; assume each instance is an
// agent and plan them one after the other (i.e. the agent order is the
// input order that the instances appear in the scenario file. 
// run until completion or failure.
void
run_sipp(warthog::scenario_manager& scenmgr, std::string alg_name)
{
    warthog::gridmap gm(scenmgr.get_experiment(0)->map().c_str());
	warthog::manhattan_heuristic heuristic(gm.header_width(), gm.header_height());
    warthog::sipp_gridmap sipp_map(&gm);
	warthog::sipp_expansion_policy expander(&sipp_map);
    warthog::pqueue_min open;

	warthog::flexible_astar<
		warthog::manhattan_heuristic,
	   	warthog::sipp_expansion_policy,
        warthog::pqueue_min> astar(&heuristic, &expander, &open);

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

        astar.get_path(pi, sol);
		std::cout
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

        // the path of the agent now becomes an obstacle for 
        // the next agent. We assume that agents reach their
        // target and then disappear after one unit of time
        for(uint32_t i = 0; i < sol.path_.size()-1; i++)
        {
            int32_t x, y;
            expander.get_xy(sol.path_.at(i), x, y);

            int32_t nx, ny;
            expander.get_xy(sol.path_.at(i+1), nx, ny);

            warthog::cbs::move direction; 
            assert(nx == x || ny == y);
            if(nx == x && ny == y) {  direction = warthog::cbs::move::WAIT; }
            if(nx == x && ny < y) {  direction = warthog::cbs::move::NORTH; }
            if(nx == x && ny > y) {  direction = warthog::cbs::move::SOUTH; }
            if(nx < x && ny == y) {  direction = warthog::cbs::move::WEST; }
            if(nx > x && ny == y) {  direction = warthog::cbs::move::EAST; }

            warthog::search_node* current = expander.generate(sol.path_.at(i));
            warthog::search_node* next = expander.generate(sol.path_.at(i+1));
            warthog::cost_t start_time = current->get_g();
            warthog::cost_t end_time = next->get_g();
            sipp_map.add_obstacle(
                (uint32_t)x, (uint32_t)y, start_time, end_time, direction);
            //std::cerr << "adding obstacle: (" << x << ", " << y << ") @ (" 
            //          << start_time << ", " << end_time << ") dir " 
            //          << direction << std::endl;
            //if(start_time < 0 || end_time < 0)
            //{
            //    std::cerr << "wtf!\n";
            //}

            if(i == sol.path_.size()-2)
            {
                sipp_map.add_obstacle( (uint32_t)nx, (uint32_t)ny, 
                    end_time, end_time+1, warthog::cbs::move::WAIT);
            //    std::cerr 
            //        << "target obstacle: (" << nx << ", " << ny 
            //        << ") @ (" << end_time << ", " << end_time+1 
            //        << ") dir WAIT" << std::endl;
                break;
            }
        }
	}
    size_t max_intervals = 0;
    size_t tot_intervals = 0;
    uint32_t map_sz = sipp_map.gm_->header_height() * sipp_map.gm_->header_width();
    for(uint32_t i = 0; i < map_sz; i++)
    {
        max_intervals = sipp_map.get_all_intervals(i).size() > max_intervals ? 
                        sipp_map.get_all_intervals(i).size() : max_intervals;
        tot_intervals += sipp_map.get_all_intervals(i).size();
    }

	std::cerr << "done. \n"
              << "total memory: "<< astar.mem() + scenmgr.mem() << "\n";
	std::cerr << "max sipp intervals per node: "<< max_intervals << "\n";
	std::cerr << "avg sipp intervals per node: "<< tot_intervals / (double)map_sz << "\n";
}

void
run_cbs_ll(warthog::scenario_manager& scenmgr, std::string alg_name)
{
    warthog::gridmap gm(scenmgr.get_experiment(0)->map().c_str());
	warthog::cbs_ll_heuristic heuristic(&gm);
	warthog::cbs_ll_expansion_policy expander(&gm, &heuristic);
    warthog::pqueue_min open;

	warthog::flexible_astar<
		warthog::cbs_ll_heuristic,
	   	warthog::cbs_ll_expansion_policy,
        warthog::pqueue_min>
            astar(&heuristic, &expander, &open);

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

        astar.get_path(pi, sol);
		std::cout
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

        // the path of the agent now becomes an obstacle for 
        // the next agent. We assume that agents reach their
        // target and then disappear after one unit of time
        for(uint32_t i = 0; i < sol.path_.size()-1; i++)
        {
            int32_t x, y;
            expander.get_xy(sol.path_.at(i), x, y);

            int32_t nx, ny;
            expander.get_xy(sol.path_.at(i+1), nx, ny);

            // block each cell occupied by the agent
            expander.get_constraint(sol.path_.at(i))->v_ = true;
            
            // block any other agent from swapping positions with the agent
            // (i.e. prevent edge collisions) 
            if(nx !=  x || ny != y)
            {
                // compute the opposite direction
                warthog::cbs::move direction; 
                if(nx == x && ny < y) {  direction = warthog::cbs::move::SOUTH; }
                if(nx == x && ny > y) {  direction = warthog::cbs::move::NORTH; }
                if(nx < x && ny == y) {  direction = warthog::cbs::move::WEST; }
                if(nx > x && ny == y) {  direction = warthog::cbs::move::EAST; }

                uint32_t xy_id = sol.path_.at(i+1) & UINT32_MAX;
                uint64_t timestep = i;
                warthog::sn_id_t block_id = (timestep << 32) | xy_id;
                expander.get_constraint(block_id)->e_ 
                    |= (uint8_t)(1 << direction);
            }
        }
	}
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
        std::cerr << "Err. Must specify a scenario file and search algorithm. "
                  << "Try --help for options.\n";
		exit(0);
	}

	warthog::scenario_manager scenmgr;
	scenmgr.load_scenario(sfile.c_str());

    if(alg == "cbs_ll")
    {
        run_cbs_ll(scenmgr, alg); 
    }
    else if(alg == "sipp")
    {
        run_sipp(scenmgr, alg);
    }
    else
    {
        std::cerr << "err; invalid search algorithm: " << alg << "\n";
    }
}
