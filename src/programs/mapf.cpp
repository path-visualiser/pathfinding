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
#include "mapf/plan.h"
#include "scenario_manager.h"
#include "timer.h"
#include "sipp_expansion_policy.h"

#include "getopt.h"

#include <fstream>
#include <functional>
#include <iomanip>
#include <sstream>
#include <string>
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
	<< "\t--plan [plan filename (existing plan describing paths of higher priority agents)]\n"
	<< "\t--verbose (optional)\n"
    << "\nRecognised values for --alg:\n"
    << "\tcbs_ll, sipp\n";
}


// run SIPP as a prioritised planning algorithm; assume each instance is an
// agent and plan them one after the other (i.e. the agent order is the
// input order that the instances appear in the scenario file. 
// run until completion or failure.
void
run_sipp(warthog::scenario_manager& scenmgr, std::string alg_name, std::string plan_file)
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

    warthog::mapf::plan theplan;

    // this function blocks (== makes obstacles from) the planned paths of 
    // higher priority agents 
    auto add_higher_priority_plan = 
    [&sipp_map, &expander](warthog::solution& sol) -> void
    {
        for(uint32_t i = 0; i < sol.path_.size(); i++)
        {
            int32_t x, y;
            expander.get_xy(sol.path_.at(i).node_id_, x, y);

            int32_t nx, ny;
            expander.get_xy(sol.path_.at(i+1).node_id_, nx, ny);

            warthog::cbs::move direction; 

            assert(nx == x || ny == y);
            if(nx == x && ny == y) {  direction = warthog::cbs::move::WAIT; }
            if(nx == x && ny < y) {  direction = warthog::cbs::move::NORTH; }
            if(nx == x && ny > y) {  direction = warthog::cbs::move::SOUTH; }
            if(nx < x && ny == y) {  direction = warthog::cbs::move::WEST; }
            if(nx > x && ny == y) {  direction = warthog::cbs::move::EAST; }

            warthog::cost_t start_time = sol.path_.at(i).cost_;
            warthog::cost_t end_time = sol.path_.at(i+1).cost_;
            sipp_map.add_obstacle(
                (uint32_t)x, (uint32_t)y, start_time, end_time, direction);

            if(verbose)
            {
                std::cerr  << " add obstacle: (" << x << ", " << y << ") @ (" 
                           << start_time << ", " << end_time << ") dir " 
                           << direction << std::endl;
            }

            if(i == (sol.path_.size()-2))
            {
                sipp_map.add_obstacle( (uint32_t)nx, (uint32_t)ny, 
                    end_time, end_time+1, warthog::cbs::move::WAIT);
                if(verbose)
                {
                    std::cerr 
                        << "add obstacle: (" << nx << ", " << ny 
                        << ") @ (" << end_time << ", " << end_time+1 
                        << ") dir WAIT" << std::endl;
                }
                break;
            }
        }
    };

    // load plans of higher priority agents (if any) and block their
    // temporal locations to avoid collisions
    if(plan_file != "")
    {
        std::ifstream ifs(plan_file);
        warthog::mapf::plan hoplan; // higher priority plans
        ifs >> hoplan;
        for(uint32_t i = 0; i < hoplan.paths_.size(); i++)
        {
            add_higher_priority_plan(hoplan.paths_.at(i));
        }
        ifs.close();
    }

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

        // record the newly found plan and make this agent
        // higher priority than all subsequent agents 
        add_higher_priority_plan(sol);
        theplan.paths_.push_back(sol);
	}

    // some extra info about sipp's performance
    size_t max_intervals = 0;
    size_t tot_intervals = 0;
    uint32_t map_sz = sipp_map.gm_->header_height() * sipp_map.gm_->header_width();
    for(uint32_t i = 0; i < map_sz; i++)
    {
        max_intervals = sipp_map.get_all_intervals(i).size() > max_intervals ? 
                        sipp_map.get_all_intervals(i).size() : max_intervals;
        tot_intervals += sipp_map.get_all_intervals(i).size();
    }

    std::string tmp_planfile = scenmgr.last_file_loaded() + "." + alg_name + ".plan";
    std::cerr  << "writing plan to " << tmp_planfile << std::endl;
    std::ofstream ofs(tmp_planfile);
    ofs << theplan;
    ofs.close();
	std::cerr << "done. \n";
	std::cerr << "total memory: "<< astar.mem() + scenmgr.mem() << "\n";
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

    warthog::mapf::plan theplan;

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

        theplan.paths_.push_back(sol);

        // the path of the agent now becomes an obstacle for 
        // the next agent. We assume that agents reach their
        // target and then disappear after one unit of time
        for(uint32_t i = 0; i < sol.path_.size(); i++)
        {
            if(sol.path_.at(i).node_id_== sol.path_.back().node_id_) 
            { break; }

            int32_t x, y;
            expander.get_xy(sol.path_.at(i).node_id_, x, y);

            int32_t nx, ny;
            expander.get_xy(sol.path_.at(i+1).node_id_, nx, ny);

            // block each cell occupied by the agent
            expander.get_constraint(sol.path_.at(i).node_id_)->v_ = true;
            if(verbose)
            {
                std::cerr << "add obstacle (" << x << ", " << y << ") @ " 
                          << i << std::endl;
            }

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

                uint32_t xy_id = sol.path_.at(i+1).node_id_ & UINT32_MAX;
                uint64_t timestep = i;
                warthog::sn_id_t block_id = (timestep << 32) | xy_id;
                expander.get_constraint(block_id)->e_ 
                    |= (uint8_t)(1 << direction);
            }
        }
        // block target cell 
        if(sol.path_.size() > 0)
        {
            int32_t x, y;
            expander.get_xy(sol.path_.back().node_id_, x, y);
            expander.get_constraint(sol.path_.back().node_id_)->v_ = true;
            if(verbose)
            {
                std::cerr << "add obstacle (" << x << ", " << y << ") @ " << i << std::endl;
            }
        }
	}
    std::string tmp_planfile = scenmgr.last_file_loaded() + "." + alg_name + ".plan";
    std::cerr  << "writing plan to " << tmp_planfile << std::endl;
    std::ofstream ofs(tmp_planfile);
    ofs << theplan;
    ofs.close();
	std::cerr << "done. total memory: "<< astar.mem() + scenmgr.mem() << "\n";
}



int 
main(int argc, char** argv)
{
	// parse arguments
	warthog::util::param valid_args[] = 
	{
		{"scen",  required_argument, 0, 0},
		{"plan",  required_argument, 0, 0},
		{"alg",  required_argument, 0, 1},
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
    std::string planfile = cfg.get_param_value("plan");

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
        run_sipp(scenmgr, alg, planfile);
    }
    else
    {
        std::cerr << "err; invalid search algorithm: " << alg << "\n";
    }
}
