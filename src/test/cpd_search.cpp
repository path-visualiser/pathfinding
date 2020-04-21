#define CATCH_CONFIG_RUNNER

#include "catch.hpp"
#include "graph_expansion_policy.h"
#include "xy_graph.h"
#include "cpd_search.h"
#include "dimacs_parser.h"
#include "graph_oracle.h"
#include "cpd_heuristic.h"

using namespace std;

int main(int argv, char* args[]) {
    Catch::Session session;
    int res = session.run(argv, args);
    return res;
}

SCENARIO("Test CPD A* on a square matrix", "[cpd][square][astar]")
{
    string map_name = "square01.map";
    warthog::graph::xy_graph g;
    warthog::gridmap d(map_name.c_str());

    warthog::graph::gridmap_to_xy_graph(&d, &g, false);
    warthog::simple_graph_expansion_policy expander(&g);
    warthog::cpd::graph_oracle oracle(&g);
    std::string cpd_filename = "square01.cpd";
    std::ifstream ifs(cpd_filename);

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

    warthog::sn_id_t start = 0;
    warthog::sn_id_t goal = 19;
    // For some reason cannot use `warthog::sn_id_t`
    int32_t x, y;
    // Cannot cut corners
    warthog::cost_t cost = warthog::ONE *
        (warthog::DBL_ONE * 6 + warthog::DBL_ROOT_TWO);

    g.get_xy(start, x, y);
    REQUIRE(x == 0);
    REQUIRE(y == 0);
    g.get_xy(goal, x, y);
    REQUIRE(x == 400000);  // (5 - 1) * GRID_TO_GRAPH_SCALE_FACTOR
    REQUIRE(y == 400000);

    GIVEN("A CPD heuristic")
    {
        warthog::cpd_heuristic h(&oracle);
        warthog::pqueue_min open;
        warthog::solution sol;
        warthog::cpd_search<
            warthog::cpd_heuristic,
            warthog::simple_graph_expansion_policy>
                astar(&h, &expander, &open);

        THEN("We should be searching quite fast")
        {
            warthog::problem_instance pi(start, goal, true);

            astar.get_path(pi, sol);

            REQUIRE(sol.sum_of_edge_costs_ == cost);
        }

        WHEN("We do not have time to search")
        {
            astar.set_max_time_cutoff(0);

            THEN("We should still find a solution because the start node " +
                 "has an incumbent.")
            {
                warthog::problem_instance pi(start, goal, true);
                astar.get_path(pi, sol);

                REQUIRE(sol.sum_of_edge_costs_ == cost);
            }
        }

        WHEN("We have a bit of time to search")
        {
            astar.set_max_ms_cutoff(1);

            THEN("We explore some of the tree")
            {
                warthog::problem_instance pi(start, goal, true);
                astar.get_path(pi, sol);

                REQUIRE(sol.nodes_expanded_ > 0);
            }
        }

        WHEN ("We limit the number of k moves")
        {
            astar.set_max_k_moves(2);

            THEN("We can still explore some")
            {
                warthog::problem_instance pi(start, goal, true);
                astar.get_path(pi, sol);

                // TODO have the actual number of nodes of some sort?
                REQUIRE(sol.nodes_expanded_ > 0);
            }
        }
    }

    GIVEN("A sub-optimal heuristic")
    {
        h.set_hscale(2.0);

        WHEN("The graph is perturbed within acceptable range")
        {
            // Modify optimal edge (0, 0) -> (0, 1), but the total cost of the
            // path is within the acceptance factor (100% in this case).
            perturb_edge(&g, 0, 1);

            THEN("The search is not modified")
            {
                warthog::problem_instance pi(start, goal, true);

                astar.get_path(pi, sol);

                REQUIRE(sol.path_ == optipath);
            }
        }
    }
}
