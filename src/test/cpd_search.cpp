#define CATCH_CONFIG_RUNNER

#include "catch.hpp"
#include "graph_expansion_policy.h"
#include "xy_graph.h"
#include "cpd_search.h"
#include "dimacs_parser.h"
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

    warthog::sn_id_t start = 0;
    warthog::sn_id_t goal = 19;
    // For some reason cannot use `warthog::sn_id_t`
    int32_t x, y;

    g.get_xy(start, x, y);
    REQUIRE(x == 0);
    REQUIRE(y == 0);
    g.get_xy(goal, x, y);
    REQUIRE(x == 400000);  // (5 - 1) * GRID_TO_GRAPH_SCALE_FACTOR
    REQUIRE(y == 400000);

    GIVEN("No CPD heuristic")
    {
        warthog::cpd_heuristic h(&g);
        warthog::pqueue_min open;
        warthog::solution sol;
        warthog::cpd_search<
            warthog::cpd_heuristic,
            warthog::simple_graph_expansion_policy>
                astar(&h, &expander, &open);

        THEN("We can still search")
        {
            warthog::problem_instance pi(start, goal, true);
            // Cannot cut corners
            warthog::cost_t cost = warthog::ONE *
                    (warthog::DBL_ONE * 6 + warthog::DBL_ROOT_TWO);

            astar.get_path(pi, sol);

            REQUIRE(sol.sum_of_edge_costs_ == cost);
            // But we have to expand all nodes but the target
            REQUIRE(sol.nodes_expanded_ == g.get_num_nodes() - 1);
        }

        WHEN("We do not have time to search")
        {
            astar.set_max_time_cutoff(0);

            THEN("We do not find a solution")
            {
                warthog::problem_instance pi(start, goal, true);
                astar.get_path(pi, sol);

                REQUIRE(sol.path_.empty());
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
    }
}
