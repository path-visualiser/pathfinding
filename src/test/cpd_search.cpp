#define CATCH_CONFIG_RUNNER

#include "catch.hpp"
#include "graph_expansion_policy.h"
#include "xy_graph.h"
#include "cpd_search.h"
#include "dimacs_parser.h"

using namespace std;
using namespace warthog;

int main(int argv, char* args[]) {
	Catch::Session session;
	int res = session.run(argv, args);
    return res;
}

SCENARIO("Test CPD A* on a square matrix", "[cpd][square][astar]")
{
    string map_name = "square01.map";
    graph::xy_graph g;
    dimacs_parser d(map_name.c_str());

    graph::dimacs_to_xy_graph(d, g, false, false, false);
    warthog::simple_graph_expansion_policy expander(&g);

    GIVEN("No CPD heuristic") {
        warthog::euclidean_heuristic h(&g);
        warthog::pqueue_min open;

        warthog::cpd_search<
            warthog::euclidean_heuristic,
            warthog::simple_graph_expansion_policy>
                astar(&h, &expander, &open);

        THEN("We can still search") {
            warthog::problem_instance pi(0, 19, true);
            warthog::solution sol;

            astar.get_path(pi, sol);
        }
    }
}
