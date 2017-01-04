#include "af_filter.h"
#include "afh_filter.h"
#include "afhd_filter.h"
#include "bbaf_filter.h"
#include "bb_filter.h"
#include "cfg.h"
#include "ch_expansion_policy.h"
#include "dimacs_parser.h"
#include "dcl_filter.h"
#include "down_distance_filter.h"
#include "fixed_graph_contraction.h"
#include "graph.h"
#include "helpers.h"
#include "lazy_graph_contraction.h"
#include "planar_graph.h"

#include <iostream>
#include <string>

int verbose=false;
int down_dist_first_id;
int down_dist_last_id;
warthog::util::cfg cfg;

void
help()
{
    std::cerr << 
        "create arc labels for " <<
        "a given (currently, DIMACS-format only) input graph\n";
	std::cerr << "valid parameters:\n"
    //<< "\t--dimacs [gr file] [co file] (IN THIS ORDER!!)\n"
	//<< "\t--order [order-of-contraction file]\n"
    << "\t--type [ downdist | dcl | af | bb | bbaf | chaf | chbb | chbbaf | afh | afhd ]\n"
    << "\t--input [ algorithm-specific input files (omit to show options) ]\n" 
	<< "\t--verbose (optional)\n";
}

void 
compute_down_distance()
{
    std::string grfile = cfg.get_param_value("input");
    std::string cofile = cfg.get_param_value("input");
    std::string orderfile = cfg.get_param_value("input");
    if(grfile == "" || cofile == "" || orderfile == "")
    {
        std::cerr << "err; insufficient input parameters."
                  << " required, in order:\n"
                  << " --input [gr file] [co file] "
                  << " [contraction order file]\n";
        return;
    }

    // load up the graph
    warthog::graph::planar_graph g;
    if(!g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true))
    {
        std::cerr << "err; could not load gr or co input files (one or both)\n";
        return;
    }

    // load up the node ordering
    std::vector<uint32_t> order;
    if(!warthog::ch::load_node_order(orderfile.c_str(), order, true))
    {
        std::cerr << "err; could not load node order input file\n";
        return;
    }

    // compute labels
    grfile.append(".ddist.arclabel");
    uint32_t firstid = 0;
    uint32_t lastid = g.get_num_nodes()-1;
    if(cfg.get_num_values("type") == 2)
    {
        std::string first = cfg.get_param_value("type");
        std::string last = cfg.get_param_value("type");
        if(strtol(first.c_str(), 0, 10) != 0)
        {
            firstid = strtol(first.c_str(), 0, 10);
        }
        if(strtol(last.c_str(), 0, 10) != 0)
        {
            lastid = strtol(last.c_str(), 0, 10);
        }
        grfile.append(".");
        grfile.append(first);
        grfile.append(".");
        grfile.append(std::to_string(lastid));
    }
    warthog::down_distance_filter ddfilter(&g);
    ddfilter.compute(firstid, lastid, &order);
    
    // save the result
    std::cerr << "saving contracted graph to file " << grfile << std::endl;
    std::fstream out(grfile.c_str(), std::ios_base::out | std::ios_base::trunc);
    if(!out.good())
    {
        std::cerr << "\nerror exporting ch to file " << grfile << std::endl;
    }
    ddfilter.print(out);
    out.close();
    std::cerr << "all done!\n";
}

void 
compute_dcl_labels()
{
    std::string grfile = cfg.get_param_value("input");
    std::string cofile = cfg.get_param_value("input");
    std::string orderfile = cfg.get_param_value("input");

    if(grfile == "" || cofile == "" || orderfile == "")
    {
        std::cerr << "err; insufficient input parameters."
                  << " required, in order:\n"
                  << " --input [gr file] [co file] "
                  << " [contraction order file]\n";
        return;
    }

    // load up (or create) the contraction hierarchy
    warthog::graph::planar_graph g;
    if(!g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true))
    {
        std::cerr << "err; could not load gr or co input files (one or both)\n";
        return;
    }

    // load up the contraction order
    std::vector<uint32_t> order;
    if(!warthog::ch::load_node_order(orderfile.c_str(), order, true))
    {
        std::cerr << "err; could not load node order input file\n";
        return;
    }

    // compute down closure
    grfile.append(".ch-dcl.arclabel");
    warthog::dcl_filter filter(&g);
    filter.compute(&order);
    
    // save the result
    std::cerr << "saving contracted graph to file " << grfile << std::endl;
    std::fstream out(grfile.c_str(), std::ios_base::out | std::ios_base::trunc);
    if(!out.good())
    {
        std::cerr << "\nerror exporting ch to file " << grfile << std::endl;
    }
    filter.print(out);
    out.close();
    std::cerr << "all done!\n";
}

void
compute_chbb_labels()
{
    std::string grfile = cfg.get_param_value("input");
    std::string cofile = cfg.get_param_value("input");
    std::string orderfile = cfg.get_param_value("input");

    if( grfile == "" || cofile == "" || orderfile == "")
    {
        std::cerr << "err; insufficient input parameters."
                  << " required, in order:\n"
                  << " --input [gr file] [co file]"
                  << " [node ordering file]\n";
        return;
    }
    std::cerr << "computing labels" << std::endl;

    // load up the graph
    warthog::graph::planar_graph g;
    if(!g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true))
    {
        std::cerr << "err; could not load gr or co input files (one or both)\n";
        return;
    }

    // load up the node ordering
    std::vector<uint32_t> order;
    if(!warthog::ch::load_node_order(orderfile.c_str(), order, true))
    {
        std::cerr << "err; could not load node order input file\n";
        return;
    }
    
    // compute labels
    uint32_t firstid = 0;
    uint32_t lastid = g.get_num_nodes()-1;
    grfile.append(".ch-bb.arclabel");
    if(cfg.get_num_values("type") == 2)
    {
        std::string first = cfg.get_param_value("type");
        std::string last = cfg.get_param_value("type");
        if(strtol(first.c_str(), 0, 10) != 0)
        {
            firstid = strtol(first.c_str(), 0, 10);
        }
        if(strtol(last.c_str(), 0, 10) != 0)
        {
            lastid = strtol(last.c_str(), 0, 10);
        }
        grfile.append(".");
        grfile.append(first);
        grfile.append(".");
        grfile.append(std::to_string(lastid));
    }
    
    // compute down_distance
    warthog::bb_filter filter(&g);
    filter.compute_ch(firstid, lastid, &order);
    
    // save the result
    std::cerr << "saving contracted graph to file " << grfile << std::endl;
    std::fstream out(grfile.c_str(), std::ios_base::out | std::ios_base::trunc);
    if(!out.good())
    {
        std::cerr << "\nerror exporting ch to file " << grfile << std::endl;
    }
    filter.print(out);
    out.close();
    std::cerr << "all done!\n";
}

void
compute_bb_labels()
{
    std::string grfile = cfg.get_param_value("input");
    std::string cofile = cfg.get_param_value("input");
    if(grfile == "" || cofile == "")
    {
        std::cerr << "err; insufficient input parameters."
                  << " required, in order:\n"
                  << " --input [gr file] [co file]\n";
        return;
    }

    // load up (or create) the graph
    warthog::graph::planar_graph g;
    if(!g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true))
    {
        std::cerr << "err; could not load gr or co input files (one or both)\n";
        return;
    }

    // compute bounding box labels

    grfile.append(".bb.arclabel");
    uint32_t firstid = 0;
    uint32_t lastid = g.get_num_nodes();

    if(cfg.get_num_values("type") == 2)
    {
        std::string first = cfg.get_param_value("type");
        std::string last = cfg.get_param_value("type");
        if(strtol(first.c_str(), 0, 10) != 0)
        {
            firstid = strtol(first.c_str(), 0, 10);
        }
        if(strtol(last.c_str(), 0, 10) != 0)
        {
            lastid = strtol(last.c_str(), 0, 10);
        }
        grfile.append(".");
        grfile.append(first);
        grfile.append(".");
        grfile.append(std::to_string(lastid));
    }
    warthog::bb_filter filter(&g);
    filter.compute(firstid, lastid);
    
    // save the result
    std::cerr << "saving contracted graph to file " << grfile << std::endl;
    std::fstream out(grfile.c_str(), std::ios_base::out | std::ios_base::trunc);
    if(!out.good())
    {
        std::cerr << "\nerror exporting ch to file " << grfile << std::endl;
    }
    filter.print(out);
    out.close();
    std::cerr << "all done!\n";
}

void
compute_chaf_labels()
{
    std::string grfile = cfg.get_param_value("input");
    std::string cofile = cfg.get_param_value("input");
    std::string partfile = cfg.get_param_value("input");
    std::string orderfile = cfg.get_param_value("input");

    if( grfile == "" || cofile == "" || partfile == "" || orderfile == "")
    {
        std::cerr << "err; insufficient input parameters."
                  << " required, in order:\n"
                  << " --input [gr file] [co file]"
                  << " [graph partition file]"
                  << " [node ordering file]\n";
        return;
    }
    std::cerr << "computing labels" << std::endl;

    // load up the graph
    warthog::graph::planar_graph g;
    g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true);

    // load up the partition info
    std::vector<uint32_t> part;
    warthog::helpers::load_integer_labels_dimacs(partfile.c_str(), part);

    // load up the node ordering
    std::vector<uint32_t> order;
    if(!warthog::ch::load_node_order(orderfile.c_str(), order, true))
    {
        std::cerr << "err; could not load node order input file\n";
        return;
    }

    // compute labels
    uint32_t firstid = 0;
    uint32_t lastid = g.get_num_nodes()-1;
    grfile.append(".ch-af.arclabel");
    if(cfg.get_num_values("type") == 2)
    {
        std::string first = cfg.get_param_value("type");
        std::string last = cfg.get_param_value("type");
        if(strtol(first.c_str(), 0, 10) != 0)
        {
            firstid = strtol(first.c_str(), 0, 10);
        }
        if(strtol(last.c_str(), 0, 10) != 0)
        {
            lastid = strtol(last.c_str(), 0, 10);
        }
        grfile.append(".");
        grfile.append(first);
        grfile.append(".");
        grfile.append(std::to_string(lastid));
    }
    warthog::af_filter filter(&g, &part);
    filter.compute_ch(firstid, lastid, &order);
    
    // save the result
    std::cerr << "saving contracted graph to file " << grfile << std::endl;
    std::fstream out(grfile.c_str(), std::ios_base::out | std::ios_base::trunc);
    if(!out.good())
    {
        std::cerr << "\nerror exporting ch to file " << grfile << std::endl;
    }
    filter.print(out);
    out.close();
    std::cerr << "all done!\n";
}

void
compute_af_labels()
{
    std::string grfile = cfg.get_param_value("input");
    std::string cofile = cfg.get_param_value("input");
    std::string partfile = cfg.get_param_value("input");

    if( grfile == "" || cofile == "" || partfile == "")
    {
        std::cerr << "err; insufficient input parameters."
                  << " required, in order:\n"
                  << " --input [gr file] [co file]"
                  << " [graph partition file]\n";
        return;
    }
    std::cerr << "computing labels" << std::endl;

    // load up (or create) the graph
    warthog::graph::planar_graph g;
    if(!g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true))
    {
        std::cerr << "err; could not load gr or co input files (one or both)\n";
        return;
    }

    // load up the partition info
    std::vector<uint32_t> part;
    if(!warthog::helpers::load_integer_labels_dimacs(partfile.c_str(), part))
    {
        std::cerr << "err; could not load partition file\n"; 
        return;
    }


    // compute labels
    uint32_t firstid = 0;
    uint32_t lastid = g.get_num_nodes()-1;
    grfile.append(".af.arclabel");
    if(cfg.get_num_values("type") == 2)
    {
        std::string first = cfg.get_param_value("type");
        std::string last = cfg.get_param_value("type");
        if(strtol(first.c_str(), 0, 10) != 0)
        {
            firstid = strtol(first.c_str(), 0, 10);
        }
        if(strtol(last.c_str(), 0, 10) != 0)
        {
            lastid = strtol(last.c_str(), 0, 10);
        }
        grfile.append(".");
        grfile.append(first);
        grfile.append(".");
        grfile.append(std::to_string(lastid));
    }
    warthog::af_filter filter(&g, &part);
    filter.compute(firstid, lastid);
    
    // save the result
    std::cerr << "saving contracted graph to file " << grfile << std::endl;
    std::fstream out(grfile.c_str(), std::ios_base::out | std::ios_base::trunc);
    if(!out.good())
    {
        std::cerr << "\nerror exporting ch to file " << grfile << std::endl;
    }
    filter.print(out);
    out.close();
    std::cerr << "all done!\n";
}

void
compute_bbaf_labels()
{
    std::string grfile = cfg.get_param_value("input");
    std::string cofile = cfg.get_param_value("input");
    std::string partfile = cfg.get_param_value("input");

    if( grfile == "" || cofile == "" || partfile == "")
    {
        std::cerr << "err; insufficient input parameters."
                  << " required, in order:\n"
                  << " --input [gr file] [co file]"
                  << " [graph partition file]\n";
        return;
    }
    std::cerr << "computing labels" << std::endl;

    // load up the graph
    warthog::graph::planar_graph g;
    if(!g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true))
    {
        std::cerr << "err; could not load gr or co input files (one or both)\n";
        return;
    }

    // load up the partition info
    std::vector<uint32_t> part;
    if(!warthog::helpers::load_integer_labels_dimacs(partfile.c_str(), part))
    {
        std::cerr << "err; could not load partition file\n"; 
        return;
    }

    // compute labels
    uint32_t firstid = 0;
    uint32_t lastid = g.get_num_nodes()-1;
    grfile.append(".bbaf.arclabel");
    if(cfg.get_num_values("type") == 2)
    {
        std::string first = cfg.get_param_value("type");
        std::string last = cfg.get_param_value("type");
        if(strtol(first.c_str(), 0, 10) != 0)
        {
            firstid = strtol(first.c_str(), 0, 10);
        }
        if(strtol(last.c_str(), 0, 10) != 0)
        {
            lastid = strtol(last.c_str(), 0, 10);
        }
        grfile.append(".");
        grfile.append(first);
        grfile.append(".");
        grfile.append(std::to_string(lastid));
    }

    // compute down_distance
    warthog::bbaf_filter filter(&g, &part);
    filter.compute(firstid, lastid);
    
    // save the result
    std::cerr << "saving contracted graph to file " << grfile << std::endl;
    std::fstream out(grfile.c_str(), std::ios_base::out | std::ios_base::trunc);
    if(!out.good())
    {
        std::cerr << "\nerror exporting ch to file " << grfile << std::endl;
    }
    filter.print(out);
    out.close();
    std::cerr << "all done!\n";
}

void
compute_chbbaf_labels()
{
    std::string grfile = cfg.get_param_value("input");
    std::string cofile = cfg.get_param_value("input");
    std::string partfile = cfg.get_param_value("input");
    std::string orderfile = cfg.get_param_value("input");

    if( grfile == "" || cofile == "" || partfile == "" || orderfile == "")
    {
        std::cerr << "err; insufficient input parameters."
                  << " required, in order:\n"
                  << " --input [gr file] [co file]"
                  << " [graph partition file]"
                  << " [node ordering file]\n";
        return;
    }
    std::cerr << "computing labels" << std::endl;

    // load up the graph
    warthog::graph::planar_graph g;
    if(!g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true))
    {
        std::cerr << "err; could not load gr or co input files (one or both)\n";
        return;
    }

    // load up the partition info
    std::vector<uint32_t> part;
    if(!warthog::helpers::load_integer_labels_dimacs(partfile.c_str(), part))
    {
        std::cerr << "err; could not load partition file\n"; 
        return;
    }

    // load up the node ordering
    std::vector<uint32_t> order;
    if(!warthog::ch::load_node_order(orderfile.c_str(), order, true))
    {
        std::cerr << "err; could not load node order input file\n";
        return;
    }

    // compute labels
    uint32_t firstid = 0;
    uint32_t lastid = g.get_num_nodes()-1;
    grfile.append(".ch-bbaf.arclabel");
    if(cfg.get_num_values("type") == 2)
    {
        std::string first = cfg.get_param_value("type");
        std::string last = cfg.get_param_value("type");
        if(strtol(first.c_str(), 0, 10) != 0)
        {
            firstid = strtol(first.c_str(), 0, 10);
        }
        if(strtol(last.c_str(), 0, 10) != 0)
        {
            lastid = strtol(last.c_str(), 0, 10);
        }
        grfile.append(".");
        grfile.append(first);
        grfile.append(".");
        grfile.append(std::to_string(lastid));
    }
    warthog::bbaf_filter filter(&g, &part);
    filter.compute_ch(firstid, lastid, &order);
    
    // save the result
    std::cerr << "saving contracted graph to file " << grfile << std::endl;
    std::fstream out(grfile.c_str(), std::ios_base::out | std::ios_base::trunc);
    if(!out.good())
    {
        std::cerr << "\nerror exporting ch to file " << grfile << std::endl;
    }
    filter.print(out);
    out.close();
    std::cerr << "all done!\n";
}

void
compute_afh_labels()
{
    std::string grfile = cfg.get_param_value("input");
    std::string cofile = cfg.get_param_value("input");
    std::string partfile = cfg.get_param_value("input");
    std::string orderfile = cfg.get_param_value("input");

    if( grfile == "" || cofile == "" || partfile == "" || orderfile == "")
    {
        std::cerr << "err; insufficient input parameters."
                  << " required, in order:\n"
                  << " --input [gr file] [co file]"
                  << " [graph partition file]"
                  << " [node ordering file]\n";
        return;
    }
    std::cerr << "computing labels" << std::endl;

    // load up the graph
    warthog::graph::planar_graph g;
    if(!g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true))
    {
        std::cerr << "err; could not load gr or co input files (one or both)\n";
        return;
    }

    // load up the partition info
    std::vector<uint32_t> part;
    if(!warthog::helpers::load_integer_labels_dimacs(partfile.c_str(), part))
    {
        std::cerr << "err; could not load partition file\n"; 
        return;
    }

    // load up the node ordering
    std::vector<uint32_t> order;
    if(!warthog::ch::load_node_order(orderfile.c_str(), order, true))
    {
        std::cerr << "err; could not load node order input file\n";
        return;
    }

    // output filename reflects the type and scope of the labeling
    grfile.append(".afh.arclabel");

    // compute labels
    warthog::afh_filter filter(&g, &part);
    filter.compute(&order);
    
    // save the result
    std::cerr << "saving contracted graph to file " << grfile << std::endl;
    std::fstream out(grfile.c_str(), std::ios_base::out | std::ios_base::trunc);
    if(!out.good())
    {
        std::cerr << "\nerror exporting ch to file " << grfile << std::endl;
    }
    filter.print(out);
    out.close();
    std::cerr << "all done!\n";
}

void
compute_afhd_labels()
{
    std::string grfile = cfg.get_param_value("input");
    std::string cofile = cfg.get_param_value("input");
    std::string partfile = cfg.get_param_value("input");
    std::string orderfile = cfg.get_param_value("input");

    if( grfile == "" || cofile == "" || partfile == "" || orderfile == "")
    {
        std::cerr << "err; insufficient input parameters."
                  << " required, in order:\n"
                  << " --input [gr file] [co file]"
                  << " [graph partition file]"
                  << " [node ordering file]\n";
        return;
    }
    std::cerr << "computing labels" << std::endl;

    // load up the graph
    warthog::graph::planar_graph g;
    if(!g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true))
    {
        std::cerr << "err; could not load gr or co input files (one or both)\n";
        return;
    }

    // load up the partition info
    std::vector<uint32_t> part;
    if(!warthog::helpers::load_integer_labels_dimacs(partfile.c_str(), part))
    {
        std::cerr << "err; could not load partition file\n"; 
        return;
    }

    // load up the node ordering
    std::vector<uint32_t> order;
    if(!warthog::ch::load_node_order(orderfile.c_str(), order, true))
    {
        std::cerr << "err; could not load node order input file\n";
        return;
    }

    // output filename reflects the type and scope of the labeling
    grfile.append(".afhd.arclabel");

    // compute down_distance
    warthog::afhd_filter filter(&g, &part);
    filter.compute(&order);
    
    // save the result
    std::cerr << "saving contracted graph to file " << grfile << std::endl;
    std::fstream out(grfile.c_str(), std::ios_base::out | std::ios_base::trunc);
    if(!out.good())
    {
        std::cerr << "\nerror exporting ch to file " << grfile << std::endl;
    }
    filter.print(out);
    out.close();
    std::cerr << "all done!\n";
}

int main(int argc, char** argv)
{

	// parse arguments
    int print_help=false;
	warthog::util::param valid_args[] = 
	{
		{"help", no_argument, &print_help, 1},
		{"verbose", no_argument, &verbose, 1},
		{"input",  required_argument, 0, 2},
        {"type", required_argument, 0, 1}
	};
	cfg.parse_args(argc, argv, "-hvd:o:p:a:", valid_args);

    if(argc == 1 || print_help)
    {
		help();
        exit(0);
    }

    std::string arclabel = cfg.get_param_value("type");
    if(arclabel.compare("downdist") == 0)
    {
        compute_down_distance();
    }
    else if(arclabel.compare("chbb") == 0)
    {
        compute_chbb_labels();
    }
    else if(arclabel.compare("bb") == 0)
    {
        compute_bb_labels();
    }
    else if(arclabel.compare("dcl") == 0)
    {
        compute_dcl_labels();
    }
    else if(arclabel.compare("chaf") == 0)
    {
        compute_chaf_labels();
    }
    else if(arclabel.compare("af") == 0)
    {
        compute_af_labels();
    }
    else if(arclabel.compare("afh") == 0)
    {
        compute_afh_labels();
    }
    else if(arclabel.compare("afhd") == 0)
    {
        compute_afhd_labels();
    }
    else if(arclabel.compare("bbaf") == 0)
    {
        compute_bbaf_labels();
    }
    else if(arclabel.compare("chbbaf") == 0)
    {
        compute_chbbaf_labels();
    }
    else
    {
        std::cerr << "invalid or missing argument: --type";
    }
    return 0;
}
