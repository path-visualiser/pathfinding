#include "apex_distance_filter.h"
#include "af_filter.h"
#include "afh_filter.h"
#include "afhd_filter.h"
#include "bb_af_filter.h"
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
    << "\t--dimacs [gr file] [co file] (IN THIS ORDER!!)\n"
	<< "\t--order [order-of-contraction file]\n"
    << "\t--type [downdist | dcl | af | bb | bbaf | chaf | chbb | afh | afhd ]\n"
	<< "\t--verbose (optional)\n";
}

void 
compute_down_distance()
{
    std::string grfile = cfg.get_param_value("dimacs");
    std::string cofile = cfg.get_param_value("dimacs");
    std::cerr << "param values " << std::endl;
    std::string orderfile = cfg.get_param_value("order");
    cfg.print_values("dimacs", std::cerr);

    std::cerr << "grfile: "<< grfile << " cofile " << cofile << std::endl;

    // load up (or create) the contraction hierarchy
    warthog::graph::planar_graph g;
    std::vector<uint32_t> order;

    if(orderfile.compare("") != 0)
    {
        g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true);
        warthog::ch::load_node_order(orderfile.c_str(), order, true);

        // compute down_distance
        grfile.append(".ddist.arclabel");
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
        warthog::down_distance_filter ddfilter(&g, &order);
        ddfilter.compute(firstid, lastid);
        
        // save the result
        std::cerr << "saving contracted graph to file " << grfile << std::endl;
        std::fstream out(grfile.c_str(), std::ios_base::out | std::ios_base::trunc);
        if(!out.good())
        {
            std::cerr << "\nerror exporting ch to file " << grfile << std::endl;
        }
        ddfilter.print(out);
        out.close();
    }
    else
    {
        std::cerr << "required: node order file. aborting.\n";
        return;
    }
    std::cerr << "all done!\n";
}

void 
compute_dcl_labels()
{
    std::string grfile = cfg.get_param_value("dimacs");
    std::string cofile = cfg.get_param_value("dimacs");
    std::cerr << "param values " << std::endl;
    std::string orderfile = cfg.get_param_value("order");
    cfg.print_values("dimacs", std::cerr);

    std::cerr << "grfile: "<< grfile << " cofile " << cofile << std::endl;

    // load up (or create) the contraction hierarchy
    warthog::graph::planar_graph g;
    std::vector<uint32_t> order;

    if(orderfile.compare("") != 0)
    {
        g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true);
        warthog::ch::load_node_order(orderfile.c_str(), order, true);

        // compute down closure
        grfile.append(".dcl.arclabel");
        warthog::dcl_filter filter(&g, &order);
        filter.compute();
        
        // save the result
        std::cerr << "saving contracted graph to file " << grfile << std::endl;
        std::fstream out(grfile.c_str(), std::ios_base::out | std::ios_base::trunc);
        if(!out.good())
        {
            std::cerr << "\nerror exporting ch to file " << grfile << std::endl;
        }
        filter.print(out);
        out.close();
    }
    else
    {
        std::cerr << "required: node order file. aborting.\n";
        return;
    }
    std::cerr << "all done!\n";
}

void 
compute_apex_distance()
{
    std::string grfile = cfg.get_param_value("dimacs");
    std::string cofile = cfg.get_param_value("dimacs");
    std::cerr << "param values " << std::endl;
    std::string orderfile = cfg.get_param_value("order");
    cfg.print_values("dimacs", std::cerr);

    std::cerr << "grfile: "<< grfile << " cofile " << cofile << std::endl;

    // load up (or create) the contraction hierarchy
    warthog::graph::planar_graph g;
    std::vector<uint32_t> order;

    if(orderfile.compare("") != 0)
    {
        g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true);
        warthog::ch::load_node_order(orderfile.c_str(), order, true);

        // compute down_distance
        grfile.append(".apex.arclabel");
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
            grfile.append(std::to_string(lastid-1));
        }
        warthog::apex_distance_filter apexfilter(&g, &order);
        apexfilter.compute(firstid, lastid);
        
        // save the result
        std::cerr << "saving contracted graph to file " << grfile << std::endl;
        std::fstream out(grfile.c_str(), std::ios_base::out | std::ios_base::trunc);
        if(!out.good())
        {
            std::cerr << "\nerror exporting ch to file " << grfile << std::endl;
        }
        apexfilter.print(out);
        out.close();
    }
    else
    {
        std::cerr << "required: node order file. aborting.\n";
        return;
    }
    std::cerr << "all done!\n";
}

void
compute_chbb_labels()
{
    std::string grfile = cfg.get_param_value("dimacs");
    std::string cofile = cfg.get_param_value("dimacs");
    std::cerr << "param values " << std::endl;
    std::string orderfile = cfg.get_param_value("order");
    cfg.print_values("dimacs", std::cerr);

    std::cerr << "grfile: "<< grfile << " cofile " << cofile << std::endl;

    // load up (or create) the contraction hierarchy
    warthog::graph::planar_graph g;
    std::vector<uint32_t> order;

    if(orderfile.compare("") != 0)
    {
        g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true);
        warthog::ch::load_node_order(orderfile.c_str(), order, true);

        // compute down_distance
        grfile.append(".chbb.arclabel");
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
            grfile.append(std::to_string(lastid-1));
        }
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
    }
    else
    {
        std::cerr << "required: node order file. aborting.\n";
        return;
    }
    std::cerr << "all done!\n";
}

void
compute_bb_labels()
{
    std::string grfile = cfg.get_param_value("dimacs");
    std::string cofile = cfg.get_param_value("dimacs");
    std::cerr << "param values " << std::endl;
    std::string orderfile = cfg.get_param_value("order");
    cfg.print_values("dimacs", std::cerr);

    std::cerr << "grfile: "<< grfile << " cofile " << cofile << std::endl;

    // load up (or create) the contraction hierarchy
    warthog::graph::planar_graph g;
    std::vector<uint32_t> order;

    g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true);
    warthog::ch::load_node_order(orderfile.c_str(), order, true);

    // compute down_distance
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
        grfile.append(std::to_string(lastid-1));
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
    std::cerr << "computing ch arcflags " << std::endl;
    std::string grfile = cfg.get_param_value("dimacs");
    std::string cofile = cfg.get_param_value("dimacs");
    std::cerr << "param values " << std::endl;
    std::string orderfile = cfg.get_param_value("order");
    std::string partfile = cfg.get_param_value("part");
    cfg.print_values("dimacs", std::cerr);

    std::cerr << "grfile: "<< grfile << " cofile " << cofile << std::endl;

    // load up (or create) the contraction hierarchy
    warthog::graph::planar_graph g;
    std::vector<uint32_t> order;
    std::vector<uint32_t> part;

    if(orderfile.compare("") != 0 && partfile.compare("") != 0)
    {
        g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true);
        warthog::ch::load_node_order(orderfile.c_str(), order, true);
        warthog::helpers::load_integer_labels_dimacs(partfile.c_str(), part);

        // output filename reflects the type and scope of the labeling
        grfile.append(".chaf.arclabel");
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

        // compute down_distance
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
    }
    else
    {
        std::cerr << "required params: --order [node order file] --part [graph partition]\n";
        return;
    }
    std::cerr << "all done!\n";
}

void
compute_af_labels()
{
    std::string grfile = cfg.get_param_value("dimacs");
    std::string cofile = cfg.get_param_value("dimacs");
    std::string partfile = cfg.get_param_value("part");

    if( grfile == "" || cofile == "" || partfile == "")
    {
        std::cerr << "err; insufficient input parameters. required:\n "
                  << " --dimacs [gr file] [co file] --part [graph partition file]\n";
        return;
    }
    std::cerr << "computing arcflags" << std::endl;


    // load up the graph
    warthog::graph::planar_graph g;
    g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true);

    // load up the partition info
    std::vector<uint32_t> part;
    warthog::helpers::load_integer_labels_dimacs(partfile.c_str(), part);

    // output filename reflects the type and scope of the labeling
    grfile.append(".af.arclabel");
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
    std::cerr << "computing bbaf labels" << std::endl;
    std::string grfile = cfg.get_param_value("dimacs");
    std::string cofile = cfg.get_param_value("dimacs");
    std::cerr << "param values " << std::endl;
    std::string orderfile = cfg.get_param_value("order");
    std::string partfile = cfg.get_param_value("part");
    cfg.print_values("dimacs", std::cerr);

    std::cerr << "grfile: "<< grfile << " cofile " << cofile << std::endl;

    // load up (or create) the contraction hierarchy
    warthog::graph::planar_graph g;
    std::vector<uint32_t> order;
    std::vector<uint32_t> part;

    if(orderfile.compare("") != 0)
    {
        g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true);
        warthog::ch::load_node_order(orderfile.c_str(), order, true);
        warthog::helpers::load_integer_labels_dimacs(partfile.c_str(), part);

        if(part.size() == 0) 
        { 
            std::cerr << "err; could not load partition file\n"; 
            return;
        }

        // output filename reflects the type and scope of the labeling
        grfile.append(".bbaf.arclabel");
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

        // compute down_distance
        warthog::bb_af_filter filter(&g, &order, &part);
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
    }
    else
    {
        std::cerr << "required: node order file. aborting.\n";
        return;
    }
    std::cerr << "all done!\n";
}

void
compute_afh_labels()
{
    std::cerr << "computing hiearchical arcflags" << std::endl;
    std::string grfile = cfg.get_param_value("dimacs");
    std::string cofile = cfg.get_param_value("dimacs");
    std::cerr << "param values " << std::endl;
    std::string orderfile = cfg.get_param_value("order");
    std::string partfile = cfg.get_param_value("part");
    cfg.print_values("dimacs", std::cerr);

    std::cerr << "grfile: "<< grfile << " cofile " << cofile << std::endl;

    // load up (or create) the contraction hierarchy
    warthog::graph::planar_graph g;
    std::vector<uint32_t> order;
    std::vector<uint32_t> part;

    if(orderfile.compare("") != 0)
    {
        g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true);
        warthog::ch::load_node_order(orderfile.c_str(), order, true);
        warthog::helpers::load_integer_labels_dimacs(partfile.c_str(), part);

        // output filename reflects the type and scope of the labeling
        grfile.append(".afh.arclabel");
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

        // compute down_distance
        warthog::afh_filter filter(&g, &order, &part);
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
    }
    else
    {
        std::cerr << "required: node order file. aborting.\n";
        return;
    }
    std::cerr << "all done!\n";
}

void
compute_afhd_labels()
{
    std::cerr << "computing afhd labels " << std::endl;
    std::string grfile = cfg.get_param_value("dimacs");
    std::string cofile = cfg.get_param_value("dimacs");
    std::cerr << "param values " << std::endl;
    std::string orderfile = cfg.get_param_value("order");
    std::string partfile = cfg.get_param_value("part");
    cfg.print_values("dimacs", std::cerr);

    std::cerr << "grfile: "<< grfile << " cofile " << cofile << std::endl;

    // load up (or create) the contraction hierarchy
    warthog::graph::planar_graph g;
    std::vector<uint32_t> order;
    std::vector<uint32_t> part;

    if(orderfile.compare("") != 0)
    {
        g.load_dimacs(grfile.c_str(), cofile.c_str(), false, true);
        warthog::ch::load_node_order(orderfile.c_str(), order, true);
        warthog::helpers::load_integer_labels_dimacs(partfile.c_str(), part);

        // output filename reflects the type and scope of the labeling
        grfile.append(".afhd.arclabel");
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

        // compute down_distance
        warthog::afhd_filter filter(&g, &order, &part);
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
    }
    else
    {
        std::cerr << "required: node order file. aborting.\n";
        return;
    }
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
		{"dimacs",  required_argument, 0, 2},
		{"order",  required_argument, 0, 1},
		{"part",  required_argument, 0, 1},
        {"type", required_argument, 0, 1}
	};
	cfg.parse_args(argc, argv, "-hvd:o:p:a:", valid_args);

    if(argc == 1 || print_help)
    {
		help();
        exit(0);
    }

    if(cfg.get_num_values("dimacs") < 2)
    {
        std::cerr << "insufficient values for param --dimacs (need gr and co files)\n";
        exit(0);
    }

    std::string arclabel = cfg.get_param_value("type");

    if(arclabel.compare("downdist") == 0)
    {
        compute_down_distance();
    }
    else if(arclabel.compare("apexdist") == 0)
    {
        compute_apex_distance();
    }
    else if(arclabel.compare("chbb") == 0)
    {
        // faster precomputation when processing
        // contraction hierarchies
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
    else
    {
        std::cerr << "invalid or missing argument: --type";
    }
    return 0;
}
