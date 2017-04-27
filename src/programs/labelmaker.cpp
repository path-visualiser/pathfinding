#include "afh_filter.h"
#include "afhd_filter.h"
#include "af_labelling.h"
#include "bb_labelling.h"
#include "bbaf_labelling.h"
#include "cfg.h"
#include "ch_expansion_policy.h"
#include "corner_point_graph.h"
#include "dimacs_parser.h"
#include "dcl_filter.h"
#include "down_distance_filter.h"
#include "fch_expansion_policy.h"
#include "fch_jpg_expansion_policy.h"
#include "graph_expansion_policy.h"
#include "graph.h"
#include "gridmap.h"
#include "helpers.h"
#include "planar_graph.h"

#include <iostream>
#include <memory>
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
    << "\t--type [ downdist | dcl | af | bb | bbaf | chaf | chbb | chbbaf "
    << "| chbb-jpg | afh | afhd ]\n"
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
    

    // gogogo
    std::cerr << "computing bounding box labelling... \n";
    std::function<warthog::ch_expansion_policy*(void)> fn_new_expander = 
        [&g, &order]() -> warthog::ch_expansion_policy*
        {
            return new warthog::ch_expansion_policy(
                    &g, &order, false, warthog::ch::DOWN);
        };
    warthog::label::bb_labelling* labelling = 
        warthog::label::bb_labelling::compute<warthog::ch_expansion_policy>(
                &g, fn_new_expander);

    // save the result
    std::string outfile(grfile);
    outfile.append(".ch-bb.arclabel");
    std::cerr << "\ndone; \nsaving to " << outfile << "\n";
    std::fstream fs_out(outfile.c_str(), 
                        std::ios_base::out | std::ios_base::trunc);
    if(!fs_out.good())
    {
        std::cerr << "\nerror opening output file " << outfile << std::endl;
    }
    else
    {
        labelling->print(fs_out);
    }
    fs_out.close();

    delete labelling;
    std::cerr << "done!\n";
}

void
compute_chbb_jpg_labels()
{
    std::string grfile = cfg.get_param_value("input");
    std::string cofile = cfg.get_param_value("input");
    std::string orderfile = cfg.get_param_value("input");
    std::string gridmapfile = cfg.get_param_value("input");

    if( grfile == "" || cofile == "" || orderfile == "" || gridmapfile == "")
    {
        std::cerr << "err; insufficient input parameters."
                  << " required, in order:\n"
                  << " --input [gr file] [co file]"
                  << " [node ordering file] [gridmap]\n";
        return;
    }
    std::cerr << "computing labels" << std::endl;

    // load up the grid
    std::shared_ptr<warthog::gridmap> map(
            new warthog::gridmap(gridmapfile.c_str()));

    // load up the graph 
    std::shared_ptr<warthog::graph::planar_graph> pg(
            new warthog::graph::planar_graph());
    if(!pg->load_dimacs(grfile.c_str(), cofile.c_str(), false, true))
    {
        std::cerr << "err; could not load gr or co input files " 
                  << "(one or both)\n";
        return;
    }
    warthog::graph::corner_point_graph cpg(map, pg);

    // load up the node ordering
    std::vector<uint32_t> order;
    if(!warthog::ch::load_node_order(orderfile.c_str(), order, true))
    {
        std::cerr << "err; could not load node order input file\n";
        return;
    }
    
    std::cerr << "creating fch-jpg-bb labelling\n";
    std::function<warthog::fch_jpg_expansion_policy*(void)> 
        fn_new_expander = [&cpg, &order] () 
            -> warthog::fch_jpg_expansion_policy*
        {
            return new warthog::fch_jpg_expansion_policy(&cpg, &order);
        };

    warthog::label::bb_labelling* labelling = 
    warthog::label::bb_labelling::compute<warthog::fch_jpg_expansion_policy>(
            pg.get(), fn_new_expander);

    std::cerr << "labelling done\n";

    // save the result
    std::string outfile(grfile);
    outfile.append(".fch-bb-jpg.arclabel");
    std::fstream fs_out(outfile.c_str(), 
                        std::ios_base::out | std::ios_base::trunc);

    std::cerr << "saving to" << outfile << "\n";
    if(!fs_out.good())
    {
        std::cerr << "\nerror opening output file " << outfile << std::endl;
    }
    else
    {
        labelling->print(fs_out);
    }

    fs_out.close();
    delete labelling;
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

    // gogogo
    std::cerr << "computing bounding box labelling... \n";
    std::function<warthog::graph_expansion_policy*(void)> fn_new_expander = 
        [&g]() -> warthog::graph_expansion_policy*
        {
            return new warthog::graph_expansion_policy(&g);
        };
    warthog::label::bb_labelling* labelling = 
        warthog::label::bb_labelling::compute
            <warthog::graph_expansion_policy>
                (&g, fn_new_expander);

    // save the result
    std::string outfile(grfile);
    outfile.append(".bb.arclabel");
    std::cerr << "\ndone; \nsaving to " << outfile << "\n";
    std::fstream fs_out(outfile.c_str(), 
                        std::ios_base::out | std::ios_base::trunc);
    if(!fs_out.good())
    {
        std::cerr << "\nerror opening output file " << outfile << std::endl;
    }
    else
    {
        labelling->print(fs_out);
    }
    fs_out.close();

    delete labelling;
    std::cerr << "done!\n";
}

void
compute_chaf_labels()
{
    std::string grfile = cfg.get_param_value("input");
    std::string cofile = cfg.get_param_value("input");
    std::string orderfile = cfg.get_param_value("input");
    std::string partfile = cfg.get_param_value("input");

    if( grfile == "" || cofile == "" || partfile == "" || orderfile == "")
    {
        std::cerr << "err; insufficient input parameters."
                  << " required, in order:\n"
                  << " --input [gr file] [co file]"
                  << " [node ordering file]"
                  << " [graph partition file]\n";
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

    uint32_t firstid = 0;
    uint32_t lastid = g.get_num_nodes()-1;
    std::string outfile(grfile);
    outfile.append(".fch-af.arclabel");
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
        outfile.append(".");
        outfile.append(first);
        outfile.append(".");
        outfile.append(std::to_string(lastid));
    }

    // gogogo
    std::cerr << "computing fch-af labelling... \n";
    std::function<warthog::fch_expansion_policy*(void)> fn_new_expander = 
        [&g, &order]() -> warthog::fch_expansion_policy*
        {
            return new warthog::fch_expansion_policy(&g, &order);
        };
    warthog::label::af_labelling* labelling = 
        warthog::label::af_labelling::compute
            <warthog::fch_expansion_policy>
                (&g, &part, fn_new_expander, firstid, lastid);

    // save the result
    std::cerr << "\ndone; \nsaving to " << outfile << "\n";
    std::fstream fs_out(outfile.c_str(), 
                        std::ios_base::out | std::ios_base::trunc);
    if(!fs_out.good())
    {
        std::cerr << "\nerror opening output file " << outfile << std::endl;
    }
    else
    {
        labelling->print(fs_out, firstid, lastid);
    }
    fs_out.close();

    delete labelling;
    std::cerr << "done!\n";
}

void
compute_chaf_jpg_labels()
{
    std::string grfile = cfg.get_param_value("input");
    std::string cofile = cfg.get_param_value("input");
    std::string orderfile = cfg.get_param_value("input");
    std::string partfile = cfg.get_param_value("input");
    std::string gridmapfile = cfg.get_param_value("input");

    if( grfile == "" || cofile == "" || partfile == "" ||
            orderfile == "" || gridmapfile == "")
    {
        std::cerr << "err; insufficient input parameters."
                  << " required, in order:\n"
                  << " --input [gr file] [co file]"
                  << " [node ordering file]"
                  << " [graph partition file]"
                  << " [gridmap]\n";
        return;
    }
    std::cerr << "computing labels" << std::endl;

    // load up the grid
    std::shared_ptr<warthog::gridmap> map(
            new warthog::gridmap(gridmapfile.c_str()));

    // load up the graph
    std::shared_ptr<warthog::graph::planar_graph> pg(
            new warthog::graph::planar_graph());
    if(!pg->load_dimacs(grfile.c_str(), cofile.c_str(), false, true))
    {
        std::cerr << "err; could not load gr or co input files " 
                  << "(one or both)\n";
        return;
    }
    warthog::graph::corner_point_graph cpg(map, pg);

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

    // gogogo
    std::cerr << "computing fch-af labelling... \n";
    std::function<warthog::fch_jpg_expansion_policy*(void)> fn_new_expander = 
        [&cpg, &order]() -> warthog::fch_jpg_expansion_policy*
        {
            return new warthog::fch_jpg_expansion_policy(&cpg, &order);
        };
    warthog::label::af_labelling* labelling = 
        warthog::label::af_labelling::compute
            <warthog::fch_jpg_expansion_policy>
                (pg.get(), &part, fn_new_expander);

    // save the result
    std::string outfile = grfile;
    outfile.append(".fch-af-jpg.arclabel");
    std::cerr << "saving arclabels file " << outfile << std::endl;
    std::fstream out(outfile.c_str(), 
            std::ios_base::out | std::ios_base::trunc);
    if(!out.good())
    {
        std::cerr << "\nerror saving arclabels file " << grfile << std::endl;
    }
    else
    {
        labelling->print(out);
    }
    out.close();

    // cleanup
    delete labelling;
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
    
    // gogogo
    std::cerr << "computing af labelling... \n";
    std::function<warthog::graph_expansion_policy*(void)> fn_new_expander = 
        [&g]() -> warthog::graph_expansion_policy*
        {
            std::cerr << "new graph_expansion_policy" << std::endl;
            return new warthog::graph_expansion_policy(&g);
        };
    warthog::label::af_labelling* labelling = 
        warthog::label::af_labelling::compute
            <warthog::graph_expansion_policy>
                (&g, &part, fn_new_expander);

    // save the result
    std::string outfile = grfile;
    outfile.append(".af.arclabel");
    std::cerr << "saving arclabels file " << outfile << std::endl;
    std::fstream out(outfile.c_str(), 
            std::ios_base::out | std::ios_base::trunc);
    if(!out.good())
    {
        std::cerr << "\nerror saving arclabels file " << grfile << std::endl;
    }
    else
    {
        labelling->print(out);
    }
    out.close();

    delete labelling;
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
    std::string outfile(grfile);
    outfile.append(".bbaf.arclabel");
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
        outfile.append(".");
        outfile.append(first);
        outfile.append(".");
        outfile.append(std::to_string(lastid));
    }

    // compute bbaf labels
    std::function<warthog::graph_expansion_policy*(void)> fn_new_expander = 
    [&g]() -> warthog::graph_expansion_policy*
    {
        return new warthog::graph_expansion_policy(&g);
    };

    std::shared_ptr<warthog::label::bbaf_labelling> 
        labelling(
            warthog::label::bbaf_labelling::compute
                <warthog::graph_expansion_policy>
                    (&g, &part, fn_new_expander, firstid, lastid));

    // save the result
    std::cerr << "\ndone; \nsaving to " << outfile << "\n";
    std::fstream fs_out(outfile.c_str(), 
                        std::ios_base::out | std::ios_base::trunc);
    if(!fs_out.good())
    {
        std::cerr << "\nerror opening output file " << outfile << std::endl;
    }
    else
    {
        labelling->print(fs_out, firstid, lastid);
    }
    fs_out.close();
    std::cerr << "done!\n";
}

void
compute_chbbaf_labels()
{
    std::string grfile = cfg.get_param_value("input");
    std::string cofile = cfg.get_param_value("input");
    std::string orderfile = cfg.get_param_value("input");
    std::string partfile = cfg.get_param_value("input");

    if( grfile == "" || cofile == "" || partfile == "" || orderfile == "")
    {
        std::cerr << "err; insufficient input parameters."
                  << " required, in order:\n"
                  << " --input [gr file] [co file]"
                  << " [node ordering file]"
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
    std::string outfile(grfile);
    outfile.append(".fch-bbaf.arclabel");
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
        outfile.append(".");
        outfile.append(first);
        outfile.append(".");
        outfile.append(std::to_string(lastid));
    }

    std::function<warthog::fch_expansion_policy*(void)> fn_new_expander = 
    [&g, &order]() -> warthog::fch_expansion_policy*
    {
        return new warthog::fch_expansion_policy(&g, &order);
    };

    std::shared_ptr<warthog::label::bbaf_labelling> 
        labelling(
            warthog::label::bbaf_labelling::compute
                <warthog::fch_expansion_policy>
                    (&g, &part, fn_new_expander, firstid, lastid));
    
    // save the result
    std::cerr << "\ndone; \nsaving to " << outfile << "\n";
    std::fstream fs_out(outfile.c_str(), 
                        std::ios_base::out | std::ios_base::trunc);
    if(!fs_out.good())
    {
        std::cerr << "\nerror opening output file " << outfile << std::endl;
    }
    else
    {
        labelling->print(fs_out, firstid, lastid);
    }
    fs_out.close();
    std::cerr << "done!\n";
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
    else if(arclabel.compare("chbb-jpg") == 0)
    {
        compute_chbb_jpg_labels();
    }
    else if(arclabel.compare("chaf-jpg") == 0)
    {
        compute_chaf_jpg_labels();
    }
    else
    {
        std::cerr << "invalid or missing argument: --type\n";
    }
    return 0;
}
