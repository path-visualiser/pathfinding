#include "dimacs_parser.h"

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <ostream>
#include <string>
#include <vector>

warthog::dimacs_parser::dimacs_parser()
{
    init();
}

warthog::dimacs_parser::dimacs_parser(const char* gr_file)
{
    init();
    load(gr_file);
}

warthog::dimacs_parser::dimacs_parser(const char* co_file, 
        const char* gr_file)
{
    init();
    load(co_file);
    load(gr_file);
}

warthog::dimacs_parser::~dimacs_parser()
{
    delete nodes_;
    delete edges_;
}

void
warthog::dimacs_parser::init()
{
    nodes_ = new std::vector<warthog::dimacs_parser::node>();
    edges_ = new std::vector<warthog::dimacs_parser::edge>();
}

bool
warthog::dimacs_parser::load(const char* filename)
{
    std::fstream* fdimacs = new std::fstream(filename, std::fstream::in);
	if(!fdimacs->is_open())
	{
		std::cerr << "err; dimacs_parser::dimacs_parser "
			"cannot open file: "<<filename << std::endl;
		exit(1);
	}

	uint32_t retval = 1;
	char* buf = new char[1024];
	const char* delim = " \t";
	uint32_t line = 0;
	while(fdimacs->good())
	{
		fdimacs->getline(buf, 1024);
		if(buf[0] == 'p')
		{
			char* type = strtok(buf, delim); // p char
			type = strtok(NULL, delim);
			if(!strcmp(type, "sp"))
			{
				retval = load_gr_file(*fdimacs);
			}
			else if(!strcmp(type, "aux"))
			{
				retval = load_co_file(*fdimacs);
			}
			else
			{
				std::cerr << "error; unrecognised problem line in dimacs file\n";
				retval = 1;
				break;
			}
		}
		line++;
	}

    delete fdimacs;
	delete buf;
	return retval;
}

bool
warthog::dimacs_parser::load_co_file(std::istream& fdimacs)
{
    nodes_->clear();

    uint32_t line = 1;
	char* buf = new char[1024];
	const char* delim = " \t";
	bool early_abort = false;
	while(fdimacs.good() && !early_abort)
	{
		char next_char = fdimacs.peek();
		switch(next_char)
		{
			case 'v':
			{
				fdimacs.getline(buf, 1024);
				char* descriptor = strtok(buf, delim);
				char* id = strtok(NULL, delim);
				char* x = strtok(NULL, delim);
				char* y = strtok(NULL, delim);
				if(!(descriptor && id && x && y))
				{
					std::cerr << "warning; badly formatted node descriptor on line "
                        <<line << std::endl;
					break;
				}
                warthog::dimacs_parser::node n;
                n.id_ = atoi(id);
                n.x_ = atoi(x);
                n.y_ = atoi(y);
                nodes_->push_back(n);

				break;
			}
			case 'p': // stop if we hit another problem line
				early_abort = true;
				break;
			default: // ignore non-node, non-problem lines
				fdimacs.getline(buf, 1024);
				break;
		}
		line++;
	}

	delete buf;
	return 0;

}

bool
warthog::dimacs_parser::load_gr_file(std::istream& fdimacs)
{
    edges_->clear();
    this->n_edges_ = 0;

    uint32_t line = 1;
	char* buf = new char[1024];
	const char* delim = " \t";
	bool early_abort = false;
	while(fdimacs.good() && !early_abort)
	{
		char next_char = fdimacs.peek();
		switch(next_char)
		{
			case 'a':
			{
				fdimacs.getline(buf, 1024);
				char* descriptor = strtok(buf, delim);
				char* from = strtok(NULL, delim);
				char* to = strtok(NULL, delim);
				char* cost = strtok(NULL, delim);
				if(!(descriptor && from && to && cost))
				{
					std::cerr << "warning; badly formatted arc descriptor on line "
                        <<line << std::endl;
					break;
				}
                warthog::dimacs_parser::edge e;
                e.tail_id_ = atoi(from);
                e.head_id_ = atoi(to);
                e.weight_ = atoi(cost);
                edges_->push_back(e);
				break;
			}
			case 'p': // another problem line. stop here
				early_abort = true;
				break;
			default: // ignore non-arc, non-problem lines
				fdimacs.getline(buf, 1024);
				break;
		}
		line++;
	}

	delete buf;
	return 0;
}

void
warthog::dimacs_parser::print(std::ostream& oss)
{
    uint32_t nnodes = nodes_->size();
    if(nnodes > 0)
    {
        oss << "p aux sp co " << nodes_->size() << std::endl;
        for(uint32_t i = 0; i < nnodes; i++)
        {
            warthog::dimacs_parser::node n = nodes_->at(i);
            oss << "v " << i+1 << " " << n.x_ << " " << n.y_ << std::endl;
        }
    }

    uint32_t nedges = edges_->size();
    if(nedges > 0)
    {
        oss << "p sp " << nnodes << " " << nedges << std::endl;
        for(uint32_t i = 0; i < nedges; i++)
        {
            warthog::dimacs_parser::edge e = edges_->at(i);
            oss << "a " << e.tail_id_ << " " << e.head_id_ << " " << e.weight_ << std::endl;
        }
    }

}
