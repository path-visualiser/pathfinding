#ifndef SCENARIOMANAGER_H
#define SCENARIOMANAGER_H

#include "Experiment.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace ScenarioManagerNS
{
	const int MAX_CONSECUTIVE_FAILS=1000;
}

using namespace ScenarioManagerNS;

class TooManyTriesException : public std::exception
{
	public:
		TooManyTriesException(int _generated, int _target) : generated(_generated), target(_target) { }
		const char* what() const throw()
		{
			std::stringstream errmsg;
			errmsg << "\n well, this sucks. can't generate enough experiments for the given map. ";
			errmsg << "I managed to create " << generated << " of "<<target;
			return errmsg.str().c_str();
		}
	
	private:
		int generated, target;
};

class mapAbstraction;

/*
 * A note on .scenario versions:
 *
 * 	v2.0: Generated by AHAScenarioManager. order: sx,sy,gx,gy,capability,size,distance,map
 * 	v2.1: Generated by AHAScenarioManager. order: sx,sy,gx,gy,capability,size,distance,map
 * 	v3.0: Generated by ScenarioManager: order: sx,sy,gx,gy,distance,map
 */
class AbstractScenarioManager 
{
	public:
		AbstractScenarioManager(){};
		virtual ~AbstractScenarioManager();

		Experiment* getNthExperiment(int which) { if(which < (int)experiments.size()) return experiments[which]; return 0; }
		void addExperiment(Experiment* newexp) { experiments.push_back(newexp); }
		int getNumExperiments() { return experiments.size(); }
		
		virtual void generateExperiments(mapAbstraction* absMap, int numexperiments) 
			throw(TooManyTriesException) = 0;
		virtual void loadScenarioFile(const char* filelocation) 
			 throw(std::invalid_argument) = 0;
		void writeScenarioFile(const char* filelocation);
		void clearExperiments() { experiments.clear(); }
		void sortExperiments(); // organise by increasing solution length
	
	protected: 
		std::vector<Experiment*> experiments;		
		int version;
};

class ScenarioManager : public AbstractScenarioManager
{
	public:
		ScenarioManager();
		virtual ~ScenarioManager();

		virtual void 
		generateExperiments(mapAbstraction* absMap, int numexperiments)
		throw(TooManyTriesException);

		virtual void loadScenarioFile(const char* filelocation)
			 throw(std::invalid_argument);

	protected:
		Experiment* generateSingleExperiment(mapAbstraction* absMap);
		void loadV1ScenarioFile(std::ifstream& infile);
		void loadV3ScenarioFile(std::ifstream& infile);
};

#endif
