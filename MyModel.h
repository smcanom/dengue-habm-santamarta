#ifndef MYMODEL
#define MYMODEL

#include <fstream>
#include <boost/mpi.hpp>
#include <chrono> 
#include <unordered_map>
#include "repast_hpc/Schedule.h"
#include "repast_hpc/Properties.h"
#include "repast_hpc/SharedContext.h"
#include "repast_hpc/SharedDiscreteSpace.h"
#include "repast_hpc/GridComponents.h"
#include "repast_hpc/ValueLayerND.h"
#include "repast_hpc/AgentRequest.h"
#include "repast_hpc/TDataSource.h"
#include "repast_hpc/SVDataSet.h"
#include "repast_hpc/Point.h"
#include "MyHuman.h"
#include "CSVWriter.h"
#include "SEIModel.h"
#include <mutex>

/* Agent Package Provider (NO ESTOY MUY SEGURA DE PARA QUE SIRVE) */
class HumanPackageProvider {
	private:
    	repast::SharedContext<Human>* agents;
	
	public:
    	HumanPackageProvider(repast::SharedContext<Human>* agentPtr);
    	void providePackage(Human * agent, std::vector<HumanPackage>& out);
    	void provideContent(repast::AgentRequest req, std::vector<HumanPackage>& out);
	
};

/* Agent Package Receiver (NO ESTOY MUY SEGURA DE PARA QUE SIRVE) */
class HumanPackageReceiver {
	private:
		repast::SharedContext<Human>* agents;
		
	public:	
		HumanPackageReceiver(repast::SharedContext<Human>* agentPtr);
		Human * createAgent(HumanPackage package);
		void updateAgent(HumanPackage package);
	
};


class RepastHPCModel {
private:
    // Model parameters from props
	//std::vector<std::vector<std::pair<double, double>>> geojson;
    // each element is one (x_centroid, y_centroid)
    std::vector<std::pair<double,double>> patchCenters;

    int stopAt;
    int countOfHumans;
    int countOfInfectedHumans;
    int xdim;
    int ydim;

    CSVWriter writer;
    repast::Properties* props;
    repast::SharedContext<Human> context;
    repast::SharedDiscreteSpace<Human, repast::StrictBorders, repast::SimpleAdder<Human>>* discreteSpace;

    // Spatial layers
    repast::ValueLayerND<int>* valueLayerSuceptibleMosquitoes;
    repast::ValueLayerND<int>* valueLayerExposedMosquitoes;
    repast::ValueLayerND<int>* valueLayerInfectedMosquitoes;
    repast::ValueLayerND<int>* valueLayerTemperature;
    repast::ValueLayerND<int>* valueLayerType;
    repast::ValueLayerND<int>* valuetotalmosquitoes0;
    repast::ValueLayerND<int>* valuetotalmosquitoes1;

    // Temperature time series
    std::vector<std::vector<int>> dataTemperatures;
    // persistent SEIModel grid and fast lookup
    std::vector<SEIModel> mosquitoPatches;
    std::unordered_map<long long, std::pair<int,int>> locToIndex;
	std::vector<std::vector<int>> residentialLocations;
	std::vector<std::vector<int>> studyLocations;
	std::vector<std::vector<int>> workLocations;
	std::vector<std::vector<int>> otherActivitiesLocations;
    std::vector<std::pair<int, int>> border_points;

    std::chrono::time_point<std::chrono::steady_clock> start_time;  // Add this
    std::chrono::time_point<std::chrono::steady_clock> end_time;    // Add this

    // Agent communication helpers
    HumanPackageProvider* provider;
    HumanPackageReceiver* receiver;

    // Mutex for thread safety
    std::mutex patchMutex;
    std::mutex valueLayerMutex;

public:
    RepastHPCModel(std::string propsFile, int argc, char** argv, boost::mpi::communicator* comm);
    ~RepastHPCModel();
    
    void init();
    void initHumans();
    void initGridValueLayers();
    void writeCsv(const std::string& filename);
    void runAllHumans();
    void runAllPatches();
    void initRecords();
    void recordResults();
    void initSchedule(repast::ScheduleRunner& runner);
    double get_border_x(int coordY);
    void printExecutionTime(); 

    // helpers
    int ageInitializer(double prob);
    int* countHumansInState();
    std::vector<int> homeLocationInitializer();
    std::vector<std::vector<int>> activitiesInitializer(int age);
    std::vector<int> placeSelector(int location_value);
};

#endif // MYMODEL_H

