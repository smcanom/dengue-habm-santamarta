#ifndef MYMODEL
#define MYMODEL

#include <fstream>
#include <boost/mpi.hpp>
#include <chrono> 
#include <unordered_map>
#include <unordered_set>
#include <mutex>
#include <random>
#include <chrono>
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

/* Parameter carrier for sensitivity analysis */
struct Params {
    double sigma_M{0.3};      // mosquito mortality sensitivity
    double sigma_H{3.0};      // human bite sensitivity
    double z{0.3};            // bite aggregation parameter
    double r{0.6};            // reproduction rate
    double C{30.0};           // carrying capacity
    double beta_mh{0.10};     // mosquito-to-human transmission
    double beta_hm{0.10};     // human-to-mosquito transmission
    int base_seed{12345};     // base random seed
    int replicate_id{0};      // replicate identifier (for CRN)
    std::string config_id{"baseline"};
    std::string perturb_param{"baseline"};
    double perturb_delta{0.0};
    std::string obs_csv{""};  // path to observed data
};

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
		const Params* params_;
		
	public:	
		HumanPackageReceiver(repast::SharedContext<Human>* agentPtr, const Params* params);
		Human * createAgent(HumanPackage package);
		void updateAgent(HumanPackage package);
	
};


class RepastHPCModel {
private:
    // Model parameters from props
    const Params params_;  // Immutable parameter set for this run
	//std::vector<std::vector<std::pair<double, double>>> geojson;
    // each element is one (x_centroid, y_centroid)
    std::vector<std::pair<double,double>> patchCenters;

    int stopAt;
    int countOfHumans;
    int countOfInfectedHumans;
    int xdim;
    int ydim;
    int run_seed;  // Computed seed for this run (base_seed + replicate_id)
    // Neighborhood case counts (per tick)
    int new_cases_oasis = 0;
    int new_cases_laquinina = 0;
    int new_cases_ondasdelcaribe = 0;
    int new_cases_pantano = 0;
    int new_cases_ochodediciembre = 0;
    int new_cases_lacoquera = 0;
    int new_cases_sanjacinto = 0;
    int new_cases_nuevabethel = 0;
    
    int total_humans_oasis = 0;
    int total_humans_laquinina = 0;
    int total_humans_ondasdelcaribe = 0;
    int total_humans_pantano = 0;
    int total_humans_ochodediciembre = 0;
    int total_humans_lacoquera = 0;
    int total_humans_sanjacinto = 0;
    int total_humans_nuevabethel = 0;


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
    // Neighborhood tagging layer: 0=None, 1=Oasis, 2=LaQuinina, 3=OndasDelCaribe,
    // 4=Pantano, 5=OchoDeDiciembre, 6=LaCoquera, 7=SanJacinto, 8=NuevaBethel
    repast::ValueLayerND<int>* valueLayerNeighborhood = nullptr;
    // HI activation mask layer (stores X_i ∈ {0,1} for demonstration purposes)
    repast::ValueLayerND<int>* valueLayerHI_activated = nullptr;
    // Membership sets for each neighborhood (packed x,y keys)
    std::unordered_set<long long> oasis_cells;
    std::unordered_set<long long> laquinina_cells;
    std::unordered_set<long long> ondasdelcaribe_cells;
    std::unordered_set<long long> pantano_cells;
    std::unordered_set<long long> ochodediciembre_cells;
    std::unordered_set<long long> lacoquera_cells;
    std::unordered_set<long long> sanjacinto_cells;
    std::unordered_set<long long> nuevabethel_cells;
    // Track which discrete cells are active (derived from patchCenters)
    std::unordered_set<long long> active_cells;

    // Temperature time series
    std::vector<std::vector<int>> dataTemperatures;
    // Dedicated RNG for daily temperature draws. Kept SEPARATE from repast::Random
    // (the global stream that drives infection/movement) so sampling temperature each
    // tick does not perturb the epidemic's Common-Random-Number streams. Seeded from
    // run_seed for per-run reproducibility.
    std::mt19937 tempRng_;
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

    // Map neighborhood_id -> HI mean in [0,1]
    std::unordered_map<int, double> hi_prob_by_neigh;

    // Pack coords into a 64-bit key for the sets
    inline long long pack_xy(int x, int y) const {
        return (static_cast<long long>(x) << 32) | (static_cast<unsigned long long>(y));
    }

public:
    RepastHPCModel(std::string propsFile, int argc, char** argv, boost::mpi::communicator* comm, const Params& params);
    ~RepastHPCModel();

    // --- Neighborhood assignment API ---
    void assignNeighborhoodCluster(const std::string& name, int nid,
                                int cx, int cy, int target_cells);
    void resetNeighborhoodCounters();
    // Convert lat/lon to grid coordinates
    std::pair<int,int> latlon_to_grid(double lat, double lon) const;    
    void init();
    void initHumans();
    void initGridValueLayers();
    void loadHIFromCSV(const std::string& csv_path);
    void writeCsv(const std::string& filename);
    void runAllHumans();
    void runAllPatches();
    void initRecords();
    void recordResults();
    void initSchedule(repast::ScheduleRunner& runner);
    double get_border_x(int coordY);
    void printExecutionTime(); 
    void validateNeighborhoodsResidentialOnly();

    // helpers
    int ageInitializer(double prob);
    int* countHumansInState();
    std::vector<int> homeLocationInitializer();
    std::vector<std::vector<int>> activitiesInitializer(int age);
    std::vector<int> placeSelector(int location_value);
};

#endif // MYMODEL_H

