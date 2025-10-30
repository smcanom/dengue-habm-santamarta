#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <random>
#include <queue>
#include <vector>
#include <iomanip>
#include <stdexcept>
#include <boost/mpi.hpp>
#include "repast_hpc/AgentId.h"
#include "repast_hpc/RepastProcess.h"
#include "repast_hpc/Utilities.h"
#include "repast_hpc/Properties.h"
#include "repast_hpc/initialize_random.h"
#include "repast_hpc/Point.h"
#include <utility>
#include "MyModel.h"
#include "SEIModel.h"
#include "MyReadData.h"



//OBLIGATORY (NO SE MUY BIEN PARA QUE ES)
HumanPackageProvider::HumanPackageProvider(repast::SharedContext<Human>* agentPtr): agents(agentPtr){ }

//OBLIGATORY (NO SE MUY BIEN PARA QUE ES)
void HumanPackageProvider::providePackage(Human * agent, std::vector<HumanPackage>& out){
	repast::AgentId id = agent->getId();
	HumanPackage package(id.id(), id.startingRank(), id.agentType(), id.currentRank(),agent->getInfectionState(), agent->getAge(), agent->getTimeSinceSuccesfullBite(), agent->getTimeSinceInfection(), agent->getHomeLocation(), agent->getActivities());
	out.push_back(package);
}

//OBLIGATORY (NO SE MUY BIEN PARA QUE ES)
void HumanPackageProvider::provideContent(repast::AgentRequest req, std::vector<HumanPackage>& out){
    std::vector<repast::AgentId> ids = req.requestedAgents();
    for(size_t i = 0; i < ids.size(); i++){
        providePackage(agents->getAgent(ids[i]), out);
    }
}

void RepastHPCModel::resetNeighborhoodCounters() {
	new_cases_oasis = 0;
	new_cases_laquinina = 0;
	new_cases_ondasdelcaribe = 0;
	new_cases_pantano = 0;
	new_cases_ochodediciembre = 0;
	new_cases_lacoquera = 0;
	new_cases_sanjacinto = 0;
	new_cases_nuevabethel = 0;
}

static inline std::string trim_csv(const std::string& s) {
    auto b = s.find_first_not_of(" \t\r\n");
    auto e = s.find_last_not_of(" \t\r\n");
    return (b == std::string::npos) ? "" : s.substr(b, e - b + 1);
}

void RepastHPCModel::loadHIFromCSV(const std::string& csv_path) {
    hi_prob_by_neigh.clear();

    std::ifstream in(csv_path.c_str());
    if (!in) throw std::runtime_error("Cannot open HI CSV at: " + csv_path);

    std::string line;
    if (!std::getline(in, line)) throw std::runtime_error("Empty HI CSV: " + csv_path);
    // Expect header: neighborhood_id,neighborhood_name,hi_mean

    while (std::getline(in, line)) {
        if (trim_csv(line).empty()) continue;
        std::istringstream ss(line);
        std::string id_str, name_str, hi_str;

        if (!std::getline(ss, id_str, ',')) continue;
        if (!std::getline(ss, name_str, ',')) continue;
        if (!std::getline(ss, hi_str, ','))  continue;

        int nid = std::stoi(trim_csv(id_str));
        double hi = std::stod(trim_csv(hi_str));     // already in [0,1]
        if (hi < 0.0) hi = 0.0;
        if (hi > 1.0) hi = 1.0;

        hi_prob_by_neigh[nid] = hi;
    }
    if (hi_prob_by_neigh.empty())
        throw std::runtime_error("No rows parsed from HI CSV: " + csv_path);
    
    std::cout << "Loaded HI data for " << hi_prob_by_neigh.size() << " neighborhoods" << std::endl;
}

std::pair<int,int> RepastHPCModel::latlon_to_grid(double lat, double lon) const {
    // Latitude: 11.075 → y=0, 11.250 → y=200
    // y = (lat - 11.075) / 0.175 * 200
    int y = static_cast<int>(std::round((lat - 11.075) / 0.175 * 200.0));
    
    // Longitude: -74.240 → x=0, -74.140 → x=135
    // x = (lon - (-74.240)) / 0.100 * 135
    int x = static_cast<int>(std::round((lon - (-74.240)) / 0.100 * 135.0));
    
    return {x, y};
}

void RepastHPCModel::assignNeighborhoodCluster(const std::string& name, int nid,
                                               int cx, int cy, int target_cells) {
    using Pt = std::pair<int,int>;
    std::queue<Pt> q;
    std::unordered_set<long long> visited;

    // try_enqueue: residential-only, active, unassigned cells
    auto try_enqueue = [&](int x, int y) {
        if (x < 0 || x >= xdim || y < 0 || y >= ydim) return;
        long long key = pack_xy(x,y);
        if (visited.count(key)) return;
        if (!active_cells.count(key)) return; // only active cells

        bool err = false;
        repast::Point<int> pt(x,y);
        int t = valueLayerType->getValueAt(pt, err);
        if (err || t != 1) return;            // residential only (Type==1)
        int cur = valueLayerNeighborhood->getValueAt(pt, err);
        if (err || cur != 0) return;          // avoid overlap
        visited.insert(key);
        q.emplace(x,y);
    };

    auto find_residential_seed = [&]() -> bool {
        const int max_radius = std::max(xdim, ydim);
        for (int r = 0; r <= max_radius; ++r) {
            for (int dx = -r; dx <= r; ++dx) {
                for (int dy = -r; dy <= r; ++dy) {
                    if (std::abs(dx) != r && std::abs(dy) != r) continue; // ring border
                    int x = cx + dx;
                    int y = cy + dy;
                    if (x < 0 || x >= xdim || y < 0 || y >= ydim) continue;
                    long long key = pack_xy(x,y);
                    if (!active_cells.count(key)) continue;
                    bool err = false;
                    repast::Point<int> pt(x,y);
                    int t = valueLayerType->getValueAt(pt, err);
                    if (!err && t == 1) {
                        int cur = valueLayerNeighborhood->getValueAt(pt, err);
                        if (!err && cur == 0) {
                            visited.insert(key);
                            q.emplace(x,y);
                            return true;
                        }
                    }
                }
            }
        }
        return false;
    };

    if (!find_residential_seed()) {
        std::cerr << "[Neighborhood " << name << "] Unable to find a residential seed near ("
                  << cx << "," << cy << ")" << std::endl;
        return;
    }

    int filled = 0;
    while (!q.empty() && filled < target_cells) {
        auto cur = q.front(); q.pop();
        int x = cur.first, y = cur.second;
        bool err = false;
        repast::Point<int> pt(x,y);
        int existing = valueLayerNeighborhood->getValueAt(pt, err);
        if (!err && existing == 0) {
            // Assign to neighborhood
            valueLayerNeighborhood->setValueAt(nid, pt, err);
            long long key = pack_xy(x,y);
            
            // Record in appropriate set
            if (nid == 1) oasis_cells.insert(key);
            else if (nid == 2) laquinina_cells.insert(key);
            else if (nid == 3) ondasdelcaribe_cells.insert(key);
            else if (nid == 4) pantano_cells.insert(key);
            else if (nid == 5) ochodediciembre_cells.insert(key);
            else if (nid == 6) lacoquera_cells.insert(key);
            else if (nid == 7) sanjacinto_cells.insert(key);
            else if (nid == 8) nuevabethel_cells.insert(key);
            
            ++filled;

            // 4-neighbor expansion (fixed order for determinism)
            try_enqueue(x+1, y);
            try_enqueue(x-1, y);
            try_enqueue(x, y+1);
            try_enqueue(x, y-1);
        }
    }

    std::cout << "[Neighborhood " << name << "] assigned " << filled << " / " << target_cells 
              << " cells" << std::endl;

    if (filled < target_cells) {
        std::cerr << "WARNING: [Neighborhood " << name << "] requested " << target_cells
                  << " cells, assigned " << filled << " (shortfall: " 
                  << (target_cells - filled) << ")" << std::endl;
    }
}

//OBLIGATORY (NO SE MUY BIEN PARA QUE ES)
HumanPackageReceiver::HumanPackageReceiver(repast::SharedContext<Human>* agentPtr, const Params* params)
    : agents(agentPtr), params_(params) {}

//OBLIGATORY (NO SE MUY BIEN PARA QUE ES)
Human * HumanPackageReceiver::createAgent(HumanPackage package){
    repast::AgentId id(package.id, package.rank, package.type, package.currentRank);
    return new Human(id, package.infectionState, package.age, package.timeSinceSuccesfullBite, package.timeSinceInfection, package.homeLocation, package.activities, params_);

}

//OBLIGATORY (NO SE MUY BIEN PARA QUE ES)
void HumanPackageReceiver::updateAgent(HumanPackage package){
    repast::AgentId id(package.id, package.rank, package.type,package.currentRank);//creo un id
    Human * agent = agents->getAgent(id);//cojo el agente que tiene ese id
    agent->setAll(package.infectionState, package.age, package.timeSinceSuccesfullBite, package.timeSinceInfection, package.homeLocation, package.activities);//le cambio las variables de estado
}

//CONSTRUCTOR 
RepastHPCModel::RepastHPCModel(std::string propsFile, int argc, char** argv, boost::mpi::communicator* comm, const Params& params)
    : context(comm), params_(params) {
	props = new repast::Properties(propsFile, argc, argv, comm); 
	//atributes that come from the model.props file
    start_time = std::chrono::steady_clock::now();
    
	// Compute run seed for CRN (Common Random Numbers)
	run_seed = params_.base_seed + params_.replicate_id;
	
	stopAt = repast::strToInt(props->getProperty("stop.at")); 
	countOfHumans= repast::strToInt(props->getProperty("count.of.humans"));
	countOfInfectedHumans= repast::strToInt(props->getProperty("count.of.infected.humans"));
	//constant_temperature = repast::strToInt(props->getProperty("constant.temperature"));

	xdim=repast::strToInt(props->getProperty("x.dim"));
	ydim=repast::strToInt(props->getProperty("y.dim"));

	// Validate grid dimensions
	if (xdim <= 0 || ydim <= 0) {
		std::cerr << "Error: Invalid grid dimensions. x.dim=" << xdim << ", y.dim=" << ydim << std::endl;
		throw std::runtime_error("Grid dimensions must be positive");
	}

	//manage randomness
	initializeRandom(*props,comm);

	/* CREATE INSTANCE OF SPATIAL PROJECTION */
	// Punto mínimo: [-59348, 4932]
	// Punto máximo: [-59230, 5044]
	//define the dimensions


	//repast::Point<double> origin(59230, 4932);
	//repast::Point<double> extent (59348, 5044);//this will extent xdim units in the x asis and ydim units in the y axis
	//repast::GridDimensions gd(origin,extent);
	
	repast::Point<double> origin(0, 0);
	repast::Point<double> extent(xdim, ydim);  // Use dimensions from properties
	repast::GridDimensions gd(origin,extent);

	//define process dimensions based on communicator size
	std::vector<int> processDims;
	int world_size = comm->size();
	int dims_x = std::max(1, (int)std::floor(std::sqrt((double)world_size)));
	while (world_size % dims_x != 0) { --dims_x; }
	int dims_y = world_size / dims_x;
	processDims.push_back(dims_x);
	processDims.push_back(dims_y);
	// e.g., 1->(1,1), 4->(2,2), 8->(2,4) etc.

	ReadData reader;
 
    patchCenters = reader.getPatchCenters("santa_marta_grid_relative.csv");
    std::cout
    << "Number of patches (centroids): " << patchCenters.size() << "\n";

    border_points = {
            {190,25},
            {188,35},
            {175,28},
            {155, 0},   // Top-middle dip
            {130, 15},  // First rightward dip
            {110, 0},   // Leftward recovery
            {75, 20},   // Second rightward dip
            {45, 10},   // Leftward movement
            {0, 30}     // Bottom position
        };

	//create a discrete space
	discreteSpace=new repast::SharedDiscreteSpace<Human,repast::StrictBorders , repast::SimpleAdder<Human> >("AgentDiscreteSpace",gd,processDims,0,comm);
	//we are chosing strict borders instead of wrap around
	//0 makes reference to the length of the buffer zone. In this case it is zero
	// comm= is the instance of the  boost::mpi::communicator*
	
	//add discreteSpace to context
	context.addProjection(discreteSpace);


	/* CREATE GRID VALUE LAYERS */
	//ValueLayerND(vector<int> processesPerDim, GridDimensions globalBoundaries, int bufferSize,bool periodic);
	int bufferSize=xdim/2;
	valueLayerSuceptibleMosquitoes=new repast::ValueLayerND<int>(processDims, gd,bufferSize,false);
	valueLayerExposedMosquitoes=new repast::ValueLayerND<int>(processDims, gd, bufferSize, false);
	valueLayerInfectedMosquitoes=new repast::ValueLayerND<int>(processDims, gd, bufferSize, false);    
	valueLayerTemperature=new repast::ValueLayerND<int>(processDims, gd, bufferSize, false); 
	//valueLayerTemperature = new repast::ValueLayerND<int>(processDims, gd, bufferSize, false);

    
	valueLayerType=new repast::ValueLayerND<int>(processDims, gd, bufferSize, false);
	// Neighborhood tagging layer: 0 = none, 1 = Oasis, 2 = La Quinina
	valueLayerNeighborhood = new repast::ValueLayerND<int>(processDims, gd, bufferSize, false);
	// HI activation mask layer
	valueLayerHI_activated = new repast::ValueLayerND<int>(processDims, gd, bufferSize, false);

	//These value layers store the total mosquitoes on t-2.
	valuetotalmosquitoes0 = new repast::ValueLayerND<int>(processDims, gd, bufferSize, false);
	valuetotalmosquitoes1  = new repast::ValueLayerND<int>(processDims, gd, bufferSize, false);

    // Initialize mosquitoPatches as a flat array:
    mosquitoPatches.resize(patchCenters.size());
    for (size_t i = 0; i < patchCenters.size(); ++i) {
        mosquitoPatches[i] = SEIModel(0,0,0,0, /*delay=*/11, /*rho=*/0.5, /*omega=*/1.0, &params_);
    }
	

	// Build locToIndex mapping from centroid → patch index:
    locToIndex.clear();
    for (size_t i = 0; i < patchCenters.size(); ++i) {
        int x = static_cast<int>(patchCenters[i].first);
        int y = static_cast<int>(patchCenters[i].second);
        long long key = (static_cast<long long>(x) << 32) | (unsigned long long)y;
        locToIndex[key] = { static_cast<int>(i), 0 };
    }

	//read the temperature file and create the list that contains the data ----->NEW
	ReadData myData;
	//std::string fileName ="eratura-precipitacion-bello.csv";
	std::string fileName ="SantaMartaChange.csv";
	int simulationTime(stopAt);
	dataTemperatures=myData.loadFromExcel(fileName,simulationTime);
	
	//NO SE MUY BIEN PARA QUE SON LAS DOS SIGUIENTES LINEAS
	provider=new HumanPackageProvider(&context);
	receiver=new HumanPackageReceiver(&context, &params_);
	
}

RepastHPCModel::~RepastHPCModel() {
    // Clear agents from context first
    // The context will manage agent cleanup
    
    // Delete agent communication helpers BEFORE deleting value layers
    if (provider) {
        delete provider;
        provider = nullptr;
    }
    if (receiver) {
        delete receiver;
        receiver = nullptr;
    }

    // Delete value layers
    if (valueLayerSuceptibleMosquitoes) {
        delete valueLayerSuceptibleMosquitoes;
        valueLayerSuceptibleMosquitoes = nullptr;
    }
    if (valueLayerExposedMosquitoes) {
        delete valueLayerExposedMosquitoes;
        valueLayerExposedMosquitoes = nullptr;
    }
    if (valueLayerInfectedMosquitoes) {
        delete valueLayerInfectedMosquitoes;
        valueLayerInfectedMosquitoes = nullptr;
    }
    if (valueLayerTemperature) {
        delete valueLayerTemperature;
        valueLayerTemperature = nullptr;
    }
    if (valueLayerType) {
        delete valueLayerType;
        valueLayerType = nullptr;
    }
    if (valueLayerNeighborhood) {
        delete valueLayerNeighborhood;
        valueLayerNeighborhood = nullptr;
    }
    if (valuetotalmosquitoes0) {
        delete valuetotalmosquitoes0;
        valuetotalmosquitoes0 = nullptr;
    }
    if (valuetotalmosquitoes1) {
        delete valuetotalmosquitoes1;
        valuetotalmosquitoes1 = nullptr;
    }
    if (valueLayerHI_activated) {
        delete valueLayerHI_activated;
        valueLayerHI_activated = nullptr;
    }

    // Delete properties
    if (props) {
        delete props;
        props = nullptr;
    }

    // Clear containers
    mosquitoPatches.clear();
    locToIndex.clear();
    patchCenters.clear();
    dataTemperatures.clear();
    residentialLocations.clear();
    studyLocations.clear();
    workLocations.clear();
    otherActivitiesLocations.clear();
    oasis_cells.clear();
    laquinina_cells.clear();
    ondasdelcaribe_cells.clear();
    pantano_cells.clear();
    ochodediciembre_cells.clear();
    lacoquera_cells.clear();
    sanjacinto_cells.clear();
    nuevabethel_cells.clear();
    active_cells.clear();
}
 void RepastHPCModel::init(){
	// Load HI data before initializing grid (hardcoded path like other input files)
	loadHIFromCSV("hi_by_neighborhood_2023.csv");
	
	initGridValueLayers();
	initHumans();
	writeCsv("locations.csv");
	initRecords();
} 
//METODO PARA INIZIALIZAR LOS HUMANOS
void RepastHPCModel::initHumans(){
	std::cout <<"incializacion Humanos\n";
	int rank = repast::RepastProcess::instance()->rank();
	std::cout<<"rank: "<<rank<<"\n"; 
	
	if(rank == 0){
		for (int i=0; i<countOfHumans; i++){
			
			//SE DEFINE EL ID 
			repast::AgentId id(i,rank,0);
			id.currentRank(rank);

			//Se define la edad con ayuda del metodo ageInitializer
			double myRand=repast::Random::instance()->nextDouble();
			int age=ageInitializer(myRand);
			//std::cout<<myRand<<" "<<rank<<"\n";
			
			//se define la homelocation con ayuda del metodo homeLocationInitializer
			std::vector<int>homeLocation=homeLocationInitializer();
			
			//se definen las actividades con ayuda del metodo activitiesInitializer
			std::vector<std::vector<int>> activities=activitiesInitializer(age); 
			
			//se define todo lo relacionado al SEIR status
			string infectionState;
			int timeSinceSuccesfullBite;
			int timeSinceInfection;
			if(i<countOfHumans-countOfInfectedHumans){
				infectionState="susceptible";
				timeSinceSuccesfullBite=-1; //equivalent to Null
				timeSinceInfection=-1; //equivalent to null
			}
			else{
				infectionState="infected";
				timeSinceSuccesfullBite=0;
				timeSinceInfection=0;
			}
			//SE CREA EL HUMANO Y SE METE AL CONTEXT 
			Human* H=new Human(id,infectionState,age,timeSinceSuccesfullBite,timeSinceInfection,homeLocation,activities,&params_);
			context.addAgent(H);

			//SE UBICA EL HUMANO EN LA HOME LOCATION
			discreteSpace->moveTo(id,homeLocation);
			//repast::RepastProcess::instance()->synchronizeAgentStatus<Human, HumanPackage,HumanPackageProvider,HumanPackageReceiver>(context, *provider,*receiver,*receiver);



			//mirar si el humano si se ubico donde era
			std::vector<int>currentLocation;
			discreteSpace->getLocation(H->getId(),currentLocation);
			//imprimir humano
			//std::cout<<"id: "<<H->getId()<<" age: "<<H->getAge()<<", infecion state: "<<H->getInfectionState()<<", timeSinceSuccesfullBite: "<<H->getTimeSinceSuccesfullBite()<<", timeSinceInfection: "<<H->getTimeSinceInfection()<<", homeLocation: ("<<H->getHomeLocation()[0]<<","<<H->getHomeLocation()[1]<<"), actividad1:("<<H->getActivities()[0][0]<<","<<H->getActivities()[0][1]<<"), actividad2:("<<H->getActivities()[1][0]<<","<<H->getActivities()[1][1]<<"),currentLocation:("<<currentLocation[0]<<";"<<currentLocation[1]<<")\n"; 
		}
	}
	
}

void RepastHPCModel::writeCsv(const std::string& filename) {
    std::ofstream file(filename);
    if (!file) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Escribe los encabezados
    file << "Category,X,Y\n";

    // Escribe datos de cada categoría
    auto writeCategoryData = [&file](const std::vector<std::vector<int>>& data, const std::string& category) {
        for (const auto& row : data) {
            if (row.size() == 2) { // Asegúrate de que cada fila tenga exactamente dos valores
                file << category << "," << row[0] << "," << row[1] << "\n";
            } else {
                std::cerr << "Error: Data row does not contain exactly two values.\n";
            }
        }
    };

    writeCategoryData(residentialLocations, "Residential");
    writeCategoryData(studyLocations, "Study");
    writeCategoryData(workLocations, "Work");
    writeCategoryData(otherActivitiesLocations, "Other Activities");

    file.close();
}

int triangularSampleInt(int a, int c, int b) {
    // a = minimum, c = mode, b = maximum
    // Must satisfy a <= c <= b
    // Returns a single integer sample from the triangular distribution

    double u = repast::Random::instance()->nextDouble();  // uniform(0,1)
    double fc = double(c - a) / double(b - a);            // fraction for the mode

    double x;
    if (u < fc) {
        // left side of the triangle
        x = a + std::sqrt(u * (b - a) * (c - a));
    } else {
        // right side of the triangle
        x = b - std::sqrt((1.0 - u) * (b - a) * (b - c));
    }

    // Convert to int using rounding. Alternatives:
    //   - static_cast<int>(std::floor(x))  // always round down
    //   - static_cast<int>(std::ceil(x))   // always round up
    //   - static_cast<int>(std::floor(x + 0.5)) // standard rounding
    return static_cast<int>(std::round(x));
}

// 3. Rewrite initGridValueLayers() to loop over patchCenters:

void RepastHPCModel::initGridValueLayers() {
    std::cout << "Initializing value layers...\n";
    size_t N = patchCenters.size();
    double minRawX = patchCenters.front().first;
    double maxRawX = patchCenters.front().first;
    double minRawY = patchCenters.front().second;
    double maxRawY = patchCenters.front().second;
    for (auto& p : patchCenters) {
    minRawX = std::min(minRawX, p.first);
    maxRawX = std::max(maxRawX, p.first);
    minRawY = std::min(minRawY, p.second);
    maxRawY = std::max(maxRawY, p.second);
    }
    size_t currentLocation = 0;
    double count_susceptibles = 0.0, count_infected = 0.0;

    for (size_t i = 0; i < N; ++i) {
        // 3.1 Fetch centroid
        double rawX = patchCenters[i].first;
        double rawY = patchCenters[i].second;

        // normalize into [0,1]
        double normX = (rawX - minRawX) / (maxRawX - minRawX);
        double normY = (rawY - minRawY) / (maxRawY - minRawY);

        // 3.2 Scale to grid
        int coordX = static_cast<int>(std::round(normX * (xdim - 1)));
        int coordY = static_cast<int>(std::round(normY * (ydim - 1)));
        std::vector<int> location = { coordX, coordY };



        // 3.3 Initialize temperature, mosquitoes, type, etc.
        bool err = false;
        // inside your initGridValueLayers() loop, replace the hard-coded temperature with this:
        int maxTempInitialDay = dataTemperatures.empty() || dataTemperatures[0].empty()
            ? 30
            : dataTemperatures[0][0];
        int minTempInitialDay = dataTemperatures.size() < 2 || dataTemperatures[1].empty()
            ? 20
            : dataTemperatures[1][0];
        double temperature = repast::Random::instance()
            ->createUniDoubleGenerator(minTempInitialDay, maxTempInitialDay)
            .next();
        valueLayerTemperature->setValueAt(temperature, location, err);

        int rand_susc = triangularSampleInt(0,15,30);
        int rand_exp  = 0;
        int rand_inf  = repast::Random::instance()
                          ->createUniIntGenerator(0,2).next();
        int total = rand_inf + rand_susc;
        valuetotalmosquitoes0->setValueAt(total, location, err);
        valuetotalmosquitoes1->setValueAt(total, location, err);
        // then synchronize once, so both layers are non‐zero at tick=1
        // (sync once after initialization, not per-cell)
        valueLayerSuceptibleMosquitoes->setValueAt(rand_susc, location, err);
        valueLayerExposedMosquitoes->setValueAt(rand_exp, location, err);
        valueLayerInfectedMosquitoes->setValueAt(rand_inf, location, err);
        count_susceptibles += rand_susc;
        count_infected += rand_inf;
        int delay = 18;
        double rho = 0.46;
        // 3.4 Seed SEIModel

        mosquitoPatches[i] = SEIModel(rand_susc, rand_exp, rand_inf,
                                      static_cast<int>(temperature),
                                      delay, rho, /*omega=*/1.0, &params_);

        // 3.5 Assign patch type
        // 3.5 Assign patch type with beach prioritization
        int type;
        
        // Parameters
        const double buffer_zone = 100;  // How far inland work locations extend
        const double max_prob = 0.85;   // Probability at the border
        const double min_prob = 0.15;   // Probability at buffer_zone distance

        // For current patch (coordX, coordY):
        double border_x = this->get_border_x(coordY);
        double distance = coordX - border_x;  // Distance from border

        if (distance >= 0 && distance <= buffer_zone) {
            // Exponential probability decay
            double work_prob = max_prob * exp(-distance / 8.0);
            
            if (repast::Random::instance()->nextDouble() < work_prob) {
                type = 3; // Work location
            } else {
                // Assign other types proportionally
                double u2 = repast::Random::instance()->nextDouble();
                type = (u2 <= 0.927 ? 1 : u2 <= 0.967 ? 2 : 4);
            }
        } else {
            // Original distribution for non-border patches
            double u = repast::Random::instance()->nextDouble();
            type = (u <= 0.896090 ? 1 :
                    u <= 0.934925 ? 2 :
                    u <= 0.968250 ? 3 : 4);
        }

        // Store type and location
        switch(type) {
            case 1: residentialLocations.push_back(location); break;
            case 2: studyLocations.push_back(location); break;
            case 3: workLocations.push_back(location); break;
            case 4: otherActivitiesLocations.push_back(location); break;
        }

        bool errType = false;
        valueLayerType->setValueAt(type, location, errType);
        // Track active cells for neighborhood assignment
        active_cells.insert(pack_xy(coordX, coordY));
        ++currentLocation;
        // Calcular y mostrar el progreso
        double progress = (static_cast<double>(currentLocation) / N) * 100.0;
        std::cout << "Progress: " << progress << "%\r";
        std::cout.flush();

    }

    // 3.6 Synchronize all layers
    valueLayerSuceptibleMosquitoes->synchronize();
    valueLayerExposedMosquitoes->synchronize();
    valueLayerInfectedMosquitoes->synchronize();
    valueLayerTemperature->synchronize();
    valueLayerType->synchronize();
    valuetotalmosquitoes0->synchronize();
    valuetotalmosquitoes1->synchronize();

    // 3.7 Assign neighborhoods after type layer is populated
    {
        std::cout << "\n=== Assigning Neighborhoods ===" << std::endl;
        
        // Existing neighborhoods (manual seeds)
        auto oasis_coords = latlon_to_grid(11.238449, -74.164696);
        assignNeighborhoodCluster("Oasis", 1, oasis_coords.first, oasis_coords.second, 108);
        auto laquinina_coords = latlon_to_grid(11.18233,-74.21992);
        assignNeighborhoodCluster("La Quinina", 2, laquinina_coords.first, laquinina_coords.second, 172);
        
        // New neighborhoods (from lat/lon centroids)
        // Ondas del Caribe: 11.241997, -74.167821
        auto ondas_coords = latlon_to_grid(11.241997, -74.167821);
        assignNeighborhoodCluster("Ondas del Caribe", 3, ondas_coords.first, ondas_coords.second, 100);
        
        // Pantano: 11.239484, -74.173353
        auto pantano_coords = latlon_to_grid(11.239484, -74.173353);
        assignNeighborhoodCluster("Pantano", 4, pantano_coords.first, pantano_coords.second, 103);
        
        // 8 de Diciembre: 11.244997, -74.164336
        auto ocho_coords = latlon_to_grid(11.244997, -74.164336);
        assignNeighborhoodCluster("8 de Diciembre", 5, ocho_coords.first, ocho_coords.second, 147);
        
        // La Coquera: 11.187200, -74.222770
        auto coquera_coords = latlon_to_grid(11.187200, -74.222770);
        assignNeighborhoodCluster("La Coquera", 6, coquera_coords.first, coquera_coords.second, 147);
        
        // San Jacinto: 11.190480, -74.221800
        auto sanjacinto_coords = latlon_to_grid(11.190480, -74.221800);
        assignNeighborhoodCluster("San Jacinto", 7, sanjacinto_coords.first, sanjacinto_coords.second, 147);
        
        // Nueva Bethel: 11.190158, -74.214619
        auto bethel_coords = latlon_to_grid(11.190158, -74.214619);
        assignNeighborhoodCluster("Nueva Bethel", 8, bethel_coords.first, bethel_coords.second, 122);

        valueLayerNeighborhood->synchronize();

        // Validate: all neighborhood cells must be residential (type=1)
        std::cout << "\n=== Validating Neighborhoods ===" << std::endl;
        validateNeighborhoodsResidentialOnly();

        // Export CSV with all neighborhoods
        try {
            std::ofstream dbg("neighborhood_cells.csv", std::ofstream::trunc);
            dbg << "neighborhood_id,neighborhood_name,x,y\n";
            
            for (const auto& key : oasis_cells) {
                int x = static_cast<int>(key >> 32);
                int y = static_cast<int>(key & 0xFFFFFFFFULL);
                dbg << "1,Oasis," << x << "," << y << "\n";
            }
            for (const auto& key : laquinina_cells) {
                int x = static_cast<int>(key >> 32);
                int y = static_cast<int>(key & 0xFFFFFFFFULL);
                dbg << "2,La Quinina," << x << "," << y << "\n";
            }
            for (const auto& key : ondasdelcaribe_cells) {
                int x = static_cast<int>(key >> 32);
                int y = static_cast<int>(key & 0xFFFFFFFFULL);
                dbg << "3,Ondas del Caribe," << x << "," << y << "\n";
            }
            for (const auto& key : pantano_cells) {
                int x = static_cast<int>(key >> 32);
                int y = static_cast<int>(key & 0xFFFFFFFFULL);
                dbg << "4,Pantano," << x << "," << y << "\n";
            }
            for (const auto& key : ochodediciembre_cells) {
                int x = static_cast<int>(key >> 32);
                int y = static_cast<int>(key & 0xFFFFFFFFULL);
                dbg << "5,8 de Diciembre," << x << "," << y << "\n";
            }
            for (const auto& key : lacoquera_cells) {
                int x = static_cast<int>(key >> 32);
                int y = static_cast<int>(key & 0xFFFFFFFFULL);
                dbg << "6,La Coquera," << x << "," << y << "\n";
            }
            for (const auto& key : sanjacinto_cells) {
                int x = static_cast<int>(key >> 32);
                int y = static_cast<int>(key & 0xFFFFFFFFULL);
                dbg << "7,San Jacinto," << x << "," << y << "\n";
            }
            for (const auto& key : nuevabethel_cells) {
                int x = static_cast<int>(key >> 32);
                int y = static_cast<int>(key & 0xFFFFFFFFULL);
                dbg << "8,Nueva Bethel," << x << "," << y << "\n";
            }
            
            std::cout << "✓ Exported neighborhood_cells.csv" << std::endl;
        } catch (...) {
            std::cerr << "Warning: failed to write neighborhood_cells.csv" << std::endl;
        }
        
        // Export summary
        try {
            std::ofstream summary("summary_neighborhood_sizes.csv", std::ofstream::trunc);
            summary << "neighborhood_id,name,target_cells,realized_cells\n";
            summary << "1,Oasis,108," << oasis_cells.size() << "\n";
            summary << "2,La Quinina,172," << laquinina_cells.size() << "\n";
            summary << "3,Ondas del Caribe,100," << ondasdelcaribe_cells.size() << "\n";
            summary << "4,Pantano,103," << pantano_cells.size() << "\n";
            summary << "5,8 de Diciembre,147," << ochodediciembre_cells.size() << "\n";
            summary << "6,La Coquera,147," << lacoquera_cells.size() << "\n";
            summary << "7,San Jacinto,147," << sanjacinto_cells.size() << "\n";
            summary << "8,Nueva Bethel,122," << nuevabethel_cells.size() << "\n";
            
            std::cout << "✓ Exported summary_neighborhood_sizes.csv" << std::endl;
        } catch (...) {
            std::cerr << "Warning: failed to write summary_neighborhood_sizes.csv" << std::endl;
        }
        
        // 3.8 Initialize HI-based stochastic activation mask
        std::cout << "\n=== Initializing HI-Based Activation Mask ===" << std::endl;
        
        // Read configuration parameters
        bool hi_activate = (props->getProperty("hi.activate") == "true");
        double eta = repast::strToDouble(props->getProperty("hi.eta"));
        int hi_seed = repast::strToInt(props->getProperty("hi.seed"));
        
        std::cout << "hi.activate = " << (hi_activate ? "true" : "false") << std::endl;
        std::cout << "hi.eta (baseline) = " << eta << std::endl;
        std::cout << "hi.seed = " << hi_seed << std::endl;
        
        // Create a dedicated RNG for HI activation (seeded independently)
        std::mt19937 hi_rng(hi_seed);
        std::uniform_real_distribution<double> uniform_01(0.0, 1.0);
        
        int activated_count = 0;
        int total_residential = 0;
        
        // Prepare CSV output (rank 0 only)
        std::vector<std::string> activation_csv_lines;
        if (repast::RepastProcess::instance()->rank() == 0) {
            activation_csv_lines.push_back("x,y,neighborhood_id,HI_bar,X_i,omega_i");
        }
        
        for (size_t i = 0; i < patchCenters.size(); ++i) {
            double rawX = patchCenters[i].first;
            double rawY = patchCenters[i].second;
            
            // Recompute grid coordinates
            double normX = (rawX - minRawX) / (maxRawX - minRawX);
            double normY = (rawY - minRawY) / (maxRawY - minRawY);
            int coordX = static_cast<int>(std::round(normX * (xdim - 1)));
            int coordY = static_cast<int>(std::round(normY * (ydim - 1)));
            repast::Point<int> pt(coordX, coordY);
            
            bool err = false;
            int type = valueLayerType->getValueAt(pt, err);
            
            // Initialize default values
            int X_i = 0;
            double omega_i = 1.0;
            
            // Only process residential patches (type == 1) if HI is activated
            if (!err && type == 1 && hi_activate) {
                total_residential++;
                
                // Get neighborhood tag
                int nid = valueLayerNeighborhood->getValueAt(pt, err);
                
                // Lookup HI probability
                double hi_bar = 0.0;
                if (auto it = hi_prob_by_neigh.find(nid); it != hi_prob_by_neigh.end()) {
                    hi_bar = it->second;
                }
                
                // Draw Bernoulli(hi_bar) using dedicated HI RNG
                X_i = (uniform_01(hi_rng) < hi_bar) ? 1 : 0;
                if (X_i == 1) activated_count++;
                
                // Compute omega = eta + (1 - eta) * X_i
                omega_i = eta + (1.0 - eta) * X_i;
                
                // Store X_i in value layer for demonstration
                valueLayerHI_activated->setValueAt(X_i, pt, err);
                
                // Record for CSV (rank 0 only)
                if (repast::RepastProcess::instance()->rank() == 0) {
                    std::ostringstream oss;
                    oss << coordX << "," << coordY << "," << nid << "," 
                        << std::fixed << std::setprecision(6) << hi_bar << "," 
                        << X_i << "," << omega_i;
                    activation_csv_lines.push_back(oss.str());
                }
            } else {
                // Non-residential or HI disabled: X_i=0, omega=1.0
                valueLayerHI_activated->setValueAt(0, pt, err);
            }
            
            // Set emergence multiplier on SEIModel
            mosquitoPatches[i].setEmergenceMultiplier(omega_i);
        }
        
        // Synchronize HI layer
        valueLayerHI_activated->synchronize();
        
        // Write activation map CSV (rank 0 only) - this is the key output for demonstration
        if (repast::RepastProcess::instance()->rank() == 0) {
            try {
                std::ofstream out("hi_activation_map.csv", std::ofstream::trunc);
                for (const auto& line : activation_csv_lines) {
                    out << line << "\n";
                }
                out.close();
                std::cout << "✓ Wrote hi_activation_map.csv with " 
                         << (activation_csv_lines.size() - 1) << " patches" << std::endl;
            } catch (...) {
                std::cerr << "Warning: failed to write hi_activation_map.csv" << std::endl;
            }
        }
        
        double activation_rate = (total_residential > 0) 
            ? (static_cast<double>(activated_count) / total_residential) 
            : 0.0;
        std::cout << "Residential patches: " << total_residential << std::endl;
        std::cout << "Activated patches (X=1): " << activated_count 
                  << " (" << (activation_rate * 100.0) << "%)" << std::endl;
        std::cout << "✓ HI activation mask initialized" << std::endl;
    }
}

void RepastHPCModel::validateNeighborhoodsResidentialOnly(){
    auto check_set = [&](const std::unordered_set<long long>& cells, const char* name){
        int non_residential = 0;
        for (const auto& key : cells) {
            int x = static_cast<int>(key >> 32);
            int y = static_cast<int>(key & 0xFFFFFFFFULL);
            bool err = false;
            repast::Point<int> pt(x,y);
            int t = valueLayerType->getValueAt(pt, err);
            if (err || t != 1) {
                ++non_residential;
                std::cerr << "ERROR: " << name << " has non-residential cell at (" 
                          << x << "," << y << "), Type=" << t << std::endl;
            }
        }
        if (non_residential > 0) {
            std::cerr << "VALIDATION FAILED: " << name << " has " << non_residential
                      << " non-residential cells!" << std::endl;
        }
    };
    check_set(oasis_cells, "Oasis");
    check_set(laquinina_cells, "La Quinina");
    check_set(ondasdelcaribe_cells, "Ondas del Caribe");
    check_set(pantano_cells, "Pantano");
    check_set(ochodediciembre_cells, "8 de Diciembre");
    check_set(lacoquera_cells, "La Coquera");
    check_set(sanjacinto_cells, "San Jacinto");
    check_set(nuevabethel_cells, "Nueva Bethel");
}
double RepastHPCModel::get_border_x(int coordY) {
    // Handle edge cases
    if (coordY >= 190) return 25;
    if (coordY <= 0) return 35;

    // Safety check for empty vector
    if (border_points.empty()) return 35;

    // Find segment containing coordY
    for (size_t i = 0; i < border_points.size() - 1; i++) {
        int y1 = border_points[i].first;
        int x1 = border_points[i].second;
        int y2 = border_points[i+1].first;
        int x2 = border_points[i+1].second;

        if (coordY <= y1 && coordY >= y2) {
            // Linear interpolation
            double fraction = static_cast<double>(y1 - coordY) / (y1 - y2);
            return x1 + fraction * (x2 - x1);
        }
    }
    return 35; // Fallback
}

void RepastHPCModel::runAllHumans() {
  try {
    std::vector<Human*> humans;
    context.selectAgents(countOfHumans, humans);

    static const std::array<std::string,3> phaseLabels = {
      "Primary activity", "Leisure activity", "Home return"
    };

	// Print current tick
	int tick_now = repast::RepastProcess::instance()->getScheduleRunner().currentTick();
	std::cout << "=== TICK " << tick_now << " ===\n";

    // Ensure per-tick unique counting for newly infected humans across all phases
    std::unordered_set<long long> counted_new_cases_this_tick;

    for (int phase = 0; phase < 3; ++phase) {
      std::cout << "=== Phase: " << phaseLabels[phase] << " ===\n";

      for (Human* h : humans) {
        repast::AgentId id = h->getId();

        // a) old location
        std::vector<int> oldLoc;
        discreteSpace->getLocation(id, oldLoc);

        // b) compute destination
        std::vector<int> dest;
        if (phase < 2) {
          // activitiesInitializer filled { primary, other }
          auto acts = h->getActivities();
          dest = acts[phase];
        } else {
          dest = h->getHomeLocation();
        }

        // c) move
        discreteSpace->moveTo(id, dest);

        // d) infection & mosquito update
        repast::Point<int> pt(dest[0], dest[1]);
        bool err = false;
        int susc    = valueLayerSuceptibleMosquitoes->getValueAt(pt, err);
        int exposed = valueLayerExposedMosquitoes->getValueAt(pt, err);
        int infd    = valueLayerInfectedMosquitoes->getValueAt(pt, err);
        int temp    = valueLayerTemperature->getValueAt(pt, err);
        if (err) continue;

        std::vector<Human*> here;
        discreteSpace->getObjectsAt(pt, here);
        int nHumans   = here.size();
        int nInfected = std::count_if(
          here.begin(), here.end(),
          [](Human* hh){ return hh->getInfectionState() == "infected"; });

        h->calculateInfectionProbabilityHuman(infd, susc, nHumans);
        // Apply infection update after each movement phase
        h->actualizeSEIRStatus(&context);

        // Count newly infected by neighborhood
        if (h->isNewlyInfected()) {
          long long human_key = (static_cast<long long>(id.id()) << 32)
                                | (unsigned long long)id.startingRank();
          if (!counted_new_cases_this_tick.count(human_key)) {
            bool err_n = false;
            int nid = valueLayerNeighborhood->getValueAt(pt, err_n);
            if (!err_n) {
              if (nid == 1) ++new_cases_oasis;
              else if (nid == 2) ++new_cases_laquinina;
              else if (nid == 3) ++new_cases_ondasdelcaribe;
              else if (nid == 4) ++new_cases_pantano;
              else if (nid == 5) ++new_cases_ochodediciembre;
              else if (nid == 6) ++new_cases_lacoquera;
              else if (nid == 7) ++new_cases_sanjacinto;
              else if (nid == 8) ++new_cases_nuevabethel;
            }
            counted_new_cases_this_tick.insert(human_key);
          }
        }

        int m1 = valuetotalmosquitoes1->getValueAt(pt, err);
        int m0 = valuetotalmosquitoes0->getValueAt(pt, err);
        SEIModel patch(susc, exposed, infd, temp, 18, 0.46, 1.0, &params_);
        patch.setHumans(nHumans);
        patch.setInfectedHumans(nInfected);
        patch.recalculateSEI(0.1, m0, m1);

        valueLayerSuceptibleMosquitoes->setValueAt(
          patch.getSuceptibleMosquitoes(), pt, err);
        valueLayerExposedMosquitoes->setValueAt(
          patch.getExposedMosquitoes(),    pt, err);
        valueLayerInfectedMosquitoes->setValueAt(
          patch.getInfectedMosquitoes(),   pt, err);

        // Update time counters once per tick (after all movement phases)
        // This should happen for all humans, not just in final phase
      }
    }

    // Update time counters for all humans once per tick
    for (Human* h : humans) {
        h->actualizeTimes();
    }

    std::cout << "=== HUMAN MOVEMENT COMPLETE ===\n";

    // Optional: one final MPI synchronization here
    // repast::RepastProcess::instance()
    //   ->synchronizeAgentStatus<Human,HumanPackage,HumanPackageProvider,HumanPackageReceiver>(
    //       context,*provider,*receiver,*receiver);

  } catch (const std::exception& e) {
    std::cerr << "Error in runAllHumans(): " << e.what() << std::endl;
  }
}

// 4. Rewrite runAllPatches() to iterate over centroids only:

void RepastHPCModel::runAllPatches() {
    // 4.1 Synchronize before reading
    valueLayerSuceptibleMosquitoes->synchronize();
    valueLayerInfectedMosquitoes->synchronize();
    valueLayerTemperature->synchronize();
    valueLayerType->synchronize();
    valuetotalmosquitoes0->synchronize();
    valuetotalmosquitoes1->synchronize();

    size_t N = patchCenters.size();
    double minRawX = patchCenters.front().first;
    double maxRawX = patchCenters.front().first;
    double minRawY = patchCenters.front().second;
    double maxRawY = patchCenters.front().second;
    for (auto& p : patchCenters) {
    minRawX = std::min(minRawX, p.first);
    maxRawX = std::max(maxRawX, p.first);
    minRawY = std::min(minRawY, p.second);
    maxRawY = std::max(maxRawY, p.second);
    }
    for (size_t i = 0; i < N; ++i) {
        // 4.2 Fetch centroid & scale
        double rawX = patchCenters[i].first;
        double rawY = patchCenters[i].second;
        // normalize into [0,1]
        double normX = (rawX - minRawX) / (maxRawX - minRawX);
        double normY = (rawY - minRawY) / (maxRawY - minRawY);

        // 3.2 Scale to grid
        int coordX = static_cast<int>(std::round(normX * (xdim - 1)));
        int coordY = static_cast<int>(std::round(normY * (ydim - 1)));
        std::vector<int> location = { coordX, coordY };
        repast::Point<int> pt(coordX, coordY);

        // 4.3 Read values, update SEIModel, and write back
        bool err = false;
        int S = valueLayerSuceptibleMosquitoes->getValueAt(pt, err);
        int E = valueLayerExposedMosquitoes->getValueAt(pt, err);
        int I = valueLayerInfectedMosquitoes->getValueAt(pt, err);
        int T = valueLayerTemperature->getValueAt(pt, err);
        // 1. Check for retrieval errors
        if (err) {
        std::cerr << "Error getting values at patch [" << i 
                    << "] with coords (" << coordX << "," << coordY << ")" 
                    << std::endl;
        continue;
        }

        // 2. Validate counts
        if (S < 0 || E < 0 || I < 0) {
        std::cerr << "Warning: Negative mosquito counts at patch [" << i
                    << "]: S=" << S << ", E=" << E << ", I=" << I 
                    << std::endl;
        S     = std::max(0, S);
        E     = std::max(0, E);
        I = std::max(0, I);
        }

        // 3. Validate temperature
        //    Suppose your model expects T between, say, 0°C and 50°C:
        if (T < 0 || T > 50) {
        std::cerr << "Warning: Implausible temperature at patch [" << i 
                    << "]: T=" << T << "°C. Clamping to [0,50]." 
                    << std::endl;
        T = std::clamp(T, 0, 50);
        }
        {
            std::lock_guard<std::mutex> lock(patchMutex);
            mosquitoPatches[i].setSuceptibleMosquitoes(S);
            mosquitoPatches[i].setExposedMosquitoes(E);
            mosquitoPatches[i].setInfectedMosquitoes(I);
            mosquitoPatches[i].setTemp(T);
        }

        // 4.4 Recalculate SEI and update total layers
        int total0 = valuetotalmosquitoes0->getValueAt(pt, err);
        int total1 = valuetotalmosquitoes1->getValueAt(pt, err);
        mosquitoPatches[i].recalculateSEI(1.0, total0, total1);

        int newTotal = mosquitoPatches[i].getSuceptibleMosquitoes()
                     + mosquitoPatches[i].getInfectedMosquitoes();
        int tick = repast::RepastProcess::instance()->getScheduleRunner().currentTick();
        if (tick % 2 == 0)
            valuetotalmosquitoes1->setValueAt(newTotal, pt, err);
        else
            valuetotalmosquitoes0->setValueAt(newTotal, pt, err);

        valueLayerSuceptibleMosquitoes->setValueAt(
            mosquitoPatches[i].getSuceptibleMosquitoes(), pt, err);
        valueLayerExposedMosquitoes->setValueAt(
            mosquitoPatches[i].getExposedMosquitoes(), pt, err);
        valueLayerInfectedMosquitoes->setValueAt(
            mosquitoPatches[i].getInfectedMosquitoes(), pt, err);
    }

    // 4.5 Final synchronizations
    valueLayerSuceptibleMosquitoes->synchronize();
    valueLayerExposedMosquitoes->synchronize();
    valueLayerInfectedMosquitoes->synchronize();
    valueLayerTemperature->synchronize();
    valueLayerType->synchronize();
    valuetotalmosquitoes0->synchronize();
    valuetotalmosquitoes1->synchronize();
}

void RepastHPCModel::initRecords(){
	int rank = repast::RepastProcess::instance()->rank();
	if(rank == 0){
		writer.initCSVFile();
		/* std::ofstream ofs("test.csv", std::ofstream::trunc);
		ofs << "tick"<< ","<<"TotalHumans"<< ","<<"Susceptible"<<","<<"Exposed"<<","<<"Infected"<<","<<"Recovered"<<"\n";
		ofs.close(); */

        // Initialize weekly neighborhood CSV
        std::ofstream ofs2("weekly_neighborhood_cases.csv", std::ofstream::trunc);
        ofs2 << "week,new_cases_oasis,new_cases_laquinina,new_cases_ondasdelcaribe,"
             << "new_cases_pantano,new_cases_ochodediciembre,new_cases_lacoquera,"
             << "new_cases_sanjacinto,new_cases_nuevabethel,"
             << "total_humans_oasis,total_humans_laquinina,total_humans_ondasdelcaribe,"
             << "total_humans_pantano,total_humans_ochodediciembre,total_humans_lacoquera,"
             << "total_humans_sanjacinto,total_humans_nuevabethel\n";
        ofs2.close();
	}
}
void RepastHPCModel::recordResults(){
	int rank = repast::RepastProcess::instance()->rank();
	if(rank == 0){
		int* totals = countHumansInState();
		writer.appendRecord(totals);
		/* std::ofstream outfile;
		outfile.open("test.csv", std::ios_base::app); // append instead of overwrite
		outfile << tick << ","<<countOfHumans<< ","<<totals[0]<<","<<totals[1]<<","<<totals[2]<<","<<totals[3]<<"\n"; */

        // Every 7 ticks (excluding tick 0), write weekly neighborhood totals and reset
        int tick = repast::RepastProcess::instance()->getScheduleRunner().currentTick();
        if (tick != 0 && (tick % 7) == 0) {
            // Count total humans in each neighborhood
            total_humans_oasis = 0;
            total_humans_laquinina = 0;
            total_humans_ondasdelcaribe = 0;
            total_humans_pantano = 0;
            total_humans_ochodediciembre = 0;
            total_humans_lacoquera = 0;
            total_humans_sanjacinto = 0;
            total_humans_nuevabethel = 0;
            
            std::vector<Human*> humans;
            context.selectAgents(countOfHumans, humans);
            for (Human* h : humans) {
                std::vector<int> loc;
                discreteSpace->getLocation(h->getId(), loc);
                repast::Point<int> pt(loc[0], loc[1]);
                bool err = false;
                int nid = valueLayerNeighborhood->getValueAt(pt, err);
                if (!err) {
                    if (nid == 1) ++total_humans_oasis;
                    else if (nid == 2) ++total_humans_laquinina;
                    else if (nid == 3) ++total_humans_ondasdelcaribe;
                    else if (nid == 4) ++total_humans_pantano;
                    else if (nid == 5) ++total_humans_ochodediciembre;
                    else if (nid == 6) ++total_humans_lacoquera;
                    else if (nid == 7) ++total_humans_sanjacinto;
                    else if (nid == 8) ++total_humans_nuevabethel;
                }
            }
            
            std::ofstream ofs("weekly_neighborhood_cases.csv", std::ios_base::app);
            int week = tick / 7;
            ofs << week << "," 
                << new_cases_oasis << "," << new_cases_laquinina << ","
                << new_cases_ondasdelcaribe << "," << new_cases_pantano << ","
                << new_cases_ochodediciembre << "," << new_cases_lacoquera << ","
                << new_cases_sanjacinto << "," << new_cases_nuevabethel << ","
                << total_humans_oasis << "," << total_humans_laquinina << ","
                << total_humans_ondasdelcaribe << "," << total_humans_pantano << ","
                << total_humans_ochodediciembre << "," << total_humans_lacoquera << ","
                << total_humans_sanjacinto << "," << total_humans_nuevabethel << "\n";
            ofs.close();
            
            // Reset all counters
            new_cases_oasis = 0;
            new_cases_laquinina = 0;
            new_cases_ondasdelcaribe = 0;
            new_cases_pantano = 0;
            new_cases_ochodediciembre = 0;
            new_cases_lacoquera = 0;
            new_cases_sanjacinto = 0;
            new_cases_nuevabethel = 0;
        }
	}
}

//aca es donde se ejecuta toda la simualcion con la ayuda del runner.scheduleEvent
void RepastHPCModel::initSchedule(repast::ScheduleRunner& runner){
	runner.scheduleEvent(1, 1, repast::Schedule::FunctorPtr(new repast::MethodFunctor<RepastHPCModel> (this, &RepastHPCModel::runAllPatches))); //aca lo que hace esto es llamar al metodo runAllPatches cada 1 tick y desde el primer tick. Este metodo recibe (startingTime, cadaCuantoSeReptie...)
	runner.scheduleEvent(1, 1, repast::Schedule::FunctorPtr(new repast::MethodFunctor<RepastHPCModel> (this, &RepastHPCModel::runAllHumans))); //aca lo que hace esto es llamar al metodo runAllHumans cada 1 tick y desde el primer tick. Este metodo recibe (startingTime, cadaCuantoSeReptie...) */
	runner.scheduleEvent(1, 1, repast::Schedule::FunctorPtr(new repast::MethodFunctor<RepastHPCModel> (this, &RepastHPCModel::recordResults))); //this event is schedulled to run after the simulation arrives at the stop end
	// Schedule printExecutionTime() 1 tick before the end
    runner.scheduleEvent(stopAt - 1, 1,  // Start at (stopAt-1), interval=1 (runs once)
        repast::Schedule::FunctorPtr(new repast::MethodFunctor<RepastHPCModel>
        (this, &RepastHPCModel::printExecutionTime)));
    runner.scheduleStop(stopAt);
}

int* RepastHPCModel::countHumansInState(){
	int *totals = new int [7];

	int current_tick = repast::RepastProcess::instance()->getScheduleRunner().currentTick();;
	int total_humans = countOfHumans;
	int total_susceptible = 0;
	int total_exposed = 0;
	int total_infected= 0;
	int total_recovered = 0;
	int new_cases = 0;

	std::vector<Human*> humans;
	context.selectAgents(countOfHumans, humans);//cojer todos humanos en el contexto y meterlos  en el vector de arriba
	std::vector<Human*>::iterator it=humans.begin(); //crear un iterador que itere sobre este vector y comienze con el primer pointer
	while(it!=humans.end()){
		//cojamos el ID
		repast::AgentId id=(*it)->getId();
		string human_infection_state = (*it)->getInfectionState();
		if (human_infection_state == "susceptible"){
			total_susceptible = total_susceptible + 1;
		}
		if (human_infection_state == "exposed"){
			total_exposed = total_exposed + 1;
			
			
		}
		if (human_infection_state == "infected"){
			total_infected= total_infected + 1;
			if ((*it)->isNewlyInfected()) {
                new_cases++;
				(*it)->resetNewlyInfected();
			}
		}	
		if (human_infection_state == "recovered"){
			total_recovered= total_recovered + 1;
		}	
		
		it++;
	}
	totals[0] = current_tick;
	totals[1] = total_humans;
	totals[2] = total_susceptible;
	totals[3] = total_exposed ;
	totals[4] = total_infected;
	totals[5] = total_recovered;
	totals[6] = new_cases;
	return totals;
}


//metodo auxiliar para la inizialicion
int RepastHPCModel::ageInitializer(double prob) {
	if(0<=prob && prob<0.0811){ //los numeros de distribuciones se crean asi
		return repast::Random::instance()->createUniIntGenerator(0, 4).next();
	}
	if (0.0811 <= prob && prob < 0.1615) { //tambien se pueden crear en un solo paso 
		return repast::Random::instance()->createUniIntGenerator(5, 9).next();
	}
	if (0.1615 <= prob && prob < 0.2435) {
		return repast::Random::instance()->createUniIntGenerator(10, 14).next();
	}
	 if (0.2435 <= prob && prob < 0.3292) {
		return repast::Random::instance()->createUniIntGenerator(15, 19).next();
	 }
	if (0.3292 <= prob && prob < 0.4249) {
		return repast::Random::instance()->createUniIntGenerator(20, 24).next();
	}
	if (0.4249 <= prob && prob < 0.5196) {
		return repast::Random::instance()->createUniIntGenerator(25, 29).next();
	}
	if (0.5196 <= prob && prob < 0.6025) {
		return repast::Random::instance()->createUniIntGenerator(30, 34).next();
	}
	if (0.6025 <= prob && prob < 0.6801) {
		return repast::Random::instance()->createUniIntGenerator(35, 39).next();
	}
	if (0.6801 <= prob && prob < 0.7506) {
		return repast::Random::instance()->createUniIntGenerator(40,44).next();
	} 
	if (0.7506 <= prob && prob < 0.8097) {
		return repast::Random::instance()->createUniIntGenerator(45, 49).next();
	}
	if (0.8097 <= prob && prob < 0.8626) {
		return repast::Random::instance()->createUniIntGenerator(50, 54).next();
	} 
	if (0.8626 <= prob && prob < 0.9062) {
		return repast::Random::instance()->createUniIntGenerator(55, 59).next();
	} 
	if (0.9062 <= prob && prob < 0.9378) {
		return repast::Random::instance()->createUniIntGenerator(60, 64).next();

	} 
	if (0.9378 <= prob && prob < 0.9601) {
		return repast::Random::instance()->createUniIntGenerator(65, 69).next();
	}
	if (0.9601 <= prob && prob < 0.9768) {
		return repast::Random::instance()->createUniIntGenerator(70, 74).next();
	} 
	if (0.9768 <= prob && prob < 0.9890) {
		return repast::Random::instance()->createUniIntGenerator(75, 79).next();
	}
    if (0.9890 <= prob && prob <= 1) {
        return repast::Random::instance()->createUniIntGenerator(80, 90).next();
    }
    // Fallback to a reasonable age if rounding errors occur
    return 30;
} 

//metodo auxiliar para la inizialicion
std::vector<int> RepastHPCModel::homeLocationInitializer(){
	// int xHome=repast::Random::instance()->createUniIntGenerator(0,xdim/2).next();
	// int yHome=repast::Random::instance()->createUniIntGenerator(0,ydim/2).next();
	// std::vector<int> home={xHome,yHome};

	std::vector<int> home;
	home = placeSelector(1);
	
	return home;
}

//metodo auxiliar para la inizialicion
std::vector<std::vector<int>> RepastHPCModel::activitiesInitializer(int age) {
    // home is handled by homeLocationInitializer(), so we only build two lists:
    
    // 1a) study or work
    std::vector<int> studyOrWork;
    if (age <= 24) {
        // under 24 → school
        studyOrWork = placeSelector(2);
    }
    else if (age <= 65) {
        // 24–65 → work
        studyOrWork = placeSelector(3);
    }
    else {
        // over 65 → we treat as “other activities”
        studyOrWork = placeSelector(4);
    }

    // 1b) other activities (leisure)
    std::vector<int> other = placeSelector(4);

    return { studyOrWork, other };
}

std::vector<int> RepastHPCModel::placeSelector(int location_value) {
    // 1) Compute raw‐coordinate bounds
    double minRawX = patchCenters[0].first;
    double maxRawX = patchCenters[0].first;
    double minRawY = patchCenters[0].second;
    double maxRawY = patchCenters[0].second;
    for (auto &p : patchCenters) {
        minRawX = std::min(minRawX, p.first);
        maxRawX = std::max(maxRawX, p.first);
        minRawY = std::min(minRawY, p.second);
        maxRawY = std::max(maxRawY, p.second);
    }

    // 2) Prepare random index generator over all patches
    // Use Repast RNG for reproducibility across ranks and runs
    auto& rng = *repast::Random::instance();
    auto uni_idx = rng.createUniIntGenerator(0, (int)patchCenters.size() - 1);

    std::vector<int> place(2);
    int act_location = -1;
    bool errorFlag = false;

    // 3) Loop until we find a patch whose type matches location_value
    while (act_location != location_value) {
        size_t idx = static_cast<size_t>(uni_idx.next());
        double rawX = patchCenters[idx].first;
        double rawY = patchCenters[idx].second;

        // 4) Normalize into [0,1]
        double normX = (rawX - minRawX) / (maxRawX - minRawX);
        double normY = (rawY - minRawY) / (maxRawY - minRawY);

        // 5) Scale into your discrete grid [0 .. xdim-1], [0 .. ydim-1]
        int coordX = static_cast<int>(std::round(normX * (xdim - 1)));
        int coordY = static_cast<int>(std::round(normY * (ydim - 1)));

        place = { coordX, coordY };

        // 6) Read the patch type
        act_location = valueLayerType->getValueAt(place, errorFlag);
        if (errorFlag) {
            std::cerr << "Error reading type at (" 
                      << coordX << "," << coordY << ")" << std::endl;
            // You might choose to break or continue here
        }
    }

    return place;
}

void RepastHPCModel::printExecutionTime() {
    int rank = repast::RepastProcess::instance()->rank();
    std::cout << "Rank " << rank << ": Entered printExecutionTime()" << std::endl; // Debug line
    
    if (rank == 0) {
        end_time = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        std::cout << "\n=== Total Execution Time: " << duration.count() << " seconds ===" << std::endl;
    }
}