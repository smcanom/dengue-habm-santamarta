#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <random>
#include <vector>
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

//OBLIGATORY (NO SE MUY BIEN PARA QUE ES)
HumanPackageReceiver::HumanPackageReceiver(repast::SharedContext<Human>* agentPtr): agents(agentPtr){}

//OBLIGATORY (NO SE MUY BIEN PARA QUE ES)
Human * HumanPackageReceiver::createAgent(HumanPackage package){
    repast::AgentId id(package.id, package.rank, package.type, package.currentRank);
    return new Human(id, package.infectionState, package.age, package.timeSinceSuccesfullBite, package.timeSinceInfection, package.homeLocation, package.activities);

}

//OBLIGATORY (NO SE MUY BIEN PARA QUE ES)
void HumanPackageReceiver::updateAgent(HumanPackage package){
    repast::AgentId id(package.id, package.rank, package.type,package.currentRank);//creo un id
    Human * agent = agents->getAgent(id);//cojo el agente que tiene ese id
    agent->setAll(package.infectionState, package.age, package.timeSinceSuccesfullBite, package.timeSinceInfection, package.homeLocation, package.activities);//le cambio las variables de estado
}

//CONSTRUCTOR 
RepastHPCModel::RepastHPCModel(std::string propsFile, int argc, char** argv, boost::mpi::communicator* comm): context(comm){
	props = new repast::Properties(propsFile, argc, argv, comm); 
	//atributes that come from the model.props file
    start_time = std::chrono::steady_clock::now();
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

	//define process dimensions 
	std::vector<int> processDims;
	processDims.push_back(2);//there will be 2 proceses in the x dimension
	processDims.push_back(2);//there will be 2 proceses in the y dimension
	//this indicates that we will be runnning the simulation with 4 proceses

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

	//These value layers store the total mosquitoes on t-2.
	valuetotalmosquitoes0 = new repast::ValueLayerND<int>(processDims, gd, bufferSize, false);
	valuetotalmosquitoes1  = new repast::ValueLayerND<int>(processDims, gd, bufferSize, false);

    // Initialize mosquitoPatches as a flat array:
    mosquitoPatches.resize(patchCenters.size());
    for (size_t i = 0; i < patchCenters.size(); ++i) {
        mosquitoPatches[i] = SEIModel(0,0,0,0, /*delay=*/11, /*rho=*/0.5);
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
	receiver=new HumanPackageReceiver(&context);
	
}

RepastHPCModel::~RepastHPCModel() {
    // Delete properties
    if (props) {
        delete props;
        props = nullptr;
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
    if (valuetotalmosquitoes0) {
        delete valuetotalmosquitoes0;
        valuetotalmosquitoes0 = nullptr;
    }
    if (valuetotalmosquitoes1) {
        delete valuetotalmosquitoes1;
        valuetotalmosquitoes1 = nullptr;
    }

    // Delete agent communication helpers
    if (provider) {
        delete provider;
        provider = nullptr;
    }
    if (receiver) {
        delete receiver;
        receiver = nullptr;
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
}
 void RepastHPCModel::init(){
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
			Human* H=new Human(id,infectionState,age,timeSinceSuccesfullBite,timeSinceInfection,homeLocation,activities);
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

        int rand_susc = triangularSampleInt(0,10,50);
        int rand_exp  = 0;
        int rand_inf  = repast::Random::instance()
                          ->createUniIntGenerator(0,2).next();
        int total = rand_inf + rand_susc;
        valuetotalmosquitoes0->setValueAt(total, location, err);
        valuetotalmosquitoes1->setValueAt(total, location, err);
        // then synchronize once, so both layers are non‐zero at tick=1
        valuetotalmosquitoes0->synchronize();
        valuetotalmosquitoes1->synchronize();
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
                                      delay, rho);

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

        // c) log & move
        std::cout << "Human " << id.id()
                  << " from (" << oldLoc[0] << "," << oldLoc[1] << ")"
                  << " → ("   << dest[0]   << "," << dest[1]   << ")\n";
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
        h->actualizeSEIRStatus(&context);

        int m1 = valuetotalmosquitoes1->getValueAt(pt, err);
        int m0 = valuetotalmosquitoes0->getValueAt(pt, err);
        SEIModel patch(susc, exposed, infd, temp);
        patch.setHumans(nHumans);
        patch.setInfectedHumans(nInfected);
        patch.recalculateSEI(0.1, m0, m1);

        valueLayerSuceptibleMosquitoes->setValueAt(
          patch.getSuceptibleMosquitoes(), pt, err);
        valueLayerExposedMosquitoes->setValueAt(
          patch.getExposedMosquitoes(),    pt, err);
        valueLayerInfectedMosquitoes->setValueAt(
          patch.getInfectedMosquitoes(),   pt, err);

        if (phase == 2) {
          // final phase: update time counters on the human
          h->actualizeTimes();
        }
      }
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
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> dist(0, patchCenters.size() - 1);

    std::vector<int> place(2);
    int act_location = -1;
    bool errorFlag = false;

    // 3) Loop until we find a patch whose type matches location_value
    while (act_location != location_value) {
        size_t idx = dist(gen);
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