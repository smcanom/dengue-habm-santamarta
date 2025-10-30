#include "MyHuman.h"
#include "MyModel.h"  // for Params
#include <cmath>
#include <random>
#include <algorithm>
//constructor that takes arguments 
Human::Human(repast::AgentId id, string InfectionState, int Age, int TimeSinceSuccesfullBite, int TimeSinceInfection, 
             std::vector<int> HomeLocation, std::vector<std::vector<int>> Activities, const Params* params)
    : id_(id),
      infectionState(std::move(InfectionState)),
      age(std::max(0, Age)),
      timeSinceSuccesfullBite(TimeSinceSuccesfullBite),
      timeSinceInfection(TimeSinceInfection),
      homeLocation(std::move(HomeLocation)),
      activities(std::move(Activities)),
      infectionProb(0.0),
      newlyInfected(false),
      params_(params)
{
    // Validate activities
    if (activities.size() != 2) {
        std::cerr << "Warning: Human " << id.id() << " has incorrect number of activities" << std::endl;
        activities.resize(2, std::vector<int>(2, 0));
    }
    
    // Validate home location
    if (homeLocation.size() != 2) {
        std::cerr << "Warning: Human " << id.id() << " has invalid home location" << std::endl;
        homeLocation = {0, 0};
    }
}

//getters
string Human::getInfectionState(){
    return infectionState;
}
int Human::getAge(){
    return age;
}
int Human::getTimeSinceSuccesfullBite(){
    return timeSinceSuccesfullBite;
}
int Human::getTimeSinceInfection(){
    return timeSinceInfection;
}
std::vector<int> Human::getHomeLocation(){
    return homeLocation;
}
std::vector<std::vector<int>> Human::getActivities(){
    return activities;
}


//setters
void Human::setAll(std::string InfectionState,int Age,int TimeSinceSuccesfullBite,int TimeSinceInfection,std::vector<int> HomeLocation,std::vector<std::vector<int>> Activities){
    infectionState=InfectionState;
    age=Age;
    timeSinceSuccesfullBite=TimeSinceSuccesfullBite;
    timeSinceInfection=TimeSinceInfection;
    homeLocation=HomeLocation;
    activities=Activities;
}
void Human::setInfectionState(string InfectionState){
    infectionState=InfectionState;
}
void Human::setAge(int Age){
    age=Age;
}
void Human::setTimeSinceSuccesfullBite(int TimeSinceSuccesfullBite){
    timeSinceSuccesfullBite=TimeSinceSuccesfullBite;
}
void Human::setTimeSinceInfection(int TimeSinceInfection){
    timeSinceInfection=TimeSinceInfection;
}
void Human::setHomeLocation(std::vector<int> HomeLocation){
    homeLocation=HomeLocation;
}
void Human::setActivities(std::vector<std::vector<int>> Activities){
    activities=Activities;
}
//Print human info
void Human::printHumanInfo(){
    std::cout<<"id: "<<id_<<"age: "<<age<<", infecion state: "<<infectionState<<", timeSinceSuccesfullBite: "<<timeSinceSuccesfullBite<<", timeSinceInfection: "<<timeSinceInfection<<", homeLocation:"<<homeLocation[0]
	<<",  "<< homeLocation[1]<<", currentLoc:"<<"\n"; 
}

// nuevos métodos para obtener si un humano se acaba de infectar
bool Human::isNewlyInfected() const {
    return newlyInfected;
}
void Human::resetNewlyInfected() {
    newlyInfected = false;
}


//ACTIONS
//metodo para actualizar los tiempos de infeccion y tiempos desde mordida 
void Human::actualizeTimes(){
    try {
        if (infectionState == "exposed") {
            timeSinceSuccesfullBite++;
            // Check for progression to infected
            if (timeSinceSuccesfullBite >= 5) {  // 5-day incubation period
                infectionState = "infected";
                timeSinceInfection = 0;
            }
        }
        else if (infectionState == "infected") {
            timeSinceInfection++;
            // Check for recovery
            if (timeSinceInfection >= 7) {  // 7-day infection period
                infectionState = "recovered";
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Error in actualizeTimes: " << e.what() << std::endl;
    }
}

//Metodo para actualizar variable de estado SEIR
void Human::actualizeSEIRStatus(repast::SharedContext<Human>* context) {
    try {
        if (!context) {
            std::cerr << "Error: Null context in actualizeSEIRStatus" << std::endl;
            return;
        }

        if (infectionState == "susceptible") {
            // Check for infection
            if (repast::Random::instance()->nextDouble() < infectionProb) {
                infectionState = "exposed";
                timeSinceSuccesfullBite = 0;
                newlyInfected = true;
            }
        }
        else if (infectionState == "exposed") {
            // State progression handled by actualizeTimes()
        }
        else if (infectionState == "infected") {
            // State progression handled by actualizeTimes()
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error in actualizeSEIRStatus: " << e.what() << std::endl;
    }
} 

//auxiliares
//ESTA POR EL MOMENTO ES TEMPORAL HASTA QUE SE DEFINA EL ESPACIO 
/**
 * This function calculates the probability of getting infected for a human.
 *
 *
 * @param infected The amount of infected mosquitoes.
 * @param susc The amount of susceptible mosquitoes.
 * @param nHumans The total humans in the actual patch.
 */
void Human::calculateInfectionProbabilityHuman(int infected, int susc, int nHumans) {
    try {
        if (nHumans <= 0) {
            infectionProb = 0.0;
            return;
        }

        // Calculate probability based on mosquito counts and human density
        double totalMosquitoes = static_cast<double>(infected + susc);
        double humanDensity = static_cast<double>(nHumans);
        
        // Prevent division by zero
        if (totalMosquitoes <= 0.0 || humanDensity <= 0.0) {
            infectionProb = 0.0;
            return;
        }

        // Calculate bite probability
        double biteProb = (mosquitoBiteDemand * totalMosquitoes) / 
                         (mosquitoBiteDemand * totalMosquitoes + maxBitesPerHuman * humanDensity);
        
        // Use beta_mh from params if available
        double beta_mh = params_ ? params_->beta_mh : 0.1;
        
        // Calculate infection probability
        infectionProb = biteProb * beta_mh * 
                       (static_cast<double>(infected) / totalMosquitoes);
        
        // Ensure probability is between 0 and 1
        infectionProb = std::max(0.0, std::min(1.0, infectionProb));
        
    } catch (const std::exception& e) {
        std::cerr << "Error calculating infection probability: " << e.what() << std::endl;
        infectionProb = 0.0;
    }
}


/* serializable Agent package Data*/
HumanPackage::HumanPackage(){ }
HumanPackage::HumanPackage(int Id,int Rank, int Type, int CurrentRank, string InfectionState, int Age, int TimeSinceSuccesfullBite,int TimeSinceInfection, std::vector<int> HomeLocation, std::vector<std::vector<int>> Activities){
    id=Id;
    rank=Rank;
    type=Type;
    currentRank=CurrentRank;
    infectionState=InfectionState;
    age=Age;
    timeSinceSuccesfullBite=TimeSinceSuccesfullBite;
    timeSinceInfection=TimeSinceInfection;
    homeLocation=HomeLocation;
    activities=Activities;
} 