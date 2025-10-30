#include "SEIModel.h"
#include "MyModel.h"  // for Params
#include <random>
#include <iostream>
#include <vector>
#include <deque>
#include <cmath>
#include "repast_hpc/AgentId.h"
#include "repast_hpc/RepastProcess.h"
#include "repast_hpc/Utilities.h"
#include "repast_hpc/Properties.h"
#include "repast_hpc/initialize_random.h"
using namespace std;

//constructor
SEIModel::SEIModel(double suceptibleM, double exposedM, double infectedM, double temp,
                   int delayTicks, double femaleFrac, double omega, const Params* params)
  : suceptibleMosquitoes(std::max(0.0, suceptibleM)),
    exposedMosquitoes(std::max(0.0, exposedM)),
    infectedMosquitoes(std::max(0.0, infectedM)),
    temperature(std::max(0.0, temp)),
    infectedHumans(0),
    nHumans(0),
    deathRate(0),
    maturationDelay(std::max(1, delayTicks)),
    femaleProportion(std::max(0.0, std::min(1.0, femaleFrac))),
    emergenceMultiplier(std::max(0.0, omega)),
    params_(params)
{
  // Calcula población inicial de adultos
  double N0 = suceptibleMosquitoes + exposedMosquitoes + infectedMosquitoes;

  // Tasa de puesta en t = 0 (misma fórmula que en recalculateSEI)
  // Use r from params if available, else fall back to 0.2
  double reproductionRate = params_ ? params_->r : 0.6;
  double carryingCapacity = params_ ? params_->C : 20.0;
  double At0 = reproductionRate * N0 *
               std::exp(1.0 - (N0 / carryingCapacity));

  // Evita valores negativos por si acaso
  At0 = std::max(0.0, At0);

  // Inicializa la cola de huevos con ese valor "biológico"
  eggQueue = std::deque<double>(maturationDelay, At0);
}



//getters
double SEIModel::getSuceptibleMosquitoes(){
    return suceptibleMosquitoes;
}

double SEIModel::getExposedMosquitoes(){
    return exposedMosquitoes;
}

double SEIModel::getInfectedMosquitoes(){
    return infectedMosquitoes;
}

double SEIModel::getTemp(){
    return temperature;
}

//setters(a setter for type is not required because patched dont change type)
void SEIModel::setSuceptibleMosquitoes(double suceptibleM){
    suceptibleMosquitoes=std::max(0.0, suceptibleM);
}

void SEIModel::setExposedMosquitoes(double exposedM){
    exposedMosquitoes=std::max(0.0, exposedM);
}

void SEIModel::setInfectedMosquitoes(double infectedM){
    infectedMosquitoes=std::max(0.0, infectedM);
}

void SEIModel::setTemp(double temp){
    temperature=std::max(0.0, temp);
}

void SEIModel::setHumans(int humans){
    nHumans = humans;
}

void SEIModel::setInfectedHumans(int infected){
    infectedHumans = infected;
}

void SEIModel::setEmergenceMultiplier(double omega){
    emergenceMultiplier = std::max(0.0, omega);
}

void SEIModel::recalculateSEI(double h, double totalM0, double totalM1) {
    try {
        // 1. Compute death rate
        calculateDeathRate();
        
        // 2. Get current state
        double S = std::max(0.0, suceptibleMosquitoes);
        double I = std::max(0.0, infectedMosquitoes);
        double N = std::max(1.0, S + I);  // Prevent division by zero

        // 3. Get lagged total for oviposition
        int tick = repast::RepastProcess::instance()->getScheduleRunner().currentTick();
        double N_lag = (tick % 2 == 0) ? std::max(0.0, totalM0) : std::max(0.0, totalM1);

        // 4. Calculate density-dependent egg laying
        // Use r and C from params if available
        double reproductionRate = params_ ? params_->r : 0.6;
        double carryingCapacity = params_ ? params_->C : 20.0;

        double At = reproductionRate * N_lag * std::exp(1.0 - (N / carryingCapacity));

        // 5. Update egg queue
        double eggsMaturing = 0.0;
        if (eggQueue.size() >= static_cast<size_t>(maturationDelay)) {
        eggsMaturing = eggQueue.front();      // grab the cohort laid delayTicks ago
        eggQueue.pop_front();                 // remove it from the buffer
        }
        // 6. Enqueue today's eggs (female fraction)
        double At_fem = At * femaleProportion;
        eggQueue.push_back(std::max(0.0, At_fem));

        // 7. Apply survival rates to the correctly‐aged cohort
        const double larvalSurvival = 0.9;
        const double pupalSurvival = 0.98;
        double phi = larvalSurvival * pupalSurvival;
        double newAdults = eggsMaturing * phi * emergenceMultiplier;

        // 8. Calculate human-mosquito interaction using z from params
        double z_param = params_ ? params_->z : 0.3;
        double Pv = 1.0 - std::pow((S / N), z_param);

        // 9. Update mosquito populations using beta_hm from params
        double beta_hm = params_ ? params_->beta_hm : 0.1;
        double S_next = S * (1.0 - beta_hm * Pv) * (1.0 - deathRate) + newAdults;
        double I_next = S * beta_hm * Pv * (1.0 - deathRate) + I * (1.0 - deathRate);

        // 10. Ensure non-negative values
        suceptibleMosquitoes = int(std::max(0.0, S_next));
        infectedMosquitoes = int(std::max(0.0, I_next));
    } catch (const std::exception& e) {
        std::cerr << "Error in recalculateSEI: " << e.what() << std::endl;
        // Set safe default values
        suceptibleMosquitoes = 0.0;
        infectedMosquitoes = 0.0;
    }
}


void SEIModel::actualizeTemp(int tick,
    const std::vector<std::vector<int>>& dataTemperatures) {
  try {
    // 1. Validate the data structure
    if (dataTemperatures.size() < 2 ||
        tick < 0 ||
        tick >= static_cast<int>(dataTemperatures[0].size()) ||
        tick >= static_cast<int>(dataTemperatures[1].size())) {
      std::cerr << "Invalid tick or temperature data in actualizeTemp: "
                << tick << std::endl;
      return;
    }

    // 2. Read min/max and ensure correct ordering
    int maxTemp = dataTemperatures[0][tick];
    int minTemp = dataTemperatures[1][tick];
    if (minTemp > maxTemp) {
      std::swap(minTemp, maxTemp);
    }

    // 3. Draw a new temperature (continuous)
    double newTemp = repast::Random::instance()
                        ->createUniDoubleGenerator(minTemp, maxTemp)
                        .next();

    // 4. Update the model’s temperature
    this->setTemp(newTemp);

  } catch (const std::exception& e) {
    std::cerr << "Error in actualizeTemp: " << e.what() << std::endl;
    setTemp(25.0);  // Safe fallback
  }
}


//auxiliaryMethods
double SEIModel::suceptible_function(double suceptible,double exposed,double infected,double birthRate, double infectionRate) {
    double s1=birthRate-infectionRate*suceptible-deathRate*suceptible;
	return s1;
}
	
double SEIModel::exposed_function(double suceptible,double exposed,double infected,double infectionRate) {
	double exposedToinfectedRate=calculateExposedToinfectedRate();
	double e1=infectionRate*suceptible-exposedToinfectedRate*exposed-deathRate*exposed;
	return e1;
}
	
double SEIModel::infected_function(double suceptible,double exposed,double infected) {
    double exposedToinfectedRate=calculateExposedToinfectedRate();
    double i1=exposedToinfectedRate*exposed-deathRate*infected;
    return i1;
}

double SEIModel::calculateBirthRate() {
    double totalMosquitoes=suceptibleMosquitoes+infectedMosquitoes+exposedMosquitoes;
    double mosquitoPopulationGrowthRate=naturalEmergenceRate-deathRate;
    double carryingCapacity = params_ ? params_->C : 20.0;
    double birthRate=totalMosquitoes*(naturalEmergenceRate-mosquitoPopulationGrowthRate*totalMosquitoes/carryingCapacity);
    return birthRate;
}

double SEIModel::calculateDeathRate() {
    try {
        double temp = std::max(0.0, getTemp());
        const double b0 = 0.8692;
        const double b1 = -0.1590;
        const double b2 = 0.01116;
        const double b3 = -0.0003408;
        const double b4 = 0.000003809;

        double sumPolGrade4 = b0 + b1 * temp + b2 * std::pow(temp, 2) + 
                             b3 * std::pow(temp, 3) + b4 * std::pow(temp, 4);

        deathRate = std::max(0.0, std::min(1.0, sumPolGrade4));
        return deathRate;
    } catch (const std::exception& e) {
        std::cerr << "Error in calculateDeathRate: " << e.what() << std::endl;
        deathRate = 0.1;  // Safe default death rate
        return deathRate;
    }
}
   
double SEIModel::calculateInfectionRate() { //depende del numero de humanos en el patch   
    int totalHumans=calculateTotalHumansInPatch();
    int humansInfected=calculateInfectedHumansInPatch();
    double totalMosquitoes=suceptibleMosquitoes+infectedMosquitoes+exposedMosquitoes;
    
    double totalSuccesfulBites=(mosquitoBiteDemand*totalMosquitoes*maxBitesPerHuman*totalHumans)/(mosquitoBiteDemand*totalMosquitoes+maxBitesPerHuman*totalHumans);
    double successfulBitesPerMosquito=totalSuccesfulBites/totalMosquitoes;
    double infectionRateMosquitoes=0;
    if (totalHumans>0) {
        double beta_hm = params_ ? params_->beta_hm : 0.1;
        infectionRateMosquitoes=successfulBitesPerMosquito*beta_hm*((humansInfected+0.0)/totalHumans);
    }
    return infectionRateMosquitoes;
}

double SEIModel::calculateExposedToinfectedRate() {
    double exposedToInfectedRate;
    if (temperature<15) {
        exposedToInfectedRate=0;
    }
    else {
        double patchIncubationPeriod = 0;
        if (15<temperature && temperature<21) {
            std::uniform_real_distribution<double> myUnifDist(10,25);//generate uniform distribution
            std::default_random_engine re;//generate random number generator
            re.seed(std::random_device{}()); //generate seed of random number generator 
            patchIncubationPeriod=myUnifDist(re);
        } else if (21<=temperature && temperature<26) {
            std::uniform_real_distribution<double> myUnifDist(7,10);//generate uniform distribution
            std::default_random_engine re;//generate random number generator
            re.seed(std::random_device{}()); //generate seed of random number generator 
            patchIncubationPeriod=myUnifDist(re);
        } else if (26<=temperature && temperature<31) {
            std::uniform_real_distribution<double> myUnifDist(4,7);//generate uniform distribution
            std::default_random_engine re;//generate random number generator
            re.seed(std::random_device{}()); //generate seed of random number generator 
            patchIncubationPeriod=myUnifDist(re);
        }
        //patchIncubationPeriod=10;
        exposedToInfectedRate=1/patchIncubationPeriod;
    }
    return exposedToInfectedRate;
}

int SEIModel::calculateTotalHumansInPatch() {	// TEMPORAL->se deben contar la cantidad de humanos en el patch
    return nHumans;
}

int SEIModel::calculateInfectedHumansInPatch() { // TEMPORAL->se deben contar la cantidad de humanos infectados en el patch
    return infectedHumans;
}	
