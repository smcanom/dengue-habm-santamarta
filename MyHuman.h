
#ifndef MyHuman
#define MyHuman
#include <string>
#include <iostream>
#include <tuple>
#include <vector>
#include "repast_hpc/AgentId.h"
#include "repast_hpc/SharedContext.h" //this is the collection in which the agent will exist during the simulation runs 
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/weibull.hpp>
#include <boost/math/distributions/lognormal.hpp>
#include <random>

using namespace std;
using namespace repast;
class Human{
    private:
        repast::AgentId id_ ; //obligatory. It is made up of: [IdinTheProcess,starting rank,type, current rank]
        string infectionState;
        int age;
        int timeSinceSuccesfullBite;
        int timeSinceInfection;
        std::vector<int> homeLocation;
        std::vector<std::vector<int>> activities;
        double infectionProb;
        bool newlyInfected;
    
    public:
        //Constant variables (they are the same for each instance)
        const double naturalEmergenceRate=0.1;
        const double deathRate=0.071428571428571;
        const double mosquitoCarryingCapacity=50;
        const double mosquitoBiteDemand=20;
        const double maxBitesPerHuman= 10;
        const double probabilityOfTransmissionHToM=0.1;
        const double probabilityOfTransmissionMToH=0.1;
        
        //constructor
        Human(repast::AgentId id,string InfectionState, int Age, int TimeSinceSuccesfullBite,int TimeSinceInfection, std::vector<int> HomeLocation, std::vector<std::vector<int>> Activities);

        //required getters
        virtual repast:: AgentId& getId(){return id_;}
        virtual const repast:: AgentId& getId() const{return id_;}

        //other getters
        string getInfectionState();
        int getAge();
        int getTimeSinceSuccesfullBite();
        int getTimeSinceInfection();
        std::vector<int> getHomeLocation();
        std::vector<std::vector<int>> getActivities();

        //setters
        void setAll(std::string InfectionState,int Age,int TimeSinceSuccesfullBite,int TimeSinceInfection,std::vector<int> HomeLocation,std::vector<std::vector<int>> Activities);
        void setInfectionState(std::string InfectionState);
        void setAge(int Age);
        void setTimeSinceSuccesfullBite(int TimeSinceSuccesfullBite);
        void setTimeSinceInfection(int TimeSinceInfection);
        void setHomeLocation(std::vector<int> HomeLocation);
        void setActivities(std::vector<std::vector<int>> Activities);


        //METHODS
        void printHumanInfo();   // print all info of a human
        
        //action that humans can perform
        void actualizeTimes();
        void actualizeSEIRStatus(repast::SharedContext<Human>* context);
        
        //auxiliares
        void calculateInfectionProbabilityHuman(int infected, int susc, int nHumans);
        bool isNewlyInfected() const;
        void resetNewlyInfected();
};

/*serializable Agent package(NO ESTOY MUY SEGURA DE ESTO COMO FUNCIONA)*/
struct HumanPackage {
    public:
        //obligatory for all models
        int id;
        int rank;
        int type;
        int currentRank;

        //state variables of this model
        string infectionState;
        int age;
        int timeSinceSuccesfullBite;
        int timeSinceInfection;
        std::vector<int> homeLocation;
        std::vector<std::vector<int>> activities;

        //constructor
        HumanPackage();//for serialization
        HumanPackage(int Id,int Rank, int Type, int CurrentRank, string InfectionState, int Age, int TimeSinceSuccesfullBite,int TimeSinceInfection, std::vector<int> HomeLocation, std::vector<std::vector<int>> Activities); 

        //for archive packaging
        template<class Archive>
        void serialize(Archive &ar, const unsigned int version){
            ar & id;
            ar & rank;
            ar & type;
            ar & currentRank;
            ar & infectionState;
            ar & age;
            ar & timeSinceSuccesfullBite;
            ar & timeSinceInfection;
            ar & homeLocation;
            ar & activities;
        }

}; 
#endif