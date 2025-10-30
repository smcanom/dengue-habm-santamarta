#ifndef SEIMODEL_H
#define SEIMODEL_H
#include <vector>
#include <deque>

// Forward declaration
struct Params;

 class SEIModel{
    private:
        double suceptibleMosquitoes;
        double exposedMosquitoes;
        double infectedMosquitoes;
        double temperature;
        double infectedHumans;
        double nHumans;
        double deathRate;
        // developmental delay: a FIFO buffer holding eggs in each of the τ immature stages 
        std::deque<double> eggQueue;
        int maturationDelay;         // τ, in ticks
        double femaleProportion;     // ρ, fraction of emergent adults that are female
        double emergenceMultiplier;  // ω_i, fixed per-patch multiplier from HI
        const Params* params_;       // pointer to global parameters (not owned)
    
    public:
        static constexpr double naturalEmergenceRate=0.1;
        static constexpr double mosquitoBiteDemand=0.3;  // More realistic
        static constexpr double maxBitesPerHuman=3;     // More realistic
        // Parameters now obtained from params_:
        // double mosquitoCarryingCapacity; -> params_->C
        // double probabilityOfTransmissionHToM; -> params_->beta_hm
        // double z; -> params_->z
        // double probabilityOfTransmissionMToH; -> params_->beta_mh

        //METHODS
        //constructors
        SEIModel() : suceptibleMosquitoes(0), exposedMosquitoes(0), infectedMosquitoes(0), 
                    temperature(0), infectedHumans(0), nHumans(0), deathRate(0), 
                    maturationDelay(11), femaleProportion(0.46), emergenceMultiplier(1.0), params_(nullptr) {}

        SEIModel(double suceptibleM, double exposedM, double infectedM, double temp,
           int delayTicks = 18, double femaleFrac = 0.46, double omega = 1.0, const Params* params = nullptr);

        // SEIModel(double suceptibleM, double exposedM, double infectedM, double temp,
        //         repast::SharedContext<AgenteHumano>* context,
        //         repast::ValueLayer<int>* valueLayerSuceptibleMosquitoes,
        //         repast::ValueLayer<int>* valueLayerExposedMosquitoes,
        //         repast::ValueLayer<int>* valueLayerInfectedMosquitoes,
        //         repast::ValueLayer<int>* valueLayerTemperature,
        //         repast::ValueLayer<int>* valueLayerType);
       
        //getters
        double getSuceptibleMosquitoes();
        double getExposedMosquitoes();
        double getInfectedMosquitoes();
        double getTemp();
        
        //setters
        void setSuceptibleMosquitoes(double suceptibleM);
        void setExposedMosquitoes(double exposedM);
        void setInfectedMosquitoes(double infectedM);
        void setTemp(double temp);
        void setHumans(int humans);
        void setInfectedHumans(int infected);
        void setEmergenceMultiplier(double omega);
        
        //most important methods
        void recalculateSEI(double timeStep, double totalM0, double totalM1);
        void actualizeTemp(int tick, const std::vector<std::vector<int>>& dataTemperatures);
        
        //auxiliary methods
        double suceptible_function(double suceptible,double exposed,double infected,double birthRate, double infectionRate);
        double exposed_function(double suceptible,double exposed,double infected,double infectionRate);
        double infected_function(double suceptible,double exposed,double infected);
        double calculateDeathRate();
        double calculateBirthRate();
        double calculateInfectionRate();
        double calculateExposedToinfectedRate();
        int calculateTotalHumansInPatch();
        int calculateInfectedHumansInPatch();
};

#endif