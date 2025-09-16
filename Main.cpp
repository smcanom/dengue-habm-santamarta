#include <boost/mpi.hpp>
#include "repast_hpc/RepastProcess.h"
#include <iostream>
#include "SEIModel.h"
#include "MyModel.h"
#include <vector>
using namespace std; 


int main(int argc, char** argv){
  //CREAR una instancia de un objeto tipe Repast HPC model 
	std::string configFile = argv[1]; // The name of the configuration file is Arg 1
	std::string propsFile  = argv[2]; // The name of the properties file is Arg 2
	boost::mpi::environment env(argc, argv); //crear un boost mpi environment
	boost::mpi::communicator world; //crear una instancia de  boost mpi communicator y llamarla world
	
	repast::RepastProcess::init(configFile); //crear un RepastProcess usanfo el metodo init. Esto siempre es obligatorio
	
		
	// crear el modelo
	RepastHPCModel* model = new RepastHPCModel(propsFile, argc, argv, &world); //Crear un modelo de repast
	//crear una instancia de un repast schedule runner 
	repast::ScheduleRunner& runner = repast::RepastProcess::instance()->getScheduleRunner(); //se uriliza el repast ptocess creado arriba
	
	//con estos tres llamados se corre la simulacion

	model->init(); //se inizializan los humanos y los grid value layers
	model->initSchedule(runner);
	runner.run();
	repast::RepastProcess::instance()->done();
	std::cout << "Fin de la simulacion " << endl;
	delete model;
}
