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
	
	// Parse parameters from model.props
	repast::Properties props(propsFile, argc, argv, &world);
	Params params;
	
	// Read parameters with defaults
	params.sigma_M = (props.getProperty("sigma_M") != "") ? repast::strToDouble(props.getProperty("sigma_M")) : 0.3;
	params.sigma_H = (props.getProperty("sigma_H") != "") ? repast::strToDouble(props.getProperty("sigma_H")) : 3.0;
	params.z = (props.getProperty("z") != "") ? repast::strToDouble(props.getProperty("z")) : 0.3;
	params.r = (props.getProperty("r") != "") ? repast::strToDouble(props.getProperty("r")) : 0.6;
	params.C = (props.getProperty("C") != "") ? repast::strToDouble(props.getProperty("C")) : 30.0;
	params.beta_mh = (props.getProperty("beta_mh") != "") ? repast::strToDouble(props.getProperty("beta_mh")) : 0.10;
	params.beta_hm = (props.getProperty("beta_hm") != "") ? repast::strToDouble(props.getProperty("beta_hm")) : 0.10;
	params.base_seed = (props.getProperty("base_seed") != "") ? repast::strToInt(props.getProperty("base_seed")) : 12345;
	params.replicate_id = (props.getProperty("replicate_id") != "") ? repast::strToInt(props.getProperty("replicate_id")) : 0;
	params.config_id = (props.getProperty("config_id") != "") ? props.getProperty("config_id") : "baseline";
	params.perturb_param = (props.getProperty("perturb_param") != "") ? props.getProperty("perturb_param") : "baseline";
	params.perturb_delta = (props.getProperty("perturb_delta") != "") ? repast::strToDouble(props.getProperty("perturb_delta")) : 0.0;
	params.obs_csv = (props.getProperty("obs_csv") != "") ? props.getProperty("obs_csv") : "";
	
	// Apply perturbation if specified
	if (params.perturb_param != "baseline" && params.perturb_delta != 0.0) {
		if (params.perturb_param == "sigma_M") params.sigma_M *= (1.0 + params.perturb_delta);
		else if (params.perturb_param == "sigma_H") params.sigma_H *= (1.0 + params.perturb_delta);
		else if (params.perturb_param == "z") params.z *= (1.0 + params.perturb_delta);
		else if (params.perturb_param == "r") params.r *= (1.0 + params.perturb_delta);
		else if (params.perturb_param == "C") params.C *= (1.0 + params.perturb_delta);
		else if (params.perturb_param == "beta_mh") params.beta_mh *= (1.0 + params.perturb_delta);
		else if (params.perturb_param == "beta_hm") params.beta_hm *= (1.0 + params.perturb_delta);
	}
		
	// crear el modelo
	RepastHPCModel* model = new RepastHPCModel(propsFile, argc, argv, &world, params); //Crear un modelo de repast
	//crear una instancia de un repast schedule runner 
	repast::ScheduleRunner& runner = repast::RepastProcess::instance()->getScheduleRunner(); //se uriliza el repast ptocess creado arriba
	
	//con estos tres llamados se corre la simulacion

	model->init(); //se inizializan los humanos y los grid value layers
	model->initSchedule(runner);
	runner.run();
	std::cout << "Fin de la simulacion " << endl;
	
	// Delete model BEFORE calling RepastProcess::done()
	delete model;
	
	// Clean up Repast process last
	repast::RepastProcess::instance()->done();
}
