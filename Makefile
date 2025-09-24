include ./env
SHELL=/bin/bash

n ?= 1

olduser ?= user

newuser ?= msuribec

all: Main.exe run

#create an executable named Main every time that Main.o, MyPatch.o or Human.o change. But where do these .o files come from?

Main.exe: Main.o MyHuman.o SEIModel.o MyReadData.o CSVWriter.o MyModel.o
	$(MPICXX) -o $@ $^ $(BOOST_LIB_DIR) $(REPAST_HPC_LIB_DIR) $(GDAL_LIB_DIR) $(REPAST_HPC_LIB) $(BOOST_LIBS) $(GDAL_LIB)
# #the line above recompile the files Main.o and Rectangle.o into a single executable named Main

#create a Main.o file  everytime Main.cpp changes
Main.o: Main.cpp
	$(MPICXX) $(REPAST_HPC_DEFINES) $(BOOST_INCLUDE) $(REPAST_HPC_INCLUDE) $(GDAL_INCLUDE) -I./include -c ./Main.cpp -o ./Main.o

#create a MyPatch.o file  everytime MyPatch.cpp or MyPatch.h changes
SEIModel.o: SEIModel.cpp SEIModel.h
	$(MPICXX) $(REPAST_HPC_DEFINES) $(BOOST_INCLUDE) $(REPAST_HPC_INCLUDE) $(GDAL_INCLUDE) -c SEIModel.cpp
	$(MPICXX) $(REPAST_HPC_DEFINES) $(BOOST_INCLUDE) $(REPAST_HPC_INCLUDE) $(GDAL_INCLUDE) -I./include -c ./SEIModel.cpp -o ./SEIModel.o

#create a Human.o file  everytime Human.cpp or Human.h changes

MyHuman.o: MyHuman.cpp MyHuman.h
	$(MPICXX) $(REPAST_HPC_DEFINES) $(BOOST_INCLUDE) $(REPAST_HPC_INCLUDE) $(GDAL_INCLUDE) -I./include -c ./MyHuman.cpp -o ./MyHuman.o


MyReadData.o: MyReadData.cpp MyReadData.h
	$(MPICXX) $(REPAST_HPC_DEFINES) $(BOOST_INCLUDE) $(REPAST_HPC_INCLUDE) $(GDAL_INCLUDE) -I./include -c ./MyReadData.cpp -o ./MyReadData.o 

CSVWriter.o: CSVWriter.cpp CSVWriter.h
	$(MPICXX) $(REPAST_HPC_DEFINES) $(BOOST_INCLUDE) $(REPAST_HPC_INCLUDE) $(GDAL_INCLUDE) -I./include -c ./CSVWriter.cpp -o ./CSVWriter.o
	
#create a Model.o file  everytime Model.cpp or Model.h changes
MyModel.o: MyModel.cpp MyModel.h
	$(MPICXX) $(REPAST_HPC_DEFINES) $(BOOST_INCLUDE) $(REPAST_HPC_INCLUDE) $(GDAL_INCLUDE) -I./include -c ./MyModel.cpp -o ./MyModel.o

#remove anything that has a .o and also the Main executable
clean:
	rm -f *.o
	rm -f Main.exe

export OMPI_MCA_btl_vader_single_copy_mechanism=none

.ONESHELL:
modifyenv:
	sed -i 's/$(olduser)/$(newuser)/g' env

.ONESHELL:
run:
	mpirun -n $(n) ./Main.exe ./props/config.props ./props/model.props