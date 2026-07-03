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

# Headers pulled in (directly or transitively) by Main.cpp via MyModel.h. These MUST be
# dependencies: otherwise editing a header rebuilds only its own .o, leaving Main.o with a
# stale class layout. That is a one-definition-rule violation -- Main.cpp's `new
# RepastHPCModel` allocates the old sizeof while the constructor writes the new layout,
# corrupting the heap (confirmed via AddressSanitizer). Always rebuild Main.o on header change.
MODEL_HEADERS = MyModel.h MyHuman.h SEIModel.h MyReadData.h CSVWriter.h

#create a Main.o file  everytime Main.cpp or any model header changes
Main.o: Main.cpp $(MODEL_HEADERS)
	$(MPICXX) $(REPAST_HPC_DEFINES) $(BOOST_INCLUDE) $(REPAST_HPC_INCLUDE) $(GDAL_INCLUDE) -I./include -c ./Main.cpp -o ./Main.o

#create a MyPatch.o file  everytime MyPatch.cpp or MyPatch.h changes
SEIModel.o: SEIModel.cpp $(MODEL_HEADERS)
	$(MPICXX) $(REPAST_HPC_DEFINES) $(BOOST_INCLUDE) $(REPAST_HPC_INCLUDE) $(GDAL_INCLUDE) -I./include -c ./SEIModel.cpp -o ./SEIModel.o

#create a Human.o file  everytime Human.cpp or any model header changes

MyHuman.o: MyHuman.cpp $(MODEL_HEADERS)
	$(MPICXX) $(REPAST_HPC_DEFINES) $(BOOST_INCLUDE) $(REPAST_HPC_INCLUDE) $(GDAL_INCLUDE) -I./include -c ./MyHuman.cpp -o ./MyHuman.o


MyReadData.o: MyReadData.cpp $(MODEL_HEADERS)
	$(MPICXX) $(REPAST_HPC_DEFINES) $(BOOST_INCLUDE) $(REPAST_HPC_INCLUDE) $(GDAL_INCLUDE) -I./include -c ./MyReadData.cpp -o ./MyReadData.o

CSVWriter.o: CSVWriter.cpp $(MODEL_HEADERS)
	$(MPICXX) $(REPAST_HPC_DEFINES) $(BOOST_INCLUDE) $(REPAST_HPC_INCLUDE) $(GDAL_INCLUDE) -I./include -c ./CSVWriter.cpp -o ./CSVWriter.o

#create a Model.o file  everytime Model.cpp or any model header changes
MyModel.o: MyModel.cpp $(MODEL_HEADERS)
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
	rm -f results_gridchange.csv weekly_neighborhood_cases.csv neighborhood_cells.csv summary_neighborhood_sizes.csv
	mpirun -n $(n) ./Main.exe ./props/config.props ./props/model.props