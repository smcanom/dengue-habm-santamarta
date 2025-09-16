
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <iterator>
#include <string>
#include <algorithm>
#include "CSVWriter.h"
//constructor that takes arguments 

CSVWriter::CSVWriter(std::string Filename,std::string Delimiter){
    fileName = Filename;
	delimeter = Delimiter;
}

void CSVWriter::initCSVFile(){
	std::ofstream ofs(fileName, std::ofstream::trunc);
	ofs << "tick"<< ","<<"TotalHumans"<< ","<<"Susceptible"<<","<<"Exposed"<<","<<"Infected"<<","<<"Recovered"<<","<<"New Cases per day"<<"\n";
	ofs.close();
}

void CSVWriter::appendRecord(int* totals){
	std::ofstream outfile;
	outfile.open(fileName, std::ios_base::app); // append instead of overwrite
	outfile << totals[0] << ","<<totals[1]<< ","<<totals[2]<<","<<totals[3]<<","<<totals[4]<<","<<totals[5]<<","<<totals[6]<<"\n";
}
