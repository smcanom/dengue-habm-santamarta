#ifndef MYREADDATA_H
#define MYREADDATA_H

#include <iostream> 
#include <fstream>
#include <string> 
#include <vector>
#include <utility>
using namespace std; 

class ReadData{
    public:
        //first method
        std::vector<std::vector<int>> loadFromExcel(std::string fileName,int simulationTime);
        // std::vector<std::vector<std::pair<double, double>>> getCoor(const char* filePath);
        std::vector<std::vector<std::pair<double, double>>> getCoor(const std::string& filePath);
        std::vector<std::pair<double, double>>getPatchCenters(const std::string& filePath);
};

#endif // MYREADDATA_H