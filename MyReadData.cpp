#include "MyReadData.h"
#include <iostream> 
#include <fstream>
#include <string> 
#include <sstream>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <numeric>  

using namespace std;

std::vector<std::vector<int>> ReadData::loadFromExcel(std::string fileName, int simulationTime) {
    std::vector<std::vector<int>> allData;
    std::vector<int> list_tempMax;
    std::vector<int> list_tempMin;
    
    std::ifstream myFile;
    try {
        myFile.open(fileName);
        if (!myFile.is_open()) {
            std::cerr << "Couldn't find file: " << fileName << std::endl;
            list_tempMax.push_back(25);
            list_tempMin.push_back(20);
            allData.push_back(list_tempMax);
            allData.push_back(list_tempMin);
            return allData;
        }
        
        std::string line;
        std::getline(myFile, line);
        
        int cont = 0;
        while (cont < simulationTime && std::getline(myFile, line)) {
            if (line.empty()) {
                continue;
            }
            
            std::vector<std::string> lineArray;
            boost::split(lineArray, line, boost::is_any_of(";"));
            
            if (lineArray.size() >= 4) {
                try {
                    std::string trimmedMax = boost::algorithm::trim_copy(lineArray[2]);
                    std::string trimmedMin = boost::algorithm::trim_copy(lineArray[3]);
                    list_tempMax.push_back(std::stoi(trimmedMax));
                    list_tempMin.push_back(std::stoi(trimmedMin));
                } catch (const std::exception& e) {
                    std::cerr << "Error parsing temperature at line " << cont + 2 << ": " << e.what() << std::endl;
                    list_tempMax.push_back(25);
                    list_tempMin.push_back(20);
                }
            } else {
                std::cerr << "Invalid data format at line " << cont + 2 << std::endl;
                list_tempMax.push_back(25);
                list_tempMin.push_back(20);
            }
            cont++;
        }
        
        myFile.close();
        
        if (list_tempMax.empty()) {
            list_tempMax.push_back(25);
            list_tempMin.push_back(20);
        }
        
        allData.push_back(list_tempMax);
        allData.push_back(list_tempMin);
        
    } catch (const std::exception& e) {
        std::cerr << "Error in loadFromExcel: " << e.what() << std::endl;
        if (list_tempMax.empty()) {
            list_tempMax.push_back(25);
            list_tempMin.push_back(20);
            allData.push_back(list_tempMax);
            allData.push_back(list_tempMin);
        }
    }
    
    return allData;
}
/*
std::vector<std::vector<std::pair<double, double>>> ReadData::getCoor(const std::string& filePath) {
    std::vector<std::vector<std::pair<double, double>>> coordinatesList;
    
    try {
        std::ifstream file(filePath);
        if (!file.is_open()) {
            std::cerr << "Cannot open the CSV file: " << filePath << std::endl;
            return coordinatesList;
        }

        std::string line;
        std::getline(file, line);

        while (std::getline(file, line)) {
            if (line.empty()) {
                continue;
            }
            
            std::vector<std::string> fields;
            boost::split(fields, line, boost::is_any_of(","));
            
            if (fields.size() < 13) {  // 4 header fields + 5 coordinate pairs
                std::cerr << "Invalid number of fields in line" << std::endl;
                continue;
            }

            std::vector<std::pair<double, double>> coordinates;
            for (size_t i = 4; i < 13; i += 2) {
                try {
                    std::string trimmedX = boost::algorithm::trim_copy(fields[i]);
                    std::string trimmedY = boost::algorithm::trim_copy(fields[i + 1]);
                    double x = std::stod(trimmedX);
                    double y = std::stod(trimmedY);
                    coordinates.push_back(std::make_pair(x, y));
                } catch (const std::exception& e) {
                    std::cerr << "Error converting coordinates to numbers: " << e.what() << std::endl;
                    continue;
                }
            }

            if (!coordinates.empty()) {
                coordinatesList.push_back(coordinates);
            }
        }

        file.close();
        
        if (coordinatesList.empty()) {
            std::cerr << "Warning: No valid coordinates found in file: " << filePath << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error in getCoor: " << e.what() << std::endl;
    }

    return coordinatesList;
}
*/

std::vector<std::pair<double,double>>
ReadData::getPatchCenters(const std::string& filePath) {
  std::ifstream file(filePath);
  std::string line;
  std::vector<std::pair<double,double>> centers;
  std::getline(file, line); // skip header
  while (std::getline(file, line)) {
    std::vector<std::string> f;
    boost::split(f, line, boost::is_any_of(","));
    std::vector<double> xs, ys;
    // fields 4,5 → first coord, 6,7 → second, …
    for (size_t i = 4; i + 1 < f.size(); i += 2) {
      if (f[i].empty() || f[i+1].empty()) break;
      xs.push_back(std::stod(boost::trim_copy(f[i])));
      ys.push_back(std::stod(boost::trim_copy(f[i+1])));
    }
    if (xs.empty()) continue;
    double cx = std::accumulate(xs.begin(), xs.end(), 0.0) / xs.size();
    double cy = std::accumulate(ys.begin(), ys.end(), 0.0) / ys.size();
    centers.emplace_back(cx, cy);
  }
  return centers;
}
