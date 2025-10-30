
#ifndef CSVWRITER_H
#define CSVWRITER_H
#include <string>
#include <iostream>
#include<tuple>
#include<vector>
#include <random>

class CSVWriter{
    private:
        std::string fileName;
        std::string delimeter;
    public:
        CSVWriter(std::string Filename = "results.csv", std::string Delimiter = ",");
        void initCSVFile();
        //CSVWriter(std::string  Filename, std::string Delimiter);
        void appendRecord(int* values);
};

#endif