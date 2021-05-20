#pragma once

#include <fstream>
#include <unordered_map>
#include <vector>
#include <string>

#include "dune/VLNets/data/structs/VarDict.h"

class CSVExporter
{
protected:
    std::ofstream ofile;

    std::vector<std::string> scalVarNames;
    std::vector<std::string> vectVarNames;

    void printHeader();
    void init(const VarDict &vars);

    bool initialized;

public:
    explicit CSVExporter(const std::string& output);

    void addScalarVar(const std::string &name);
    void addVectorVar(const std::string &name);

    void setPrecision(int precision);
    void exportVars(const VarDict &vars);
};

