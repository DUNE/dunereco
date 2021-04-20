#pragma once

#include <unordered_map>
#include <vector>
#include <string>

// TODO: Port from double to std::variant when new c++ standard is out
struct VarDict
{
    std::unordered_map<std::string, double>               scalar;
    std::unordered_map<std::string, std::vector<double>>  vector;
};

