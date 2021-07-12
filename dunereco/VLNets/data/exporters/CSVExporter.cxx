#include "CSVExporter.h"

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <functional>

template<typename T, typename S>
bool printSeparatedValues(
    std::ostream         &output,
    const std::vector<T> &values,
    bool                 firstValuePrinted,
    S                    separator,
    std::function<void (std::ostream &, const T&)> printValue
)
{
    if (values.empty()) {
        return firstValuePrinted;
    }

    size_t startIdx = 0;

    if (! firstValuePrinted)
    {
        printValue(output, values[0]);
        startIdx = 1;
    }

    for (size_t i = startIdx; i < values.size(); i++)
    {
        output << separator;
        printValue(output, values[i]);
    }

    return true;
}

CSVExporter::CSVExporter(const std::string& output)
  : ofile(output), initialized(false)
{
    if (! ofile) {
        throw std::runtime_error("Failed to open output file");
    }

    ofile << std::setprecision(std::numeric_limits<double>::max_digits10);
}

void CSVExporter::setPrecision(int precision)
{
    ofile << std::setprecision(precision);
}

void CSVExporter::printHeader()
{
    bool firstValuePrinted = false;

    firstValuePrinted = printSeparatedValues<std::string, char>(
        ofile, scalVarNames, firstValuePrinted, ',',
        [] (auto &output, const auto &name) { output << name; }
    );

    printSeparatedValues<std::string, char>(
        ofile, vectVarNames, firstValuePrinted, ',',
        [] (auto &output, const auto &name) { output << name; }
    );

    ofile << std::endl;
}

void CSVExporter::init(const VarDict &vars)
{
    if (scalVarNames.empty() && vectVarNames.empty()) {
        for (const auto &kval : vars.scalar) {
            scalVarNames.push_back(kval.first);
        }

        for (const auto &kval : vars.vector) {
            vectVarNames.push_back(kval.first);
        }

        std::sort(scalVarNames.begin(), scalVarNames.end());
        std::sort(vectVarNames.begin(), vectVarNames.end());
    }

    printHeader();
    initialized = true;
}

void CSVExporter::addScalarVar(const std::string &name)
{
    if (! initialized) {
        scalVarNames.push_back(name);
    }
}

void CSVExporter::addVectorVar(const std::string &name)
{
    if (! initialized) {
        vectVarNames.push_back(name);
    }
}

void CSVExporter::exportVars(const VarDict &vars)
{
    if (! initialized) {
        init(vars);
    }

    bool firstValuePrinted = false;

    firstValuePrinted = printSeparatedValues<std::string, char>(
        ofile, scalVarNames, firstValuePrinted, ',',
        [&vars] (auto &output, const auto &name)
        {
            auto it = vars.scalar.find(name);
            if (it != vars.scalar.end()) {
                output << it->second;
            }
        }
    );

    firstValuePrinted = printSeparatedValues<std::string, char>(
        ofile, vectVarNames, firstValuePrinted, ',',
        [&vars] (auto &output, const auto &name)
        {
            output << "\"";
            auto it = vars.vector.find(name);

            if (it != vars.vector.end()) {
                printSeparatedValues<double, char>(
                    output, it->second, false, ',',
                    [] (auto &output, const auto x) { output << x; }
                );
            }

            output << "\"";
        }
    );

    ofile << std::endl;
}


