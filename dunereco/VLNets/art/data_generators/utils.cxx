#include "utils.h"

#include <stdexcept>
#include <string.h>

namespace VLN {

Flavor parseFlavor(const std::string &flavStr)
{
    if (flavStr == "numu") {
        return Flavor::NuMu;
    }

    if (flavStr.empty() || (flavStr == "any")) {
        return Flavor::Any;
    }

    throw std::invalid_argument(
        "Unknown flavor: " + flavStr
        + ". Supported Flavors: 'numu', 'any', ''"
    );
}

Format parseFormat(const std::string &formatStr)
{
    if (formatStr == "csv") {
        return Format::CSV;
    }

    throw std::invalid_argument(
        "Unknown format: " + formatStr + ". Supported Formats: 'csv'"
    );
}

std::string convertFilename(
    const std::string &path, const std::string &root, Format format
)
{
    std::string filename(basename(path.c_str()));
    const size_t dotIdx = filename.find_last_of('.');

    if (dotIdx != std::string::npos) {
        filename.resize(dotIdx);
    }

    switch (format) {
    case Format::CSV:
        return filename + ".csv";
    default:
        throw std::invalid_argument("Unknown format");
    }
}

}

