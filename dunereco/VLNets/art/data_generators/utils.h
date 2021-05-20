#pragma once
#include <string>

namespace VLN {

enum class Flavor : int { Any = 0, NuMu = 14 };
enum class Format { CSV };

Flavor parseFlavor(const std::string &flavStr);
Format parseFormat(const std::string &formatStr);

std::string convertFilename(
    const std::string &path, const std::string &root, Format format
);

}
