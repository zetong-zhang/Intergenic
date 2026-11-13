#include "Utils.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <iostream>

[[noreturn]] void logger(const std::string &msg) {
    std::cerr << msg << '\n';
    std::exit(1);
}

std::string lower(std::string msg) {
    std::transform(msg.begin(), msg.end(), msg.begin(),
        [](unsigned char c) { return (char) std::tolower(c); });
    return msg;
}

bool startswith(const std::string &msg, const std::string &prefix) {
    return msg.compare(0, prefix.size(), prefix) == 0;
}

bool endswith(const std::string &msg, const std::string &suffix) {
    return suffix.size() <= msg.size() &&
           msg.compare(msg.size()-suffix.size(), suffix.size(), suffix) == 0;
}

std::string filename(std::string path) {
    const auto pos = path.find_last_of("/\\");
    return pos == std::string::npos ? path : path.substr(pos+1);
}
