#include "Utils.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <iostream>

[[noreturn]] void logger(const std::string &s) {
    std::cerr << s << '\n';
    std::exit(1);
}

std::string lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
        [](unsigned char c) { return (char) std::tolower(c); });
    return s;
}

bool startswith(const std::string &s, const std::string &p) {
    return s.compare(0, p.size(), p) == 0;
}

bool endswith(const std::string &s, const std::string &p) {
    return p.size() <= s.size() &&
           s.compare(s.size()-p.size(), p.size(), p) == 0;
}

std::string filename(std::string p) {
    const auto pos = p.find_last_of("/\\");
    return pos == std::string::npos ? p : p.substr(pos+1);
}
