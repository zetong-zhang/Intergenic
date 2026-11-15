#pragma once

#include <string>

[[noreturn]] void logger(const std::string &s);
std::string lower(std::string s);
bool startswith(const std::string &s, const std::string &p);
bool endswith(const std::string &s, const std::string &p);
std::string filename(std::string p);
