#pragma once

#include <string>

[[noreturn]] void logger(const std::string &msg);
std::string lower(std::string msg);
bool startswith(const std::string &msg, const std::string &prefix);
bool endswith(const std::string &msg, const std::string &suffix);
std::string filename(std::string path);
