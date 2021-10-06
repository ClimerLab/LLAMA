#ifndef LLAMA_H
#define LLAMA_H

#include "ConfigParser.h"
#include <string>
#include <vector>

bool checkUserSettings(ConfigParser &);
void readNN_File(const std::string &, std::vector<std::size_t>&);

#endif