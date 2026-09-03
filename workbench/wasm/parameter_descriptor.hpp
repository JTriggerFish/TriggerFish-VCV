#pragma once

#include <string>

namespace tfworkbench {

enum class ParameterScale : int { Linear, Logarithmic, Boolean };

struct ParameterDescriptor {
  std::string key;
  std::string name;
  std::string unit;
  float minimum;
  float maximum;
  float defaultValue;
  ParameterScale scale{ParameterScale::Linear};
};

} // namespace tfworkbench
