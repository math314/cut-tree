#pragma once
#include "edge_type.h"
#include "basic_graph.h"

namespace agl {
template<typename WeightType = double>
using weighted_graph = basic_graph<weighted_edge<WeightType>>;
}  // namespace agl
