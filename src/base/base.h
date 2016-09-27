#pragma once

#include <cstddef>
#include <limits>

namespace agl {
static constexpr size_t kInfSize = std::numeric_limits<size_t>::max();
}  // namespace agl

#include "macros.h"
#include "functions.h"
#include "jlog.h"
#include "random.h"
#include "irange.h"
#include "type.h"
#include "data_structures.h"
#include "geometry2d.h"
