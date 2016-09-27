#pragma once
#include "two_edge_cc_filter.h"
#include "cut_tree_with_2ecc.h"
#include "dinitz.h"
#include "bi_dinitz.h"
#include "cut_tree_query_handler.h"
#include "plain_gomory_hu/gomory_hu_dinitz.h"
#include "plain_gomory_hu/gomory_hu_bi_dinitz.h"

namespace agl {
using cut_tree = agl::cut_tree_internal::two_edge_cc_filter<cut_tree_with_2ecc>; // fastest

using gomory_hu_bi_dinitz = agl::cut_tree_internal::connected_components_filter<agl::cut_tree_internal::plain_gomory_hu::gomory_hu_bi_dinitz>; //faster than plain_gusfield_dinitz
using gomory_hu_dinitz = agl::cut_tree_internal::connected_components_filter<agl::cut_tree_internal::plain_gomory_hu::gomory_hu_dinitz>;
} //namespace agl
