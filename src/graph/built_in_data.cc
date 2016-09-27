#include "built_in_data.h"
#include "built_in_data_inc.h"  // Real data is here
#include <cstring>
using namespace agl;

namespace {
unweighted_graph::edge_list_type array_to_edge_list(V a[][2], size_t n) {
  unweighted_graph::edge_list_type el(n);
  for (auto i : make_irange(n)) {
    el[i].first = a[i][0];
    el[i].second = a[i][1];
  }
  return el;
}
}  // namespace

namespace agl {

#define MATCH_AND_RETURN(id) \
    do { \
      if (strcasecmp(name, #id) == 0) {  \
        return array_to_edge_list(id, sizeof(id) / sizeof(id[0]));  \
      }  \
    } while (false);

unweighted_graph::edge_list_type built_in_edge_list(const char *name) {
  using namespace built_in;

  MATCH_AND_RETURN(karate_club);
  MATCH_AND_RETURN(dolphin);
  MATCH_AND_RETURN(ca_grqc);
  FAIL_MSG("unknown built-in data name specified");
}

#undef MATCH_AND_RETURN

unweighted_graph built_in_graph(const char *name) {
  return unweighted_graph(built_in_edge_list(name));
}
}  // namespace agl
