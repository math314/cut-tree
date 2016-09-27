#pragma once
#include "separator.h"

namespace agl {
namespace cut_tree_internal {
namespace plain_gomory_hu {
class gomory_hu_bi_dinitz {
  void separate_all(separator& sep) {
    const disjoint_cut_set& dcs = sep.get_disjoint_cut_set();
    for(int group_id = 0; group_id < num_vertices_; group_id++) {
      while (dcs.has_two_elements(group_id)) {
        V s, t; std::tie(s, t) = dcs.get_two_elements(group_id);
        sep.mincut(s, t);
      }
    }
  }

public:

  gomory_hu_bi_dinitz(std::vector<std::pair<V, V>>&& edges, int num_vs) :
    num_vertices_(num_vs),
    gh_builder_(num_vs) {
    //G g(edges);

    if(num_vs > 10000) fprintf(stderr, "gomory_hu_bi_dinitz::constructor start : memory %ld MB\n", jlog_internal::get_memory_usage() / 1024);

    disjoint_cut_set dcs(num_vs);

    //dinicの初期化
    bi_dinitz dz_base(std::move(edges), num_vs);
    if(num_vs > 10000) fprintf(stderr, "gomory_hu_bi_dinitz::bi_dinitz after init : memory %ld MB\n", jlog_internal::get_memory_usage() / 1024);

    separator sep(dz_base, dcs, gh_builder_);

    // 残った頂点groupをcutする、gomory_hu treeの完成
    JLOG_ADD_BENCHMARK_IF("time.separate_all", num_vertices_ > 10000) {
      separate_all(sep);
    }

    sep.output_debug_infomation();

    gh_builder_.build();
    //gh_builder_.test(g);
  }

  int query(V u, V v) const {
    return gh_builder_.query(u, v);
  }

  const std::vector<std::pair<V, int>>& parent_weight() const {
    return gh_builder_.parent_weight();
  }

private:
  const int num_vertices_;
  //todo .ccと.hに分ける時に、std::unique_ptrで囲む
  gomory_hu_tree_builder gh_builder_;
};
} // namespace plain_gomory_hu
} // namespace cut_tree_internal
} // namespace agl
