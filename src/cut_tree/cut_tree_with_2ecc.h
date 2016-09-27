#pragma once
#include <base/base.h>
#include <graph/graph.h>

DECLARE_int32(cut_tree_try_greedy_tree_packing);
DECLARE_int32(cut_tree_try_large_degreepairs);
DECLARE_int32(cut_tree_separate_near_pairs_d);
DECLARE_int32(cut_tree_contraction_lower_bound);
DECLARE_bool(cut_tree_enable_greedy_tree_packing);
DECLARE_bool(cut_tree_enable_adjacent_cut);
DECLARE_bool(cut_tree_enable_goal_oriented_search);

namespace agl {
namespace cut_tree_internal {
class disjoint_cut_set;
class separator;
class gomory_hu_tree_builder;
} // cut_tree_internal

// 2ecc = two-edge connected components
class cut_tree_with_2ecc {
  void find_cuts_by_tree_packing(std::vector<std::pair<V, V>>& edges, cut_tree_internal::disjoint_cut_set* dcs, const std::vector<int>& degree);
  void contract_degree2_vertices(std::vector<std::pair<V, V>>& edges, std::vector<int>& degree);

  //次数の大きい頂点対をcutする
  void separate_high_degreepairs(cut_tree_internal::separator* sep);

  //隣接頂点同士を見て、まだ切れていなかったらcutする
  void separate_adjacent_pairs(cut_tree_internal::separator* sep);

  void separate_all(cut_tree_internal::separator* sep);

  void separate_near_pairs(cut_tree_internal::separator* sep);

  //次数の最も高い頂点に対して、出来る限りの頂点からflowを流してmincutを求める
  void find_cuts_by_goal_oriented_search(cut_tree_internal::separator* sep);

public:

  cut_tree_with_2ecc(std::vector<std::pair<V, V>>&& edges, int num_vs);
  ~cut_tree_with_2ecc();

  int query(V u, V v) const;
  const std::vector<std::pair<V, int>>& parent_weight() const;

private:
  const int num_vertices_;
  std::unique_ptr<cut_tree_internal::gomory_hu_tree_builder> gh_builder_;
};

} // namespace agl