#pragma once
#include <base/base.h>
#include <graph/graph.h>

DECLARE_int32(cut_tree_gtp_dfs_edge_max);

namespace agl {

// ある頂点からdfsをして、貪欲にtree packingを求める
class greedy_treepacking {
  class vecE {
  public:
    vecE() {
      idx_ = 0;
    }

    void add(int to) {
      to_.emplace_back(to);
    }

    bool empty() const { return to_.empty(); }
    int size() const { return int(to_.size()); }
    void advance() {
      idx_++;
      if (int(to_.size()) == idx_) idx_ = 0;
    }
    const V current() const { return to_[idx_]; }
    void remove_current() {
      std::swap(to_[idx_], to_.back());
      to_.pop_back();
      if (int(to_.size()) == idx_) idx_ = 0;
    }
  public:
    int rem_size_, idx_;
    std::vector<V> to_;
  };

  void dfs(int v);
public:
  greedy_treepacking(const std::vector<std::pair<V, V>>& edges, int num_vs);
  void arborescence_packing(int from);

  int inedge_count(int v) const { return inedge_count_[v]; }

private:
  int n_;
  std::vector<vecE> edges_;
  std::vector<int> inedge_count_; // dfs tree packingで、その頂点に何本のin-edgeがあるか
  std::vector<int> used_revision_;
  int vertices_revision_;
};
} // namespace agl