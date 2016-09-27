#include "greedy_treepacking.h"

DEFINE_int32(cut_tree_gtp_dfs_edge_max, 1000000000, "greedy tree packing breadth limit");
using namespace std;

namespace agl {

void greedy_treepacking::dfs(int v) {
  used_revision_[v] = vertices_revision_;
  auto& to_edges = edges_[v];
  const int rem_edges = min(to_edges.size(), FLAGS_cut_tree_gtp_dfs_edge_max);
  for (int i = 0; i < rem_edges; i++) {
    V to = to_edges.current();
    // logging::gtp_edge_count++;
    if (used_revision_[to] == vertices_revision_) {
      to_edges.advance();
      // logging::gtp_edge_miss++;
    } else {
      to_edges.remove_current();
      inedge_count_[to]++; // in-edgeの本数が増える
      dfs(to);
      // logging::gtp_edge_use++;
    }
  }
}
greedy_treepacking::greedy_treepacking(const vector<pair<V, V>>& edges, int num_vs) :
  n_(num_vs), edges_(n_), inedge_count_(n_), used_revision_(n_), vertices_revision_(1) {
  for (auto& e : edges) {
    edges_[e.first].add(e.second);
    edges_[e.second].add(e.first);
  }

  for (auto& e : edges_) {
    auto& edges_ = this->edges_;
    sort(e.to_.begin(), e.to_.end(), [&edges_](const V l, const V r) {
      int a = edges_[l].size();
      int b = edges_[r].size();
      return a < b;
    });
  }
}

void greedy_treepacking::arborescence_packing(int from) {
  for (int i = 0; i < edges_[from].size(); i++) {
    dfs(from);
    vertices_revision_++;
  }
}

} // namespace agl