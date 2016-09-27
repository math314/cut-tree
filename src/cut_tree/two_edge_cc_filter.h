#pragma once
#include "connected_components_filter.h"

namespace agl {
namespace cut_tree_internal {
// two-edge connected components filter
template<class handler_t>
class two_edge_cc_filter {
public:

  void lowlink_dfs(int v, int par, int& cur_ord) {
    lowlink_[v] = order_[v] = cur_ord++;
    for (int dir = 0; dir < 2; dir++) for (auto to : g_.edges(v, D(dir))) {
      if (to == par) continue;
      if (order_[to] == -1) {
        lowlink_dfs(to, v, cur_ord);
        lowlink_[v] = std::min(lowlink_[v], lowlink_[to]);
        if (order_[v] < lowlink_[to]) bridge_.emplace_back(v, to);
        else biconnected_graphs_edges_.emplace_back(v, to);
      } else {
        lowlink_[v] = std::min(lowlink_[v], lowlink_[to]);
        if (v < to) biconnected_graphs_edges_.emplace_back(v, to);
      }
    }
  }

  std::vector<std::vector<int>> get_local_id2global_id() const {
    const int cc_size = int(biconnected_graph_handler_->handlers().size());
    std::vector<std::vector<int>> local_id2global_id(cc_size);

    const auto& local_indices = biconnected_graph_handler_->local_indices();
    for (int v = 0; v < n_; v++) {
      auto& l2g = local_id2global_id[biconnected_graph_handler_->handlers_index(v)];
      if (int(l2g.size()) < local_indices[v] + 1) l2g.resize(local_indices[v] + 1);
      l2g[local_indices[v]] = v;
    }

    return local_id2global_id;
  }

  two_edge_cc_filter(G& g) : n_(g.num_vertices()), g_(g), uf_(n_), lowlink_(n_, -1), order_(n_, -1) {

    G new_g;
    for (int v = 0; v < n_; v++) for (auto& e : g_.edges(v)) {
      uf_.unite(v, to(e));
    }

    const int num_edges = g_.num_edges();

    for (int v = 0; v < n_; v++) if (uf_.root(v) == v) {
      int cur_ord = 0;
      lowlink_dfs(v, -1, cur_ord);
    }

    CHECK(bridge_.size() + biconnected_graphs_edges_.size() == std::size_t(num_edges));

    g_.clear_and_shrink_to_fit();
    new_g = G(biconnected_graphs_edges_, n_);

    //dealloc
    lowlink_.clear(); lowlink_.shrink_to_fit();
    order_.clear(); order_.shrink_to_fit();
    // bridge_.clear(); bridge_.shrink_to_fit();
    biconnected_graphs_edges_.clear(); biconnected_graphs_edges_.shrink_to_fit();

    biconnected_graph_handler_.reset(new connected_components_filter<handler_t>(new_g));
  }

public:
  int query(V u, V v) {
    int ans = biconnected_graph_handler_->query(u, v);
    if (ans == 0) {
      if (uf_.is_same(u, v)) return 1; // 橋で間接的につながっている
      else return 0;
    }
    return ans;
  }

  void print_gomory_hu_tree(std::ostream& os) {
    std::vector<int> roots;
    for (int v = 0; v < n_; v++) if (uf_.root(v) == v) roots.push_back(v);
    for (int i = 0; i < int(roots.size()) - 1; i++) os << roots[0] << " " << roots[i + 1] << " 0\n";
    for (auto& e : bridge_) os << e.first << " " << e.second << " 1\n";

    std::vector<std::vector<int>> local_id2global_id = get_local_id2global_id();

    //weight2以上
    for (int i = 0; i < int(biconnected_graph_handler_->handlers().size()); i++) {
      const auto& l2g = local_id2global_id[i];
      const auto& gusfield_core = biconnected_graph_handler_->handlers()[i];
      for (int v = 0; v < int(gusfield_core->parent_weight().size()); v++) {
        const auto& kv = gusfield_core->parent_weight()[v];
        int u = kv.first;
        if (u == -1) continue; // 親への辺が存在しない
        int weight = kv.second;
        CHECK(weight >= 2);
        os << l2g[v] << " " << l2g[u] << " " << weight << "\n";
      }
    }
  }

private:
  const int n_;
  G& g_;
  union_find uf_;
  std::vector<int> lowlink_, order_;
  std::vector<std::pair<V, V>> bridge_, biconnected_graphs_edges_;

  std::unique_ptr<connected_components_filter<handler_t>> biconnected_graph_handler_;
};
} // namespace cut_tree_internal
} // namespace agl
