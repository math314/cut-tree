#pragma once
#include <base/base.h>
#include <graph/graph.h>
#include <vector>
#include <queue>

namespace agl {
namespace cut_tree_internal {
template<class handler_t>
class connected_components_filter {
public:
  connected_components_filter(const G& g)
    : n_(g.num_vertices()), uf_(n_), local_indices_(n_), handlers_indices_(n_), num_connected_components_(0) {

    for (int v = 0; v < n_; v++) for (auto e : g.edges(v)) {
      V u = to(e);
      uf_.unite(u, v);
    }
    num_connected_components_ = 0;
    for (int v = 0; v < n_; v++) {
      if (uf_.root(v) != v) local_indices_[v] = ++local_indices_[uf_.root(v)];
      else handlers_indices_[v] = num_connected_components_++;
    }

    std::vector<bool> used(n_);
    for (int v = 0; v < n_; v++) {
      if (uf_.root(v) != v) continue;
      used[v] = true;

      const int num_vs = local_indices_[uf_.root(v)] + 1;
      local_indices_[uf_.root(v)] = 0;
      std::vector<std::pair<V, V>> edges;
      std::queue<int> q;
      q.push(v);
      while (!q.empty()) {
        V u = q.front(); q.pop();
        for (int dir = 0; dir < 2; dir++) for (auto& e : g.edges(u, D(dir))) {
          V w = to(e);
          if (!used[w]) {
            used[w] = true;
            q.push(w);
          }
          if (dir == 0) {
            edges.emplace_back(local_indices_[u], local_indices_[w]);
          }
        }
      }

      edges.shrink_to_fit();
      handlers_.emplace_back(new handler_t(std::move(edges), num_vs));
    }
  }

  int query(V u, V v) {
    if (!uf_.is_same(u, v)) return 0;
    int lu = local_indices_[u], lv = local_indices_[v];
    CHECK(lu != lv);
    auto& handler = handlers_[handlers_indices_[uf_.root(u)]];
    return handler->query(lu, lv);
  }

  std::vector<std::vector<int>> get_local_id2global_id() {
    const int cc_size = int(handlers().size());
    std::vector<std::vector<int>> local_id2global_id(cc_size);

    const auto& local_indices = this->local_indices();
    for(int v = 0; v < n_; v++) {
      auto& l2g = local_id2global_id[handlers_index(v)];
      if (int(l2g.size()) < local_indices[v] + 1) l2g.resize(local_indices[v] + 1);
      l2g[local_indices[v]] = v;
    }

    return local_id2global_id;
  }

  void print_gomory_hu_tree(std::ostream& os) {
    std::vector<int> roots;
    for(int v = 0; v < n_; v++) if (uf_.root(v) == v) roots.push_back(v);
    for(int i = 0; i < int(roots.size()) - 1; i++) os << roots[0] << " " << roots[i + 1] << " 0\n";

    std::vector<std::vector<int>> local_id2global_id = get_local_id2global_id();

    for(int i = 0; i < int(handlers().size()); i++) {
      const auto& l2g = local_id2global_id[i];
      const auto& gusfield_core = handlers()[i];
      for(int v = 0; v < int(gusfield_core->parent_weight().size()); v++) {
        const auto& kv = gusfield_core->parent_weight()[v];
        int u = kv.first;
        if (u == -1) continue; // 親への辺が存在しない
        int weight = kv.second;
        os << l2g[v] << " " << l2g[u] << " " << weight << "\n";
      }
    }
  }

  int num_connected_components() const { return num_connected_components_; }
  const std::vector<std::unique_ptr<handler_t>>& handlers() const { return handlers_; }
  const std::vector<int>& local_indices() const { return local_indices_; }
  const int handlers_index(int v) { return handlers_indices_[uf_.root(v)]; }

private:
  const int n_;
  union_find uf_;
  std::vector<int> local_indices_, handlers_indices_;
  std::vector<std::unique_ptr<handler_t>> handlers_;
  int num_connected_components_;
};
} // namespace cut_tree_internal
} //namespace agl
