#pragma once
#include <base/base.h>
#include <graph/graph.h>

namespace agl {
namespace cut_tree_internal {
class dinitz {
  struct E {
    int to_, rev_, cap_;
    E(int to, int rev_, int cap_) : to_(to), rev_(rev_), cap_(cap_) {}
  };

  void bfs(int s);
  int dfs(int v, int t, int f);
  void reset_graph();

public:
  dinitz(const G& g);
  int max_flow(int s, int t);

  std::vector<E>& edges(V v) { return e_[v]; }
  V to(const E& e) const { return e.to_; }
  int cap(E& e) { return e.cap_; }
  E& rev(const E& e_in) { return e_[e_in.to_][e_in.rev_]; }
  void add_undirected_edge(int f, int t, int c);

  void reconnect_edge(E& rm, int sside_vtx, int tside_vtx) {
    const int to = rm.to_;
    const int to_rev = rm.rev_;
    E& rm_rev = e_[to][to_rev];
    const int from = rm_rev.to_;
    const int from_rev = rm_rev.rev_;

    add_undirected_edge(sside_vtx, tside_vtx, 1);
    E& se = e_[sside_vtx].back();
    E& te = e_[tside_vtx].back();

    rm.to_ = sside_vtx;
    rm.rev_ = te.rev_;
    rm_rev.to_ = tside_vtx;
    rm_rev.rev_ = se.rev_;

    se.to_ = from;
    se.rev_ = from_rev;
    te.to_ = to;
    te.rev_ = to_rev;
  }

  void add_vertex() {
    level_.emplace_back();
    iter_.emplace_back();
    e_.emplace_back();
    n_++;
  }

  int n() const { return n_; }

private:
  G g;
  std::vector<int> level_, iter_;
  std::vector<std::vector<E>> e_;
  V n_;
};

} //namespace cut_tree_internal
} //namespace agl
