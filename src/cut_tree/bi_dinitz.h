#pragma once
#include <base/base.h>
#include <graph/graph.h>

DECLARE_int32(cut_tree_goal_oriented_dfs_aster_ub);

namespace agl {

class bi_dinitz {
public:
  // two sided bfsが終了した理由
  enum reason_for_finishing_bfs_t {
    kQsIsEmpty,
    kQtIsEmpty,
  };

  class E {
  public:
    E(int to, int rev, int cap) :
      revision_(0), to_(to), cap_(cap), rev_(rev) {
      CHECK(init_cap_ == cap);
    }

    const int cap(int currenct_revision) {
      if (revision_ != currenct_revision) {
        revision_ = currenct_revision;
        cap_ = init_cap_;
      }
      // logging::getcap_counter++;
      return cap_;
    }
    void add_cap(int val, int currenct_revision) {
      if (revision_ != currenct_revision) {
        revision_ = currenct_revision;
        cap_ = init_cap_;
      }
      // logging::addcap_counter++;
      cap_ += val;
    }

    void reset() {
      revision_ = 0;
      cap_ = init_cap_;
    }

  private:
    static const int init_cap_ = 1;
    int revision_;
    int to_;
    int cap_ : 3;
    unsigned int rev_ : 29;

    friend class bi_dinitz;
  };

private:
  bool bi_dfs(int s, int t);
  int dfs(int v, int t, bool use_slevel_, int f);
  void add_undirected_edge(int f, int t, int c);
  void reset_revision();
  int goal_oriented_dfs_inner(int v, int flow, int astar_cost);

  //v -> goal_oriented_bfs_root_ にflowを出来る限り送る
  int goal_oriented_dfs(int v);

public:
  bi_dinitz() : n_(0) {}
  bi_dinitz(const G& g)
    : n_(g.num_vertices()), level_(n_), iter_(n_), bfs_revision_(n_), dfs_revision_(n_), e_(n_),
    s_side_bfs_revision_(2), t_side_bfs_revision_(3), graph_revision_(0), goal_oriented_bfs_root_(-1) {
    for (int v = 0; v < n_; v++) for (auto& e : g.edges(v)) {
      add_undirected_edge(v, agl::to(e), 1);
    }
  }
  bi_dinitz(const std::vector<std::pair<V, V>>& edges, int num_vs)
    : n_(num_vs), level_(n_), iter_(n_), bfs_revision_(n_), dfs_revision_(n_), e_(n_),
    s_side_bfs_revision_(2), t_side_bfs_revision_(3), graph_revision_(0), goal_oriented_bfs_root_(-1) {
    for (auto& uv : edges) {
      add_undirected_edge(uv.first, uv.second, 1);
    }
  }
  bi_dinitz(std::vector<std::pair<V, V>>&& edges, int num_vs)
    : n_(num_vs), level_(n_), iter_(n_), bfs_revision_(n_), dfs_revision_(n_), e_(n_),
    s_side_bfs_revision_(2), t_side_bfs_revision_(3), graph_revision_(0), goal_oriented_bfs_root_(-1) {
    edges.shrink_to_fit();

    //こまめに解放しながら辺を追加していく
    while (edges.size() >= 1) {
      std::size_t loop = std::max(edges.size() / 2, size_t(10000));
      loop = std::min(loop, edges.size());
      while (loop--) {
        auto& uv = edges.back();
        add_undirected_edge(uv.first, uv.second, 1);
        edges.pop_back();
      }
      edges.shrink_to_fit();
    }

    for (auto& e : e_) e.shrink_to_fit();
  }

  void add_vertex();
  void reconnect_edge(E& rm, int sside_vtx, int tside_vtx);

  int max_flow_core(int s, int t);
  int max_flow(int s, int t);

  bool path_dont_exists_to_t(const int v) const;
  bool path_dont_exists_from_s(const int v) const;

  void reset_graph();

  //rootを起点にbfsをして、 s -> rootのflowを高速化する
  void goal_oriented_bfs_init(int root);

  int n() const { return n_; }
  reason_for_finishing_bfs_t reason_for_finishing_bfs() const { return reason_for_finishing_bfs_; }

  std::vector<E>& edges(V v) { return e_[v]; }
  const std::vector<E>& edges(V v) const { return e_[v]; }
  V to(const E& e) const { return e.to_; }
  int cap(E& e) { return e.cap(graph_revision_); }
  E& rev(const E& e_in) { return e_[e_in.to_][e_in.rev_]; }

private:
  int n_;
  std::vector<std::pair<int, int>> level_;
  std::vector<int> iter_;
  std::vector<int> bfs_revision_, dfs_revision_;
  std::vector<std::vector<E>> e_;
  int s_side_bfs_revision_, t_side_bfs_revision_;
  int graph_revision_;
  reason_for_finishing_bfs_t reason_for_finishing_bfs_;

  int goal_oriented_bfs_root_;
  std::vector<int> goal_oriented_bfs_depth_;
};
} //namespace agl