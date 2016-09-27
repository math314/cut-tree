#include "bi_dinitz.h"
#include <queue>

DEFINE_int32(cut_tree_goal_oriented_dfs_aster_ub, 2, "bi_dinitz's goal oriented search relaxation");

using namespace std;

namespace agl {

bool bi_dinitz::bi_dfs(int s, int t) {
  queue<int> qs, qt;
  qs.push(s); qt.push(t);
  level_[s].first = level_[t].second = 0;
  bfs_revision_[s] = s_side_bfs_revision_;
  bfs_revision_[t] = t_side_bfs_revision_;

  size_t qs_next_get_cap = e_[s].size();
  size_t qt_next_get_cap = e_[t].size();
  int slevel_ = 0, tlevel_ = 0;
  while (qs.size() != 0 && qt.size() != 0) {
    bool path_found = false;
    if (qs_next_get_cap <= qt_next_get_cap) {
      int size = int(qs.size());
      for (int _ = 0; _ < size; _++) {
        const int v = qs.front(); qs.pop();
        qs_next_get_cap -= e_[v].size();
        for (auto& t : e_[v]) {
          if (t.cap(graph_revision_) == 0 || bfs_revision_[t.to_] == s_side_bfs_revision_) continue;
          if (bfs_revision_[t.to_] == t_side_bfs_revision_) {
            path_found = true;
            continue;
          }
          bfs_revision_[t.to_] = s_side_bfs_revision_;
          level_[t.to_].first = slevel_ + 1;
          qs_next_get_cap += e_[t.to_].size();
          qs.push(t.to_);
        }
      }
      slevel_++;
    } else {
      int size = int(qt.size());
      for (int _ = 0; _ < size; _++) {
        const int v = qt.front(); qt.pop();
        qt_next_get_cap -= e_[v].size();
        for (auto& t : e_[v]) {
          if (e_[t.to_][t.rev_].cap(graph_revision_) == 0 || bfs_revision_[t.to_] == t_side_bfs_revision_) continue;
          if (bfs_revision_[t.to_] == s_side_bfs_revision_) {
            path_found = true;
            continue;
          }
          bfs_revision_[t.to_] = t_side_bfs_revision_;
          level_[t.to_].second = tlevel_ + 1;
          qt_next_get_cap += e_[t.to_].size();
          qt.push(t.to_);
        }
      }
      tlevel_++;
    }
    // fprintf(stderr, "slevel_ : %d, tlevel_ : %d\n",slevel_, tlevel_);
    if (path_found) return true;
  }

  reason_for_finishing_bfs_ = (qs.empty()) ? kQsIsEmpty : kQtIsEmpty;
  return false;
}

int bi_dinitz::dfs(int v, int t, bool use_slevel, int f) {
  // cut_tree_goal_oriented_dfs_aster_ub >= 3 を設定すると、dfs中に同じ辺を使ってしまい、辺のコストが破綻して f < 0 となることがある
  CHECK(f >= 0);

  if (v == t) return f;
  if (dfs_revision_[v] != bfs_revision_[v]) {
    dfs_revision_[v] = bfs_revision_[v];
    iter_[v] = 0;
  }
  for (int &i = iter_[v]; i < int(e_[v].size()); i++) {
    E& _e = e_[v][i];
    const int cap = _e.cap(graph_revision_);
    if (cap == 0 || bfs_revision_[_e.to_] / 2 != s_side_bfs_revision_ / 2) continue;

    bool rec;
    if (use_slevel) rec = bfs_revision_[_e.to_] == t_side_bfs_revision_ || level_[v].first < level_[_e.to_].first;
    else rec = bfs_revision_[_e.to_] == t_side_bfs_revision_ && level_[v].second > level_[_e.to_].second;
    if (!rec) continue;

    bool next_slevel_ = use_slevel && bfs_revision_[_e.to_] == s_side_bfs_revision_;
    int d = dfs(_e.to_, t, next_slevel_, min(f, cap));
    if (d > 0) {
      _e.add_cap(-d, graph_revision_);
      e_[_e.to_][_e.rev_].add_cap(d, graph_revision_);
      return d;
    }
  }
  return 0;
}

void bi_dinitz::add_undirected_edge(int f, int t, int c) {
  e_[f].push_back(E(t, int(e_[t].size()), c));
  e_[t].push_back(E(f, int(e_[f].size()) - 1, c));
}

void bi_dinitz::reset_revision() {
  for (int v = 0; v < n_; v++) for (auto& e : e_[v]) e.reset();
  memset(bfs_revision_.data(), 0, sizeof(bfs_revision_[0]) * bfs_revision_.size());
  memset(dfs_revision_.data(), 0, sizeof(dfs_revision_[0]) * dfs_revision_.size());
  s_side_bfs_revision_ = 2;
  t_side_bfs_revision_ = 3;
  graph_revision_ = 0;
}

int bi_dinitz::goal_oriented_dfs_inner(int v, int flow, int astar_cost) {
  if (v == goal_oriented_bfs_root_)
    return flow;

  if (dfs_revision_[v] != s_side_bfs_revision_) {
    dfs_revision_[v] = s_side_bfs_revision_;
    iter_[v] = 0;
  }
  for (int &i = iter_[v]; i < int(e_[v].size()); i++) {
    E& to_edge = e_[v][i];
    int to = to_edge.to_;
    int add_aster_cost = goal_oriented_bfs_depth_[to] - goal_oriented_bfs_depth_[v] + 1;
    if (add_aster_cost == 2) return 0; //コストの増える頂点は辿らない
    int n_astar_cost = astar_cost + add_aster_cost;
    if (n_astar_cost > FLAGS_cut_tree_goal_oriented_dfs_aster_ub) {
      //sort済なので, これより後で自身の深さよりも浅い頂点は存在しない
      return 0;
    }
    int cap = to_edge.cap(graph_revision_);
    if (cap == 0) continue;
    int d = goal_oriented_dfs_inner(to, min(flow, cap), n_astar_cost);
    if (d > 0) {
      to_edge.add_cap(-d, graph_revision_);
      e_[to_edge.to_][to_edge.rev_].add_cap(d, graph_revision_);
      return d;
    }
  }

  return 0;
}

//v -> goal_oriented_bfs_root_ にflowを出来る限り送る
int bi_dinitz::goal_oriented_dfs(int v) {
  int flow = 0;
  s_side_bfs_revision_ += 2;
  t_side_bfs_revision_ += 2;

  for (auto it = e_[v].begin(); it != e_[v].end(); ++it) {
    auto& to_edge = *it;
    while (to_edge.cap(graph_revision_) > 0) {
      int add = goal_oriented_dfs_inner(to_edge.to_, to_edge.cap(graph_revision_), 0);
      if (add == 0) break;
      flow += add;
      to_edge.add_cap(-add, graph_revision_);
      e_[to_edge.to_][to_edge.rev_].add_cap(add, graph_revision_);
    }
  }
  return flow;
}

void bi_dinitz::add_vertex() {
  level_.emplace_back();
  iter_.emplace_back();
  bfs_revision_.emplace_back();
  dfs_revision_.emplace_back();
  e_.emplace_back();
  n_++;
}

void bi_dinitz::reconnect_edge(E& rm, int sside_vtx, int tside_vtx) {
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

int bi_dinitz::max_flow_core(int s, int t) {
  assert(s != t);
  reset_graph();

  int flow = 0;
  int preflow = 0;
  if (goal_oriented_bfs_root_ == t) {
    preflow = goal_oriented_dfs(s);
  }

  if (preflow == int(e_[s].size())) {
    reason_for_finishing_bfs_ = kQsIsEmpty;
  } else {
    s_side_bfs_revision_ += 2;
    t_side_bfs_revision_ += 2;

    int bfs_counter = 0;
    for (;; s_side_bfs_revision_ += 2, t_side_bfs_revision_ += 2) {
      // fprintf(stderr, "bfs_start\n");
      bool path_found = bi_dfs(s, t);
      bfs_counter++;
      if (!path_found) break;
      while (true) {
        int f = dfs(s, t, true, numeric_limits<int>::max());
        if (f == 0) break;
        flow += f;
      }
    }
    // fprintf(stderr, "bfs_counter : %d\n", bfs_counter);
  }
  // fprintf(stderr, "(%d,%d) : preflow = %d, flow = %d\n", s, t, preflow, flow);
  if (flow == 0 && preflow > 0) {
    if (int(e_[s].size()) == preflow) {
      // logging::preflow_eq_degree++;
    }
    // logging::flow_eq_0++;
  } else {
    // fprintf(stderr, "// logging::getcap_counter : %lld\n", // logging::getcap_counter);
    // fprintf(stderr, "// logging::addcap_counter : %lld\n", // logging::addcap_counter);
  }
  return flow + preflow;
}

int bi_dinitz::max_flow(int s, int t) {
  int ans = max_flow_core(s, t);
  return ans;
}

bool bi_dinitz::path_dont_exists_to_t(const int v) const {
  if (reason_for_finishing_bfs_ == kQsIsEmpty) {
    //sから到達可能な頂点のbfs_revision_には、必ずs_side_bfs_revision_が代入されている
    return bfs_revision_[v] == s_side_bfs_revision_;
  } else {
    //tから到達不可能
    return bfs_revision_[v] != t_side_bfs_revision_;
  }
}

bool bi_dinitz::path_dont_exists_from_s(const int v) const {
  if (reason_for_finishing_bfs_ == kQsIsEmpty) {
    //sから到達不可能
    return bfs_revision_[v] != s_side_bfs_revision_;
  } else {
    //tから到達可能な頂点のbfs_revision_には、必ずt_side_bfs_revision_が代入されている
    return bfs_revision_[v] == t_side_bfs_revision_;
  }
}

//フローを流す前に実行する
void bi_dinitz::reset_graph() {
  graph_revision_++;
  if (s_side_bfs_revision_ >= numeric_limits<decltype(s_side_bfs_revision_)>::max() / 2) {
    reset_revision();
  }
}

//rootを起点にbfsをして、 s -> rootのflowを高速化する
void bi_dinitz::goal_oriented_bfs_init(int root) {
  goal_oriented_bfs_root_ = root;
  goal_oriented_bfs_depth_.clear();
  goal_oriented_bfs_depth_.resize(n_, n_); // bfsの深さをn(=INF)で初期化

  // set goal_oriented_bfs_depth_
  queue<int> q;
  q.push(root);
  goal_oriented_bfs_depth_[root] = 0;
  while (!q.empty()) {
    const int v = q.front(); q.pop();
    const int ndepth = goal_oriented_bfs_depth_[v] + 1;
    for (auto& to_edge : e_[v]) {
      if (goal_oriented_bfs_depth_[to_edge.to_] > ndepth) {
        goal_oriented_bfs_depth_[to_edge.to_] = ndepth;
        q.push(to_edge.to_);
      }
    }
  }

  //sort edges by depth order
  for (int v = 0; v < n_; v++) {
    auto& dep = goal_oriented_bfs_depth_;
    sort(e_[v].begin(), e_[v].end(), [&dep](const E& l, const E& r) {
      return dep[l.to_] < dep[r.to_];
    });

    //reset rev_ edge's "rev_" value
    for (int i = 0; i < int(e_[v].size()); i++) {
      const E& to_edge = e_[v][i];
      E& rev = e_[to_edge.to_][to_edge.rev_];
      rev.rev_ = i;
    }
  }
}

} //namespace agl