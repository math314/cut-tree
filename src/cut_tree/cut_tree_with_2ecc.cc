#include "cut_tree_with_2ecc.h"
#include "bi_dinitz.h"
#include "greedy_treepacking.h"
#include <queue>
#include <unordered_set>

DEFINE_int32(cut_tree_try_greedy_tree_packing, 1, "number of tree packing");
DEFINE_int32(cut_tree_try_large_degreepairs, 10, "number of separate large degree pairs");
DEFINE_int32(cut_tree_separate_near_pairs_d, 1, "separate near pairs radius");
DEFINE_int32(cut_tree_contraction_lower_bound, 2, "contraction upper bound");
DEFINE_bool(cut_tree_enable_greedy_tree_packing, true, "");
DEFINE_bool(cut_tree_enable_adjacent_cut, true, "");
DEFINE_bool(cut_tree_enable_goal_oriented_search, true, "");

using namespace std;
using namespace agl::cut_tree_internal;

namespace agl {
namespace cut_tree_internal {
class disjoint_cut_set {
  struct Node {
    int pv, nt;
    int root;
  };

  void erase(int node_id) {
    group_size_[nodes[node_id].root]--;

    int pv = nodes[node_id].pv, nt = nodes[node_id].nt;
    if (pv != -1) {
      nodes[pv].nt = nt;
    }
    if (nt != -1) {
      nodes[nt].pv = pv;
    }
    if (pv == -1) {
      root[nodes[node_id].root] = nt;
    }
  }

  void add(int node_id, int group_id) {
    group_size_[group_id]++;

    int nt = root[group_id];
    nodes[node_id].root = group_id;
    nodes[node_id].pv = -1;
    nodes[node_id].nt = nt;
    root[group_id] = node_id;
    if (nt != -1) {
      nodes[nt].pv = node_id;
    }
  }

public:
  disjoint_cut_set(int n) : root(n, -1), nodes(n), group_num(1), group_size_(n) {
    root[0] = 0;
    nodes[0].pv = -1;
    nodes[n - 1].nt = -1;
    for (int i = 0; i < n - 1; i++) {
      nodes[i].nt = i + 1;
      nodes[i + 1].pv = i;
    }
    for (int i = 0; i < n; i++) nodes[i].root = 0;
    group_size_[0] = n;
  }

  const int node_num() const {
    return int(nodes.size());
  }

  void create_new_group(int id) {
    erase(id);
    add(id, group_num++);
  }

  bool is_same_group(int a, int b) const {
    if (a >= int(nodes.size()) || b >= int(nodes.size())) return false;
    return nodes[a].root == nodes[b].root;
  }

  void move_other_group(int src, int dst) {
    erase(src);
    add(src, nodes[dst].root);
  }

  int other_id_in_same_group(int id) const {
    const int grp_id = nodes[id].root;
    const int rt = root[grp_id];
    CHECK(rt != -1);
    if (rt != id) return rt;
    const int nxt = nodes[rt].nt;
    CHECK(nxt != -1);
    return nxt;
  }

  pair<int, int> get_two_elements(int group_id) const {
    const int rt = root[group_id];
    CHECK(rt != -1);
    const int nxt = nodes[rt].nt;
    if (nxt == -1) return make_pair(-1, -1);
    return make_pair(rt, nxt);
  }

  bool has_two_elements(int group_id) const {
    auto uv = get_two_elements(group_id);
    return uv.first != -1;
  }

  vector<int> get_group(int group_id) const {
    vector<int> ret;
    int cur = root[group_id];
    while (cur != -1) {
      ret.push_back(cur);
      cur = nodes[cur].nt;
    }
    return ret;
  }

  int group_id(int id) const {
    return nodes[id].root;
  }

  int group_size(int grp_id) const {
    return group_size_[grp_id];
  }

  int debug_group_num() const {
    return group_num;
  }

private:
  vector<int> root;
  vector<Node> nodes;
  int group_num;
  vector<int> group_size_;
};

class gomory_hu_tree_builder {
  void dfs(V v, V par = -1) {
    int dep = (par == -1) ? 0 : depth_[par] + 1;
    depth_[v] = dep;
    for (auto& to : edges_[v]) {
      if (get<0>(to) == par) continue;
      parent_cost_[get<0>(to)] = make_pair(v, get<1>(to));
      dfs(get<0>(to), v);
    }
  }

public:
  gomory_hu_tree_builder(int n) : n_(n), edges_(n), depth_(n), parent_cost_(n) {
    add_edge_count_ = 0;
  }

  void add_degree2_edge(V u, V v) {
    add_edge_count_++;
    degree2_edges_.emplace_back(u, v);
  }

  void contraction(V s, V t, V sside_new_vtx, V tside_new_vtx) {
    edges_.resize(edges_.size() + 2);
    CHECK(get<0>(edges_[s].back()) == t);
    CHECK(get<0>(edges_[t].back()) == s);
    int f = get<1>(edges_[s].back());

    get<0>(edges_[s].back()) = sside_new_vtx;
    get<2>(edges_[s].back()) = 0;
    edges_[sside_new_vtx].emplace_back(s, f, edges_[s].size() - 1);

    get<0>(edges_[t].back()) = tside_new_vtx;
    get<2>(edges_[t].back()) = 0;
    edges_[tside_new_vtx].emplace_back(t, f, edges_[t].size() - 1);

    edges_[sside_new_vtx].emplace_back(tside_new_vtx, f, 1);
    edges_[tside_new_vtx].emplace_back(sside_new_vtx, f, 1);
  }

  void add_edge(V u, V v, int cost, const vector<V>& vs, const disjoint_cut_set* dcs) {
    CHECK(u != v);
    add_edge_count_++;
    for (V w : vs) {
      if (!dcs->is_same_group(u, w)) continue;
      for (size_t i = 0; i < edges_[w].size(); ++i) {
        V t = get<0>(edges_[w][i]);
        int r = get<2>(edges_[w][i]);
        edges_[v].emplace_back(edges_[w][i]);
        get<0>(edges_[t][r]) = v;
        get<2>(edges_[t][r]) = edges_[v].size() - 1;
      }
      edges_[w].clear();
    }
    for (V w : vs) {
      if (w == u) continue;
      for (size_t i = 0; i < edges_[w].size(); ++i) {
        V t = get<0>(edges_[w][i]);
        if (!dcs->is_same_group(t, v)) continue;
        int r = get<2>(edges_[w][i]);
        edges_[t][r] = edges_[t].back();
        get<2>(edges_[get<0>(edges_[t][r])][get<2>(edges_[t][r])]) = r;
        edges_[t].pop_back();
        get<0>(edges_[w][i]) = u;
        get<2>(edges_[w][i]) = edges_[u].size();
        edges_[u].emplace_back(w, get<1>(edges_[w][i]), i);
      }
    }
    edges_[u].emplace_back(v, cost, edges_[v].size());
    edges_[v].emplace_back(u, cost, edges_[u].size() - 1);
  }

  void build() {
    for (auto e : degree2_edges_) {
      edges_[e.first].emplace_back(e.second, 2, edges_[e.second].size());
      edges_[e.second].emplace_back(e.first, 2, edges_[e.first].size() - 1);
    }
    degree2_edges_.clear();
    for (V v = n_; v < (V)(edges_.size()); ++v) {
      CHECK(edges_[v].size() == 2);
      for (auto& e : edges_[get<0>(edges_[v][0])]) {
        if (get<0>(e) == v) {
          get<0>(e) = get<0>(edges_[v][1]);
        }
      }
      for (auto& e : edges_[get<0>(edges_[v][1])]) {
        if (get<0>(e) == v) {
          get<0>(e) = get<0>(edges_[v][0]);
        }
      }
    }
    CHECK(add_edge_count_ == n_ - 1);
    parent_cost_[0] = make_pair(-1, 0);
    dfs(0);
    edges_.clear(); edges_.shrink_to_fit();
  }

  int query(V u, V v) const {
    CHECK(u != v);
    CHECK(u < n_ && v < n_);
    int ans = numeric_limits<int>::max();
    while (u != v) {
      if (depth_[u] > depth_[v]) {
        ans = min(ans, parent_cost_[u].second);
        u = parent_cost_[u].first;
      } else {
        ans = min(ans, parent_cost_[v].second);
        v = parent_cost_[v].first;
      }
    }
    return ans;
  }

  const vector<pair<V, int>>& parent_weight() const {
    return parent_cost_;
  }

  int debug_add_edge_count() const {
    return add_edge_count_;
  }

  void test(const G& g) {
    if (n_ == 1) {
      return;
    }
    queue<V> que;
    vector<vector<V>> edge(n_), children(n_);
    vector<V> cnt(n_);
    vector<unordered_set<V>> d(n_);
    for (V v : irange<V>(n_)) {
      for (auto e : g.neighbors(v)) {
        edge[v].emplace_back(to(e));
        edge[to(e)].emplace_back(v);
      }
    }
    for (V v : irange<V>(n_)) {
      if (parent_cost_[v].first != -1) {
        children[parent_cost_[v].first].emplace_back(v);
        ++cnt[parent_cost_[v].first];
      }
    }
    for (V v : irange<V>(n_)) {
      if (cnt[v] == 0) {
        que.emplace(v);
      }
    }
    while (!que.empty()) {
      V v = que.front();
      que.pop();
      for (V u : children[v]) {
        for (V t : d[u]) {
          d[v].emplace(t);
        }
      }
      d[v].emplace(v);
      V cut = 0;
      for (V u : d[v]) {
        for (V t : edge[u]) {
          if (d[v].count(t) == 0) {
            ++cut;
          }
        }
      }
      CHECK(cut == parent_cost_[v].second);
      if (parent_cost_[v].first >= 0) {
        --cnt[parent_cost_[v].first];
        if (cnt[parent_cost_[v].first] == 0) {
          que.emplace(parent_cost_[v].first);
        }
      }
    }
  }

private:
  int add_edge_count_;
  int n_;
  vector<vector<tuple<V, int, int>>> edges_;
  vector<int> depth_;
  vector<pair<V, int>> parent_cost_;
  vector<pair<V, V>> degree2_edges_;
};

class separator {

  const int used_flag_value() const {
    return max_flow_times_;
  }

  // 一定期間置きに進捗を出力する
  void print_progress_at_regular_intervals(V s, V t, int cost) {
    if (max_flow_times_ % 10000 == 0) {
      stringstream ss;
      ss << "max_flow_times_ = " << max_flow_times_ << ", (" << s << "," << t << ") cost = " << cost;
      JLOG_ADD("separator.progress", ss.str());
      fprintf(stderr, "cut details : ");
      for (auto& kv : debug_count_cut_size_for_a_period_) fprintf(stderr, "(%d,%d), ", kv.first, kv.second);
      fprintf(stderr, "\n");
      debug_count_cut_size_for_a_period_.clear();
    }
  }

  int max_flow(const V s, const V t) {
    int cost = dz_.max_flow(s, t);
    debug_last_max_flow_cost_ = cost;

    // fprintf(stderr, "(%d,%d) : %d\n", s, t, cost);
    //debug infomation

    max_flow_times_++;
    print_progress_at_regular_intervals(s, t, cost);

    cross_other_mincut_count_ = 0;
    auto check_crossed_mincut = [this](const V add) {
      if (add >= int(this->mincut_group_revision_.size())) return;
      const int group_id = this->dcs_->group_id(add);
      const int group_size = this->dcs_->group_size(group_id);
      if (group_size == 1) return;

      const int F = this->used_flag_value();
      if (this->mincut_group_revision_[group_id] != F) {
        this->mincut_group_revision_[group_id] = F;
        this->mincut_group_counter_[group_id] = 0;
      }

      if (this->mincut_group_counter_[group_id] == 0) this->cross_other_mincut_count_++;
      this->mincut_group_counter_[group_id]++;
      if (this->mincut_group_counter_[group_id] == group_size) this->cross_other_mincut_count_--;
    };

    //s側の頂点とt側の頂点に分類する
    queue<int> q;
    vector<V> vs;
    const int F = used_flag_value();
    int one_side = 0;
    if (dz_.reason_for_finishing_bfs() == bi_dinitz::kQsIsEmpty) {
      //s側に属する頂点の親を新しいgroupに移動する
      q.push(s);
      vs.emplace_back(s);
      grouping_used_[s] = F;
      dcs_->create_new_group(s);
      while (!q.empty()) {
        V v = q.front(); q.pop();
        one_side++;
        for (auto& e : dz_.edges(v)) {
          const int cap = dz_.cap(e);
          if (cap == 0 || grouping_used_[dz_.to(e)] == F) continue;
          grouping_used_[dz_.to(e)] = F;
          q.push(dz_.to(e));
          vs.emplace_back(dz_.to(e));
          if (dcs_->is_same_group(t, dz_.to(e))) {
            dcs_->move_other_group(dz_.to(e), s); //tと同じgroupにいた頂点を、s側のgroupに移動
          } else {
            check_crossed_mincut(dz_.to(e));
          }
        }
      }
      gh_builder_->add_edge(s, t, cost, vs, dcs_);
    } else {
      //t側に属する頂点の親を,tに変更する
      q.push(t);
      vs.emplace_back(t);
      grouping_used_[t] = F;
      dcs_->create_new_group(t);
      while (!q.empty()) {
        V v = q.front(); q.pop();
        one_side++;
        for (auto& e : dz_.edges(v)) {
          const int cap = dz_.cap(dz_.rev(e));
          if (cap == 0 || grouping_used_[dz_.to(e)] == F) continue;
          grouping_used_[dz_.to(e)] = F;
          q.push(dz_.to(e));
          vs.emplace_back(dz_.to(e));
          if (dcs_->is_same_group(s, dz_.to(e))) {
            dcs_->move_other_group(dz_.to(e), t); //sと同じgroupにいた頂点を、t側のgroupに移動
          } else {
            check_crossed_mincut(dz_.to(e));
          }
        }
      }
      gh_builder_->add_edge(t, s, cost, vs, dcs_);
    }

    return one_side;
  }

  void contraction(const V s, const V t) {
    contraction_count_++;
    //gomory_hu algorithm
    //縮約後の頂点2つを追加する
    const int sside_new_vtx = dz_.n();
    const int tside_new_vtx = sside_new_vtx + 1;
    for (int _ = 0; _ < 2; _++) {
      dz_.add_vertex();
      grouping_used_.emplace_back();
      contraction_used_.emplace_back();
    }

    gh_builder_->contraction(s, t, sside_new_vtx, tside_new_vtx);

    queue<int> q;
    const int F = used_flag_value();
    int num_reconnected = 0; //枝を繋ぎ直した回数
    if (dz_.reason_for_finishing_bfs() == bi_dinitz::kQsIsEmpty) {
      q.push(s);
      contraction_used_[s] = F;
      while (!q.empty()) {
        V v = q.front(); q.pop();
        for (auto& e : dz_.edges(v)) {
          const int cap = dz_.cap(e);
          if (contraction_used_[dz_.to(e)] == F) continue;
          if (cap == 0) {
            if (grouping_used_[dz_.to(e)] != F) {
              //辺を上手に張り替える
              dz_.reconnect_edge(e, sside_new_vtx, tside_new_vtx);
              num_reconnected++;
            }
          } else {
            contraction_used_[dz_.to(e)] = F;
            q.push(dz_.to(e));
          }
        }
      }
    } else {
      q.push(t);
      contraction_used_[t] = F;
      while (!q.empty()) {
        V v = q.front(); q.pop();
        for (auto& e : dz_.edges(v)) {
          const int cap = dz_.cap(dz_.rev(e));
          if (contraction_used_[dz_.to(e)] == F) continue;
          if (cap == 0) {
            if (grouping_used_[dz_.to(e)] != F) {
              //辺を上手に張り替える
              dz_.reconnect_edge(e, tside_new_vtx, sside_new_vtx);
              num_reconnected++;
            }
          } else {
            contraction_used_[dz_.to(e)] = F;
            q.push(dz_.to(e));
          }
        }
      }
    }
    CHECK(num_reconnected == debug_last_max_flow_cost_); // 枝を繋ぎ直した回数 == maxflow
  }

public:

  separator(bi_dinitz& dz, disjoint_cut_set* dcs, unique_ptr<gomory_hu_tree_builder>& gh_builder)
    : dz_(dz), dcs_(dcs), gh_builder_(gh_builder),
    max_flow_times_(0), contraction_count_(0), grouping_used_(dz.n()), contraction_used_(dz.n()),
    mincut_group_counter_(dcs->node_num()), mincut_group_revision_(dcs->node_num()) {
  }

  void goal_oriented_bfs_init(const V goal) {
    dz_.goal_oriented_bfs_init(goal);
  }

  void mincut(V s, V t, bool enable_contraction = true) {
    if (dz_.edges(s).size() > dz_.edges(t).size()) swap(s, t);

    const int one_side = max_flow(s, t);
    if (enable_contraction) {
      const int other_side_estimated = dz_.n() - one_side;
      if (cross_other_mincut_count_ != 0) {
        fprintf(stderr, "(%d,%d) couldn't separate (crossed).\n", s, t);
      }

      const bool contract = cross_other_mincut_count_ == 0 &&
        min(one_side, other_side_estimated) >= FLAGS_cut_tree_contraction_lower_bound;
      if (contract) {
        contraction(s, t);
      }
    }

    // debug infomation
    debug_count_cut_size_all_time_[one_side]++;
    debug_count_cut_size_for_a_period_[one_side]++;
  }

  void output_debug_infomation() const {
    if (debug_count_cut_size_all_time_.size() > 10) {
      stringstream ss;
      for (auto& kv : debug_count_cut_size_all_time_) ss << "(" << kv.first << "," << kv.second << "), ";
      JLOG_ADD("separator.debug_count_cut_size_all_time_", ss.str());
      JLOG_ADD("separator.contraction_count", contraction_count_);
    }
  }


  void debug_verify() const {
    if (dcs_->node_num() > 10000) fprintf(stderr, "separator::debug_verify... ");
    union_find uf(dz_.n());
    for (int i = 0; i < dz_.n(); i++) for (const auto& to_edge : dz_.edges(i)) {
      uf.unite(i, dz_.to(to_edge));
    }
    for (int g = 0; g < dcs_->debug_group_num(); g++) {
      auto v = dcs_->get_group(g);
      CHECK(int(v.size()) == dcs_->group_size(g));
      for (int i = 0; i < int(v.size()) - 1; i++) {
        int u = v[i], x = v[i + 1];
        CHECK(uf.is_same(u, x));
      }
    }
    if (dcs_->node_num() > 10000) fprintf(stderr, "OK\n");
  }

  const bi_dinitz& get_bi_dinitz() const { return dz_; }
  const disjoint_cut_set* get_disjoint_cut_set() const { return dcs_; }

  const int contraction_count() { return contraction_count_; }

private:

  bi_dinitz& dz_;
  disjoint_cut_set* dcs_;
  unique_ptr<gomory_hu_tree_builder>& gh_builder_;

  int max_flow_times_; //maxflowを流した回数
  int contraction_count_; // contractionが呼ばれた回数
  vector<int> grouping_used_; // 'maxflowを流した後、mincutを求めるbfs'で使うused
  vector<int> contraction_used_; // 'mincutを元に、頂点縮約を行うbfs'で使うused

  vector<int> mincut_group_counter_;
  vector<int> mincut_group_revision_;
  int cross_other_mincut_count_; // 今回のmincutが、他のmincutと何回交わったか

  map<int, int> debug_count_cut_size_all_time_;
  map<int, int> debug_count_cut_size_for_a_period_;
  int debug_last_max_flow_cost_;
};
} // cut_tree_internal

//class cut_tree_with_2ecc
void cut_tree_with_2ecc::find_cuts_by_tree_packing(vector<pair<V, V>>& edges, disjoint_cut_set* dcs, const vector<int>& degree) {
  vector<int> current_parent(num_vertices_, -1);
  vector<int> current_weight(num_vertices_, -1);
  greedy_treepacking packing_base(edges, num_vertices_);

  auto set_solved = [&](V v, V parent, int weight) {
    current_parent[v] = parent;
    current_weight[v] = weight;
  };

  //degreeの最も大きな頂点をrootに
  int temp_root = 0;
  for (int v = 0; v < num_vertices_; v++) {
    if (degree[temp_root] < degree[v]) {
      temp_root = v;
    }
  }

  //次数2のcutを設定
  for (int v = 0; v < num_vertices_; v++) {
    if (v == temp_root) continue;
    if (degree[v] == 2) set_solved(v, temp_root, 2); //二重連結成分分解後なので自明なcut
  }

  if (FLAGS_cut_tree_enable_greedy_tree_packing) {
    //次数の高い頂点から順に、 一定回数 greedy tree packingを行って、flowの下界を求める
    const int iteration = min(FLAGS_cut_tree_try_greedy_tree_packing, num_vertices_);
    vector<int> idx(num_vertices_);
    for (int i = 0; i < num_vertices_; i++) idx[i] = i;
    partial_sort(idx.begin(), idx.begin() + iteration, idx.end(), [&degree](int l, int r) {
      return degree[l] > degree[r];
    });

    for (int trying = 0; trying < iteration; trying++) {
      const V v = idx[trying];
      if (degree[v] == 2) continue; // 自明なcutがある

      //debug infomation
      auto packing = packing_base;
      packing.arborescence_packing(v);

      if (current_parent[v] != -1) {
        current_parent[v] = v; // 閉路が出来上がるのを防ぐために、親を自分自身であると登録しておく
      }
      for (int to = 0; to < num_vertices_; to++) {
        if (to == v) continue;
        if (current_parent[to] != -1) continue;
        //tree packingの結果がdegreeと一致するなら、flowは流さなくてよい
        if (packing.inedge_count(to) == degree[to]) {
          set_solved(to, v, packing.inedge_count(to));
          // fprintf(stderr, "(%d, %d), cost = %d\n",v, to, packing.inedge_count(to));
        }
      }
    }

    //閉路が出来上がるのを防ぐためcurrent_parentに代入していた値を、元に戻す
    for (int trying = 0; trying < iteration; trying++) {
      const V v = idx[trying];
      if (current_parent[v] == v) current_parent[v] = -1;
    }

  }

  // cutの求まっていない頂点達について、gusfieldでcutを求める
  int pruned = 0;
  for (int v = 0; v < num_vertices_; v++) {
    if (v == temp_root) continue;
    if (current_parent[v] != -1) {
      // cutがもとまっている
      dcs->create_new_group(v);
      if (degree[v] != 2) {
        vector<V> vs;
        vs.emplace_back(v);
        gh_builder_->add_edge(v, current_parent[v], current_weight[v], vs, dcs);
      }
      pruned++;
    }
  }

  if (num_vertices_ > 10000) {
    JLOG_OPEN("prune") {
      JLOG_ADD("num_vs", num_vertices_);
      JLOG_ADD("pruned", pruned);
    }
  }
}

void cut_tree_with_2ecc::contract_degree2_vertices(vector<pair<V, V>>& edges, vector<int>& degree) {
  const int n = int(degree.size());
  vector<vector<int>> e(n);

  for (auto& uv : edges) {
    int u, v; tie(u, v) = uv;
    e[u].push_back(v);
    e[v].push_back(u);
  }

  int mi = distance(degree.begin(), max_element(degree.begin(), degree.end()));
  queue<int> que;
  unordered_set<int> visited;
  que.emplace(mi);
  visited.emplace(mi);
  while (!que.empty()) {
    int v = que.front();
    que.pop();
    for (int u : e[v]) {
      if (visited.count(u)) {
        continue;
      }
      visited.emplace(u);
      que.emplace(u);
      if (degree[u] == 2) {
        unordered_set<V> vs;
        vs.emplace(u);
        gh_builder_->add_degree2_edge(u, v);
      }
    }
  }

  for (int i = 0; i < n; i++) if (e[i].size() == 2) {
    int a = e[i][0], b = e[i][1];
    e[a].erase(find(e[a].begin(), e[a].end(), i));
    e[b].erase(find(e[b].begin(), e[b].end(), i));
    e[a].push_back(b);
    e[b].push_back(a);
    e[i].clear();
  }

  edges.clear();

  for (int i = 0; i < n; i++) {
    for (auto to : e[i]) {
      if (i < to) {
        edges.emplace_back(i, to);
      }
    }
  }
}

//次数の大きい頂点対をcutする
void cut_tree_with_2ecc::separate_high_degreepairs(separator* sep) {
  const disjoint_cut_set* dcs = sep->get_disjoint_cut_set();
  const bi_dinitz& dz = sep->get_bi_dinitz();

  vector<int> vtxs;
  for (int v = 0; v < num_vertices_; v++) {
    if (dcs->group_size(v) >= 2) vtxs.push_back(v);
  }
  const int tries = max(min(FLAGS_cut_tree_try_large_degreepairs, int(vtxs.size()) - 1), 0);
  partial_sort(vtxs.begin(), vtxs.begin() + tries, vtxs.end(), [&dz](V l, V r) {
    return dz.edges(l).size() > dz.edges(r).size();
  });

  int cut_large_degreecount = 0;
  for (int i = 1; i <= tries; i++) {
    for (int pari = 0; pari < i; pari++) {
      V s = vtxs[i];
      V t = vtxs[pari];
      if (!dcs->is_same_group(s, t)) continue;
      sep->mincut(s, t);
      cut_large_degreecount++;
    }
  }

  if (cut_large_degreecount > 0) {
    JLOG_PUT("separate_high_degreepairs.cut_large_degreecount", cut_large_degreecount);
    JLOG_PUT("separate_high_degreepairs.separated_count", sep->contraction_count());
  }
}

//隣接頂点同士を見て、まだ切れていなかったらcutする
void cut_tree_with_2ecc::separate_adjacent_pairs(separator* sep) {
  const bi_dinitz& dz = sep->get_bi_dinitz();
  const disjoint_cut_set* dcs = sep->get_disjoint_cut_set();

  for (int s = 0; s < num_vertices_; s++) {
    for (const auto& to_edge : dz.edges(s)) {
      const V t = dz.to(to_edge);
      if (!dcs->is_same_group(s, t)) continue;
      sep->mincut(s, t);
    }
  }
}

void cut_tree_with_2ecc::separate_all(separator* sep) {
  const disjoint_cut_set* dcs = sep->get_disjoint_cut_set();
  for (int group_id = 0; group_id < num_vertices_; group_id++) {
    while (dcs->has_two_elements(group_id)) {
      V s, t; tie(s, t) = dcs->get_two_elements(group_id);
      sep->mincut(s, t);
    }
  }
}

void cut_tree_with_2ecc::separate_near_pairs(separator* sep) {
  const disjoint_cut_set* dcs = sep->get_disjoint_cut_set();
  const bi_dinitz& dz = sep->get_bi_dinitz();

  vector<int> used(num_vertices_ * 2, -1);
  int used_revision = 0;

  for (int s = 0; s < num_vertices_; s++) {
    if (dcs->group_size(s) <= 1) continue;
    queue<V> q;
    q.push(s);
    used[s] = used_revision;
    for (int depth = 0; depth < FLAGS_cut_tree_separate_near_pairs_d; depth++) {
      const int loop_num = int(q.size());
      for (int _ = 0; _ < loop_num; _++) {
        const V v = q.front(); q.pop();
        for (auto& to_edge : dz.edges(v)) {
          const V t = dz.to(to_edge);
          if (used[t] == used_revision) continue;
          used[t] = used_revision;
          q.push(t);
          if (s != t && dcs->is_same_group(s, t)) {
            sep->mincut(s, t);
            used.resize(dz.n(), -1); // dinic中に頂点数が変わる場合がある
          }
        }
      }
    }
    used_revision++;
  }
}

//次数の最も高い頂点に対して、出来る限りの頂点からflowを流してmincutを求める
void cut_tree_with_2ecc::find_cuts_by_goal_oriented_search(separator* sep) {
  const bi_dinitz& dz = sep->get_bi_dinitz();
  const disjoint_cut_set* dcs = sep->get_disjoint_cut_set();

  int max_degreevtx = 0;
  for (int v = 0; v < num_vertices_; v++)
    if (dz.edges(max_degreevtx).size() < dz.edges(v).size()) max_degreevtx = v;

  sep->goal_oriented_bfs_init(max_degreevtx);
  for (int v = 0; v < num_vertices_; v++) {
    if (v == max_degreevtx) continue;
    if (!dcs->is_same_group(v, max_degreevtx)) continue;
    //graphの形状が変わると損なので、ここでは enable_contraction = false する
    sep->mincut(v, max_degreevtx, false);
  }
}

cut_tree_with_2ecc::cut_tree_with_2ecc(vector<pair<V, V>>&& edges, int num_vs) :
  num_vertices_(num_vs),
  gh_builder_(new gomory_hu_tree_builder(num_vs)) {
  vector<int> degree(num_vertices_);
  for (auto& e : edges) degree[e.first]++, degree[e.second]++;

  //次数2の頂点と接続を持つ辺を削除して、探索しやすくする
  JLOG_ADD_BENCHMARK_IF("time.contract_degree2_vertices", num_vertices_ > 10000) {
    contract_degree2_vertices(edges, degree);
  }

  unique_ptr<disjoint_cut_set> dcs(new disjoint_cut_set(num_vs));

  JLOG_ADD_BENCHMARK_IF("time.find_cuts_by_tree_packing", num_vertices_ > 10000) {
    find_cuts_by_tree_packing(edges, dcs.get(), degree);
  }

  //dinicの初期化
  bi_dinitz dz_base(std::move(edges), num_vs);

  separator sep(dz_base, dcs.get(), gh_builder_);

  if (FLAGS_cut_tree_enable_goal_oriented_search) {
    JLOG_ADD_BENCHMARK_IF("time.find_cuts_by_goal_oriented_search", num_vertices_ > 10000) {
      find_cuts_by_goal_oriented_search(&sep);
    }
  }

  // 次数の高い頂点対をcutする
  // グラフをなるべく2分するcutを見つけられると有用
  JLOG_ADD_BENCHMARK_IF("time.separate_high_degreepairs", num_vertices_ > 10000) {
    separate_high_degreepairs(&sep);
  }

  //まず隣接頂点対からcutしていく
  if (FLAGS_cut_tree_enable_adjacent_cut) {
    CHECK(FLAGS_cut_tree_separate_near_pairs_d >= 1);
    if (FLAGS_cut_tree_separate_near_pairs_d == 1) {
      JLOG_ADD_BENCHMARK_IF("time.separate_adjacent_pairs", num_vertices_ > 10000) {
        separate_adjacent_pairs(&sep);
      }
    } else {
      JLOG_ADD_BENCHMARK_IF("time.separate_near_pairs", num_vertices_ > 10000) {
        separate_near_pairs(&sep);
      }
    }
  }

  // sep.debug_verify();

  // 残った頂点groupをcutする、gomory_hu treeの完成
  JLOG_ADD_BENCHMARK_IF("time.separate_all", num_vertices_ > 10000) {
    separate_all(&sep);
  }

  sep.output_debug_infomation();

  gh_builder_->build();
}

cut_tree_with_2ecc::~cut_tree_with_2ecc() = default;

int cut_tree_with_2ecc::query(V u, V v) const {
  return gh_builder_->query(u, v);
}

const vector<pair<V, int>>& cut_tree_with_2ecc::parent_weight() const {
  return gh_builder_->parent_weight();
}
} // namespace agl
