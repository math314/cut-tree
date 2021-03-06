#pragma once
#include "gomory_hu_tree_builder.h"

namespace agl {
namespace cut_tree_internal {
namespace plain_gomory_hu {
class dinitz_separator {

  const int used_flag_value() const {
    return max_flow_times_;
  }

  // 一定期間置きに進捗を出力する
  void print_progress_at_regular_intervals(V s, V t, int cost) {
    if (max_flow_times_ % 10000 == 0) {
      std::stringstream ss;
      ss << "max_flow_times_ = " << max_flow_times_ << ", (" << s << "," << t << ") cost = " << cost;
      JLOG_ADD("separator.progress", ss.str());

      fprintf(stderr, "cut details : ");
      for(auto& kv : debug_count_cut_size_for_a_period_) fprintf(stderr, "(%d,%d), ", kv.first, kv.second);
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
      if(add >= int(this->mincut_group_revision_.size())) return ;
      const int group_id = this->dcs_.group_id(add);
      const int group_size = this->dcs_.group_size(group_id);
      if(group_size == 1) return;

      const int F = this->used_flag_value();
      if(this->mincut_group_revision_[group_id] != F){
        this->mincut_group_revision_[group_id] = F;
        this->mincut_group_counter_[group_id] = 0;
      }

      if(this->mincut_group_counter_[group_id] == 0) this->cross_other_mincut_count_++;
      this->mincut_group_counter_[group_id]++;
      if(this->mincut_group_counter_[group_id] == group_size) this->cross_other_mincut_count_--;
    };

    //s側の頂点とt側の頂点に分類する
    std::queue<int> q;
    std::vector<int> vs;
    const int F = used_flag_value();
    int one_side = 0;
    //s側に属する頂点の親を新しいgroupに移動する
    q.push(s);
    vs.emplace_back(s);
    grouping_used_[s] = F;
    dcs_.create_new_group(s);
    while (!q.empty()) {
      V v = q.front(); q.pop();
      one_side++;
      for (auto& e : dz_.edges(v)) {
        const int cap = e.cap_;
        if (cap == 0 || grouping_used_[e.to_] == F) continue;
        grouping_used_[e.to_] = F;
        q.push(e.to_);
        vs.emplace_back(e.to_);
        if (dcs_.is_same_group(t, e.to_)) {
          dcs_.move_other_group(e.to_, s); //tと同じgroupにいた頂点を、s側のgroupに移動
        } else {
          check_crossed_mincut(e.to_);
        }
      }
    }
    gh_builder_.add_edge(s, t, cost, vs, dcs_); //cutした結果をgomory_hu treeの枝を登録

    return one_side;
  }

  void contraction(const V s,const V t) {
    contraction_count_++;
    //gomory_hu algorithm
    //縮約後の頂点2つを追加する
    const int sside_new_vtx = dz_.n();
    const int tside_new_vtx = sside_new_vtx + 1;
    for(int _ = 0; _ < 2; _++) {
      dz_.add_vertex();
      grouping_used_.emplace_back();
      contraction_used_.emplace_back();
    }

    gh_builder_.contraction(s, t, sside_new_vtx, tside_new_vtx);

    std::queue<int> q;
    const int F = used_flag_value();
    int num_reconnected = 0; //枝を繋ぎ直した回数
    q.push(s);
    contraction_used_[s] = F;
    while (!q.empty()) {
      V v = q.front(); q.pop();
      for (auto& e : dz_.edges(v)) {
        const int cap = e.cap_;
        if (contraction_used_[e.to_] == F) continue;
        if (cap == 0) {
          if (grouping_used_[e.to_] != F) {
            //辺を上手に張り替える
            dz_.reconnect_edge(e, sside_new_vtx, tside_new_vtx);
            num_reconnected++;
          }
        } else {
          contraction_used_[e.to_] = F;
          q.push(e.to_);
        }
      }
    }
    CHECK(num_reconnected == debug_last_max_flow_cost_); // 枝を繋ぎ直した回数 == maxflow
  }

public:

  dinitz_separator(dinitz& dz, disjoint_cut_set& dcs, gomory_hu_tree_builder& gh_builder) 
    : dz_(dz), dcs_(dcs), gh_builder_(gh_builder),
      max_flow_times_(0), contraction_count_(0), grouping_used_(dz.n()), contraction_used_(dz.n()), 
      mincut_group_counter_(dcs.node_num()), mincut_group_revision_(dcs.node_num()) {
  }

  void mincut(V s, V t, bool enable_contraction = true) {
    if (dz_.edges(s).size() > dz_.edges(t).size()) std::swap(s, t);

    const int one_side = max_flow(s, t);
    if (enable_contraction) {
      const int other_side_estimated = dz_.n() - one_side;
      if(cross_other_mincut_count_ != 0) {
        fprintf(stderr, "(%d,%d) couldn't separate (crossed).\n", s, t);
      }

      const bool contract = cross_other_mincut_count_ == 0 &&
        std::min(one_side, other_side_estimated) >= 2/*FLAGS_contraction_lower_bound*/;
      if(contract) {
        contraction(s, t);
      }
    }

    // debug infomation
    debug_count_cut_size_all_time_[one_side]++;
    debug_count_cut_size_for_a_period_[one_side]++;
  }

  void output_debug_infomation() const {
    if (debug_count_cut_size_all_time_.size() > 10) {
      std::stringstream ss;
      for (auto& kv : debug_count_cut_size_all_time_) ss << "(" << kv.first << "," << kv.second << "), ";
      JLOG_ADD("separator.debug_count_cut_size_all_time_", ss.str());
      JLOG_ADD("separator.contraction_count", contraction_count_);
    }
  }


  void debug_verify() const {
    if(dcs_.node_num() > 10000) fprintf(stderr, "separator::debug_verify... ");
    union_find uf(dz_.n());
    for(int i = 0; i < dz_.n(); i++) for (const auto& to_edge : dz_.edges(i)) {
      uf.unite(i, to_edge.to_);
    }
    for(int g = 0; g < dcs_.debug_group_num(); g++) {
      auto v = dcs_.get_group(g);
      CHECK(int(v.size()) == dcs_.group_size(g));
      for(int i = 0; i < int(v.size()) - 1; i++) {
        int u = v[i], x = v[i + 1];
        CHECK(uf.is_same(u, x));
      }
    }
    if(dcs_.node_num() > 10000) fprintf(stderr, "OK\n");
  }

  const dinitz& get_dinitz() const { return dz_; }
  const disjoint_cut_set& get_disjoint_cut_set() const { return dcs_; }

  const int contraction_count() { return contraction_count_; }

private:

  dinitz& dz_;
  disjoint_cut_set& dcs_;
  gomory_hu_tree_builder& gh_builder_;

  int max_flow_times_; //maxflowを流した回数
  int contraction_count_; // contractionが呼ばれた回数
  std::vector<int> grouping_used_; // 'maxflowを流した後、mincutを求めるbfs'で使うused
  std::vector<int> contraction_used_; // 'mincutを元に、頂点縮約を行うbfs'で使うused

  std::vector<int> mincut_group_counter_;
  std::vector<int> mincut_group_revision_;
  int cross_other_mincut_count_; // 今回のmincutが、他のmincutと何回交わったか

  std::map<int, int> debug_count_cut_size_all_time_;
  std::map<int, int> debug_count_cut_size_for_a_period_;
  int debug_last_max_flow_cost_;
};

class gomory_hu_dinitz {
  void separate_all(dinitz_separator& sep) {
    const disjoint_cut_set& dcs = sep.get_disjoint_cut_set();
    for(int group_id = 0; group_id < num_vertices_; group_id++) {
      while (dcs.has_two_elements(group_id)) {
        V s, t; std::tie(s, t) = dcs.get_two_elements(group_id);
        sep.mincut(s, t);
      }
    }
  }

public:

  gomory_hu_dinitz(std::vector<std::pair<V, V>>&& edges, int num_vs) :
    num_vertices_(num_vs),
    gh_builder_(num_vs) {

    if(num_vs > 10000) fprintf(stderr, "gomory_hu_dinitz::constructor start : memory %ld MB\n", jlog_internal::get_memory_usage() / 1024);

    disjoint_cut_set dcs(num_vs);

    //dinicの初期化
    G g(edges);
    dinitz dz_base(g);
    for (auto& e : edges) {
      dz_base.add_undirected_edge(e.first, e.second, 1);
    }
    if(num_vs > 10000) fprintf(stderr, "gomory_hu_dinitz::dinitz after init : memory %ld MB\n", jlog_internal::get_memory_usage() / 1024);

    dinitz_separator sep(dz_base, dcs, gh_builder_);

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
