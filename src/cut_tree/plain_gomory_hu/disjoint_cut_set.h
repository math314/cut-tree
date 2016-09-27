#pragma once

namespace agl {
namespace cut_tree_internal {
namespace plain_gomory_hu {
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
    for(int i = 0; i < n - 1; i++) {
      nodes[i].nt = i + 1;
      nodes[i + 1].pv = i;
    }
    for(int i = 0; i < n; i++) nodes[i].root = 0;
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

  std::pair<int, int> get_two_elements(int group_id) const {
    const int rt = root[group_id];
    CHECK(rt != -1);
    const int nxt = nodes[rt].nt;
    if (nxt == -1) return std::make_pair(-1, -1);
    return std::make_pair(rt, nxt);
  }

  bool has_two_elements(int group_id) const {
    auto uv = get_two_elements(group_id);
    return uv.first != -1;
  }

  std::vector<int> get_group(int group_id) const {
    std::vector<int> ret;
    int cur = root[group_id];
    while(cur != -1){
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
  std::vector<int> root;
  std::vector<Node> nodes;
  int group_num;
  std::vector<int> group_size_;
};
} // namespace plain_gomory_hu
} // namespace cut_tree_internal
} // namespace agl
