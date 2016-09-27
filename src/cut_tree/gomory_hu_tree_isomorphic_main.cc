#include <easy_cui.h>

DEFINE_string(s1, "s1.txt", "gomory-hu tree 1");
DEFINE_string(s2, "s2.txt", "gomory-hu tree 2");

#include <cut_tree/bi_dinitz.h>

map<int, vector<pair<V, V>>> load(const string& path) {
  FILE* fp = fopen(path.c_str(), "r");
  map<int, vector<pair<V, V>>> ret;
  int s, v, weight;
  while (~fscanf(fp, "%d%d%d", &s, &v, &weight)) {
    ret[weight].emplace_back(s, v);
  }
  return ret;
}

void check(G& g, map<int, vector<pair<V, V>>>& l, map<int, vector<pair<V, V>>>& r) {
  using ull = unsigned long long;
  int n = 1;
  for (auto& kv : l) n += int(kv.second.size());
  vector<ull> zobrist_hash(n);
  for (int i = 0; i < n; i++) zobrist_hash[i] = (ull(agl::random()) << 32) | ull(agl::random());

  union_find ufl(n), ufr(n);
  vector<ull> hl = zobrist_hash, hr = zobrist_hash;

  vector<int> weight;
  for (auto& kv : l) weight.push_back(kv.first);
  reverse(weight.begin(), weight.end());

  auto on_mismatched = [&](V u) {
    set<int> onlyl, onlyr, intersect;
    for (int a = 0; a < n; a++) if (ufl.is_same(u, a)) onlyl.insert(a);
    for (int a = 0; a < n; a++) if (ufr.is_same(u, a)) {
      if (onlyl.count(a)) onlyl.erase(a), intersect.insert(a);
      else onlyr.insert(a);
    }
    int cnt;
    printf("left: ");
    cnt = 0;
    for (auto x : onlyl) {
      printf("%d, ", x);
      cnt++;
      if (cnt >= 10) { printf("..."); break; }
    }
    puts("");
    printf("right: ");
    cnt = 0;
    for (auto x : onlyr) {
      printf("%d, ", x);
      cnt++;
      if (cnt >= 10) { printf("..."); break; }
    }
    puts("");
    printf("intersect: ");
    cnt = 0;
    for (auto x : intersect) {
      printf("%d, ", x);
      cnt++;
      if (cnt >= 10) { printf("..."); break; }
    }
    puts("");

    auto check_cost = [&](int x, map<int, vector<pair<int, int>>>& mp, string name) {
      bi_dinitz dc(g);
      for (auto& kv : mp) {
        int w = kv.first;
        for (auto& uv : kv.second) {
          int u, v; tie(u, v) = uv;
          if (u != x && v != x) continue;
          dc.reset_graph();
          int cost = dc.max_flow(u, v);
          if (cost != w) {
            printf("(%d,%d) invalid cost. answer = %d, but '%s' indicate %d\n", u, v, cost, name.c_str(), w);
            exit(0);
          }
        }
      }
    };

    for (auto x : onlyl) {
      check_cost(x, l, FLAGS_s1);
      check_cost(x, r, FLAGS_s2);
    }

    for (auto y : onlyr) {
      check_cost(y, l, FLAGS_s1);
      check_cost(y, r, FLAGS_s2);
    }
  };

  for (const int w : weight) {
    if (l[w].size() != r[w].size()) printf("*");
    printf("w = %d, L : %d , R : %d", w, int(l[w].size()), int(r[w].size()));
    puts("");
  }

  for (const int w : weight) {
    for (auto& uv : l[w]) {
      int u, v; tie(u, v) = uv;
      u = ufl.root(u), v = ufl.root(v);
      CHECK(u != v);
      ufl.unite(u, v);
      ull new_hash = hl[u] ^ hl[v];
      hl[ufl.root(u)] = new_hash;
    }

    for (auto& uv : r[w]) {
      int u, v; tie(u, v) = uv;
      u = ufr.root(u), v = ufr.root(v);
      CHECK(u != v);
      ufr.unite(u, v);
      ull new_hash = hr[u] ^ hr[v];
      hr[ufr.root(u)] = new_hash;
    }

    for (auto& uv : l[w]) {
      int u, v; tie(u, v) = uv;
      if (hl[ufl.root(u)] != hr[ufr.root(u)]) {
        on_mismatched(u);
        exit(-1);
      }
      if (hl[ufl.root(v)] != hr[ufr.root(v)]) {
        on_mismatched(v);
        exit(-1);
      }
    }
    for (auto& uv : r[w]) {
      int u, v; tie(u, v) = uv;
      if (hl[ufl.root(u)] != hr[ufr.root(u)]) {
        on_mismatched(u);
        exit(-1);
      }
      if (hl[ufl.root(v)] != hr[ufr.root(v)]) {
        on_mismatched(v);
        exit(-1);
      }
    }
  }
}

G to_directed_graph(G&& g) {
  vector<pair<V, V>> ret;
  for (auto& e : g.edge_list()) {
    if (e.first < to(e.second)) ret.emplace_back(e.first, to(e.second));
    else if (to(e.second) < e.first) ret.emplace_back(to(e.second), e.first);
  }
  sort(ret.begin(), ret.end());
  ret.erase(unique(ret.begin(), ret.end()), ret.end());
  return G(ret);
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  auto g = to_directed_graph(read_graph_binary<G>(FLAGS_graph.c_str()));
  auto gomory_hu1 = load(FLAGS_s1);
  auto gomory_hu2 = load(FLAGS_s2);
  check(g, gomory_hu1, gomory_hu2);

  fprintf(stderr, "ok. %s == %s\n", FLAGS_s1.c_str(), FLAGS_s2.c_str());
  return 0;
}
