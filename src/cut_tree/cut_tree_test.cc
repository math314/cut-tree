#include "cut_tree.h"
#include <gtest/gtest.h>

#include <sstream>
#include <vector>
#include <map>

using namespace agl;
using namespace agl::cut_tree_internal;
using namespace std;

namespace {

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

template<class cut_tree_t>
void output_tree_varify(G g) {
  cut_tree_t gf(g);
  stringstream ss;
  gf.print_gomory_hu_tree(ss);
  auto query = cut_tree_query_handler::from_file(ss);
  const int n = g.num_vertices();
  for(int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      int a1 = gf.query(i, j);
      int a2 = query.query(i, j);
      ASSERT_EQ(a1, a2);
    }
  }
}

typedef testing::Types<cut_tree, gomory_hu_bi_dinitz, gomory_hu_dinitz> CutTreeTestTypes;

template<typename T>
class cut_tree_test : public testing::Test {};
TYPED_TEST_CASE(cut_tree_test, CutTreeTestTypes);

TYPED_TEST(cut_tree_test, output_tree) {
  using cut_tree_t = TypeParam;
  output_tree_varify<cut_tree_t>(to_directed_graph(built_in_graph("karate_club")));
  output_tree_varify<cut_tree_t>(to_directed_graph(built_in_graph("dolphin")));
  // output_tree_varify<cut_tree_t>(to_directed_graph(built_in_graph("ca_grqc")));
    for (int trial = 0; trial < 3; ++trial) {
      V M = 3;
      V N = M + agl::random(1000);
      auto es = generate_ba(N, M);
      G g = to_directed_graph(G(es));
      output_tree_varify<cut_tree_t>(g);
    }
}

TEST(cut_tree_test, dinitz_eq_bi_dinitz) {
  auto check = [](vector<pair<V,V>>&& es) {
    G g = to_directed_graph(G(es));
    dinitz dz1(g);
    bi_dinitz dz2(g);
    for (int i = 0; i < 1000; i++) {
      V s = agl::random() % g.num_vertices();
      V t = agl::random() % (g.num_vertices() - 1);
      if (s <= t) t++;
      int a = dz1.max_flow(s, t);
      int b = dz2.max_flow(s, t);
      ASSERT_EQ(a, b);
    }
  };
  
  check(built_in_graph("ca_grqc").edge_list());
  check(generate_path(100));
  for (int trial = 0; trial < 3; ++trial) {
    check(generate_erdos_renyi(500,10));
  }
  for (int trial = 0; trial < 3; ++trial) {
    check(generate_grid(30,3));
    check(generate_grid(10,10));
  }
  check(generate_barbell(100));
  check(generate_cycle(100));
  for (int trial = 0; trial < 3; ++trial) {
    V M = 3;
    V N = M + agl::random(1000);
    check(generate_ba(N, M));
  }
}

TEST(cut_tree_test, cut_tree_eq_bi_dinitz) {
  auto check = [](vector<pair<V,V>>&& es) {
    G g = to_directed_graph(G(es));
    const int n = g.num_vertices();
    bi_dinitz dt(g);
    //cut_tree constructor break down 'g'
    cut_tree ct(g);
    ASSERT_EQ(g.num_vertices(), 0);
    for (int i = 0; i < 10000; i++) {
      V s = agl::random() % n;
      V t = agl::random() % (n - 1);
      if (s <= t) t++;
      int a = ct.query(s, t);
      int b = dt.max_flow(s, t);
      ASSERT_EQ(a, b);
    }
  };
  check(built_in_graph("ca_grqc").edge_list());
  check(generate_path(100));
  for (int trial = 0; trial < 3; ++trial) {
    check(generate_erdos_renyi(500,10));
  }
  for (int trial = 0; trial < 3; ++trial) {
    check(generate_grid(30,3));
    check(generate_grid(10,10));
  }
  check(generate_barbell(100));
  check(generate_cycle(100));
  for (int trial = 0; trial < 3; ++trial) {
    V M = 3;
    V N = M + agl::random(1000);
    check(generate_ba(N, M));
  }
}

TYPED_TEST(cut_tree_test, corner_case_small_graph) {
  using cut_tree_t = TypeParam;
  for(int vertex = 0; vertex <= 2; vertex++){
    vector<pair<V,V>> es;
    for(int i = 0; i < vertex; i++) {
      for(int j = i + 1; j < vertex; j++) {
        es.emplace_back(i, j);
      }
    }
    G g(es, vertex);
    cut_tree_t ct(g);
    stringstream ss;
    ct.print_gomory_hu_tree(ss);
  }
}

}
