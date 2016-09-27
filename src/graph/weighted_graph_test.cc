#include "graph.h"
#include "gtest/gtest.h"
using namespace agl;
using namespace std;

using graph_type = weighted_graph<double>;

TEST(weighted_graph_test, instantiation) {
  graph_type g;

  g.assign(graph_type::edge_list_type{
                                              {0, {1, 1.5}},
                                              {1, {2, 20}},
                                              {2, {3, 5}},
                                              {3, {1, 0.01}},
                                      });

  pretty_print(g);
  pretty_print(g);

}
TEST(weighted_graph, read_graph_tsv) {
  graph_type g(add_random_weight<graph_type>(generate_erdos_renyi(10, 2)));
  pretty_print(g);
}
