#include <easy_cui.h>
#include "cut_tree_query_handler.h"

DEFINE_int32(cut_tree_num_query, 10000000, "");
DEFINE_int64(cut_tree_node_pair_random_seed, 922337203685477583LL, "");
DEFINE_string(cut_tree_path, "", "");

void from_file() {
  cut_tree_query_handler tq;
  JLOG_PUT_BENCHMARK("initialize_time") {
    tq = cut_tree_query_handler::from_file(FLAGS_cut_tree_path);
  }

  agl::random_type random(FLAGS_cut_tree_node_pair_random_seed);
  JLOG_PUT_BENCHMARK("query_time") {
    for (int i = 0; i < FLAGS_cut_tree_num_query; i++) {
      V s = random() % tq.num_vertices();
      V t = random() % (tq.num_vertices() - 1);
      if (s <= t) t++;
      tq.query(s, t);
    }
  }
}

int main(int argc, char** argv) {
  JLOG_INIT(&argc, argv);
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  from_file();

  return 0;
}
