#include <easy_cui.h>
#include "cut_tree_query_handler.h"

DEFINE_string(cut_tree_path, "", "");
DEFINE_string(query_path, "", "");
DEFINE_string(output_path, "", "");

void from_file() {
  cut_tree_query_handler tq;
  tq = cut_tree_query_handler::from_file(FLAGS_cut_tree_path);
  std::ifstream ifs(FLAGS_query_path.c_str());
  ostream* ost;
  if (FLAGS_output_path == "") {
    ost = new iostream(cout.rdbuf());
  } else {
    ost = new ofstream(FLAGS_output_path.c_str(), ios_base::out);
  }
  V s, t;
  while (ifs >> s >> t) {
    *ost << tq.query(s, t) << endl;
  }
}

int main(int argc, char** argv) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  from_file();

  return 0;
}
