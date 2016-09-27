#include <easy_cui.h>

DEFINE_string(cut_tree_input_path, "", "");
DEFINE_string(cut_tree_output_path, "", "");

class union_find_with_size {
private:
  int n;
  vector<int> a;
public:
  union_find_with_size(int n) : n(n), a(n, -1) {}
  int find(int x) { return a[x] < 0 ? x : (a[x] = find(a[x])); }
  bool same(int x, int y) { return find(x) == find(y); }
  bool unite(int x, int y) {
    x = find(x), y = find(y);
    if (x == y) return false;
    if (a[x] > a[y]) swap(x, y);
    a[x] += a[y];
    a[y] = x;
    n--;
    return true;
  }
  int size(int x) { return -a[find(x)]; }
};

map<int, vector<pair<V, V>>> load(const string& path) {
  ifstream is(path.c_str());
  CHECK_MSG(!is.fail(), "input gomory_hu tree file cannot open.");

  map<int, vector<pair<V, V>>> ret;
  int s, v, weight;
  while (is >> s >> v >> weight) {
    ret[weight].emplace_back(s, v);
  }
  return ret;
}

void save(const string& path, const map<int, long long>& connectivity_distribution_sum) {
  ofstream os(path.c_str());
  CHECK_MSG(!os.fail(), "output connectivity distribution file cannot open.");
  for (const auto& cd : connectivity_distribution_sum) {
    os << cd.first << " " << cd.second << endl;
  }
}

string get_cut_tree_output_path() {
  string cut_tree_output_path = FLAGS_cut_tree_output_path;
  if (cut_tree_output_path == "") {
    cut_tree_output_path = FLAGS_cut_tree_input_path;
    auto idx = cut_tree_output_path.find_last_of('.');
    if (idx != string::npos) {
      cut_tree_output_path = cut_tree_output_path.substr(0, idx);
    }
    cut_tree_output_path += ".cdist";
  }
  return cut_tree_output_path;
}

int main(int argc, char** argv) {
  google::ParseCommandLineFlags(&argc, &argv, true);
  const string cut_tree_output_path = get_cut_tree_output_path();

  auto gomory_hu_tree = load(FLAGS_cut_tree_input_path);
  //graph size
  const int n = [&gomory_hu_tree]() {
    int edges = 0;
    for (const auto& wuv : gomory_hu_tree) edges += (int)wuv.second.size();
    return edges + 1; // on tree, num_vertices = edge + 1
  }();

  union_find_with_size uf(n);

  map<int, long long> connectivity_distribution_sum;
  long long cur = 0;
  for (auto it = gomory_hu_tree.rbegin(); it != gomory_hu_tree.rend(); ++it) {
    const int w = it->first;
    for (const auto& uv : it->second) {
      V u, v; tie(u, v) = uv;
      CHECK(!uf.same(u, v));
      const long long us = uf.size(u);
      const long long vs = uf.size(v);
      cur += us * vs;
      uf.unite(u, v);
    }
    connectivity_distribution_sum[w] = cur;
  }

  save(cut_tree_output_path, connectivity_distribution_sum);
}
