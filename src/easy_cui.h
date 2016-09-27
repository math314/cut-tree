/*
 * 試作やベンチマーク目的でグラフを扱う CUI アプリケーションを作りたい場合，
 * こいつを使うと便利です．ズボラな人用．必ず本体から 1 度のみ include されるようにして下さい．
 */

#pragma once
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <stack>
#include <queue>
#include <set>
#include <map>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cmath>
#include <cassert>
#include <unordered_set>
#include <unordered_map>
#include "agl.h"
using namespace std;
using namespace agl;

DEFINE_string(type, "auto", "auto, tsv, agl, built_in, gen");
DEFINE_string(graph, "-", "input graph");
DEFINE_bool(force_undirected, false, "Automatically add reverse edges?");

template<typename GraphType = G>
string guess_type() {
  FILE* fp = fopen(FLAGS_graph.c_str(), "r");
  if(fp == NULL) return "built_in";
  
  char buf[20];
  if(fgets(buf, sizeof(buf), fp) == NULL){
    string mes = "An error occured on read '" + FLAGS_graph + "'.";
    FAIL_MSG(mes.c_str());
  }
  fclose(fp);

  string header(buf);
  if(header == "AGL_BINARY\n") {
    return "agl";
  } else {
    return "tsv";
  }

}

template<typename GraphType = G>
GraphType easy_cui_init(int argc, char **argv) {
  JLOG_INIT(&argc, argv);
  google::ParseCommandLineFlags(&argc, &argv, true);

  if(FLAGS_type == "auto")
    FLAGS_type = guess_type();

  //agl形式は隣接リストを介すると効率が悪いため、直接GraphTypeを生成する
  if (FLAGS_type == "agl") {
    auto g = read_graph_binary<G>(FLAGS_graph.c_str());
    if(FLAGS_force_undirected){
      g = GraphType(make_undirected(g.edge_list()));
    }
    pretty_print(g);
    return g;
  }

  G::edge_list_type es;
  if (FLAGS_type == "tsv") {
    es = read_edge_list_tsv(FLAGS_graph.c_str());
  } else if (FLAGS_type == "built_in") {
    es = built_in_edge_list(FLAGS_graph.c_str());
  } else if (FLAGS_type == "gen") {
    istringstream iss(FLAGS_graph);
    string family;
    iss >> family;
    if (family == "barbell") {
      V n;
      if (!(iss >> n)) n = 4;
      es = generate_barbell(n);
    } else if (family == "grid") {
      size_t r, c;
      if (!(iss >> r)) r = 4;
      if (!(iss >> c)) c = r;
      es = generate_grid(r, c);
    } else if (family == "erdos_renyi") {
      V n;
      double d;
      if (!(iss >> n)) n = 10;
      if (!(iss >> d)) d = 3.0;
      es = generate_erdos_renyi(n, d);
    } else if (family == "random_planar") {
      V n;
      size_t e;
      if (!(iss >> n)) n = 10;
      if (!(iss >> e)) e = 25;
      es = generate_random_planar(n, e);
    } else if (family == "cycle") {
      V n;
      if (!(iss >> n)) n = 10;
      es = generate_cycle(n);
    } else if (family == "ba") {
      V n, m;
      if (!(iss >> n)) n = 10;
      if (!(iss >> m)) m = 5;
      es = generate_ba(n, m);
    } else if (family == "dms") {
      V n, m, k0;
      if (!(iss >> n)) n = 10;
      if (!(iss >> m)) m = 5;
      if (!(iss >> k0)) k0 = -2;
      es = generate_dms(n, m, k0);
    } else if (family == "hk") {
      V n, m;
      double p;
      if (!(iss >> n)) n = 10;
      if (!(iss >> m)) m = 5;
      if (!(iss >> p)) p = 0.5;
      es = generate_hk(n, m, p);
    } else if (family == "ws") {
      V n, d;
      double p;
      if (!(iss >> n)) n = 10;
      if (!(iss >> d)) d = 4;
      if (!(iss >> p)) p = 0.5;
      es = generate_ws(n, d, p);
    } else if (family == "kronecker") {
      int scale, n;
      size_t avg_deg;
      if (!(iss >> scale)) scale = 5;
      if (!(iss >> avg_deg)) avg_deg = 16;
      vector<vector<double>> mat;
      if (!(iss >> n)) {
        n = 2;
        mat = vector<vector<double>>(n, vector<double>(n));
        mat[0][0] = 0.57;
        mat[0][1] = 0.19;
        mat[1][0] = 0.19;
        mat[1][1] = 0.05;
      }
      for (int i = 0; i < n * n; ++i) {
        if (!(iss >> mat[i / n][i % n])) mat[i / n][i % n] = 1.0 / (n * n);
      }
      es = generate_kronecker(scale, avg_deg, mat);
    } else if (family == "flower") {
      V required, u, v;
      if (!(iss >> required)) required = 44;
      if (!(iss >> u)) u = 2;
      if (!(iss >> v)) v = 2;
      es = generate_uv_flower(required, u, v);
    } else if (family == "shm") {
      V required_num, initial_num;
      int t;
      double P;
      if (!(iss >> required_num)) required_num = 101;
      if (!(iss >> initial_num)) initial_num = 5;
      if (!(iss >> t)) t = 2;
      if (!(iss >> P)) P = 0.0;
      es = generate_shm(required_num, initial_num, t, P);
    } else {
      FAIL_MSG("Unknown generator family: " + family);
    }
  }

  if (FLAGS_force_undirected) es = make_undirected(es);

  GraphType g(es);
  pretty_print(g);
  return g;
}
