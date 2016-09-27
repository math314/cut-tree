#include "generator.h"
#include <set>
#include "base/base.h"
using namespace std;

namespace agl {
unweighted_edge_list generate_path(V num_vertices) {
  unweighted_edge_list es;
  for (V v = 0; v + 1 < num_vertices; ++v) {
    es.emplace_back(v, v + 1);
  }
  return es;
}

unweighted_edge_list generate_erdos_renyi(V num_vertices, double avg_deg) {
  avg_deg = min(avg_deg, static_cast<double>(max(0, num_vertices - 1)));
  set<pair<V, V>> es;
  std::uniform_int_distribution<V> rng(0, num_vertices - 1);
  while (es.size() < num_vertices * avg_deg) {
    V u = rng(agl::random), v = rng(agl::random);
    if (u == v) continue;
    es.insert(make_pair(u, v));
  }
  return vector<pair<V, V>>(es.begin(), es.end());
};

unweighted_edge_list generate_grid(size_t num_rows, size_t num_cols) {
  auto vid = [=](int i, int j) { return i * num_cols + j; };

  vector<pair<V, V>> es;
  for (auto i : make_irange(num_rows)) {
    for (auto j : make_irange(num_cols)) {
      if (j + 1 < num_cols) es.emplace_back(vid(i, j), vid(i, j + 1));
      if (i + 1 < num_rows) es.emplace_back(vid(i, j), vid(i + 1, j));
    }
  }
  return es;
}

unweighted_edge_list generate_barbell(V size_clique) {
  unweighted_edge_list out;
  for (V i : make_irange(2)) {
    for (V v : make_irange(size_clique)) {
      for (V u : make_irange(v)) {
        out.emplace_back(i * size_clique + u, i * size_clique + v);
      }
    }
  }
  out.emplace_back(0, size_clique);
  return out;
}

unweighted_edge_list generate_cycle(V num_vertices) {
  unweighted_edge_list es;
  if (num_vertices < 2) {
    return es;
  }
  for (V v = 0; v + 1 < num_vertices; ++v) {
    es.emplace_back(v, v + 1);
  }
  es.emplace_back(num_vertices - 1, 0);
  return es;
}

/**
 * Generate a random scale-free network by the Barabasi-Albert (BA) model.
 * The degree distribution resulting from the BA model is scale free
 * with power-law coefficient &gamma; = 3.
 * \param initial_num is a number of nodes of the initial connected network.
 * \param final_num is a number of finally generated network.
 */
unweighted_edge_list generate_ba(V final_num, V initial_num) {
  CHECK(initial_num >= 2);
  unweighted_edge_list es;
  for (int v = 0; v < initial_num; ++v) {
    for (int u = 0; u < v; ++u) {
      es.emplace_back(u, v);
    }
  }

  for (int v = initial_num; v < final_num; ++v) {
    set<V> next;
    std::uniform_int_distribution<size_t> rng(0, es.size() - 1);
    while (next.size() < (size_t)initial_num) {
      size_t e = rng(agl::random);
      V u = agl::random() % 2 ? es[e].first : es[e].second;
      next.insert(u);
    }
    for (auto u : next) {
      es.emplace_back(u, v);
    }
  }
  return es;
}

/**
 * Generate a random scale-free network by the Dorogovtsev-Mendes-Samukhin (DMS) model.
 * Depending on K0, the power-law coefficient &gamma; takes values from 2 to &infin;.
 * \param initial_num is a number of nodes of the initial connected network.
 * \param final_num is a number of finally generated network.
 * \param K0 is a constant value and forms &gamma; = 3 - K0/initial_num .
 */
unweighted_edge_list generate_dms(V final_num, V initial_num, V K0) {
  CHECK(initial_num > K0);
  unweighted_edge_list es;
  vector<V> vs;
  for (int v = 0; v <= initial_num; ++v) {
    for (int u = 0; u < v; ++u) {
      es.emplace_back(u, v);
    }
    for (int i = 0; i < initial_num - K0; ++i) {
      vs.emplace_back(v);
    }
  }

  for (int v = initial_num + 1; v < final_num; ++v) {
    set<V> next;
    std::uniform_int_distribution<size_t> rng(0, vs.size() - 1);
    while (next.size() < (size_t)initial_num) {
      V u = vs[rng(agl::random)];
      next.insert(u);
    }
    for (auto u : next) {
      es.emplace_back(u, v);
      vs.emplace_back(u);
    }
    for (int i = 0; i < initial_num - K0; ++i) {
      vs.emplace_back(v);
    }
  }

  return es;
}

/**
 * Generate a random scale-free network by the Holme-Kim (HK) model.
 * The degree distribution resulting from the HK model is scale free
 * with power-law coefficient &gamma; = 3.
 * The clustering coefficient C can be increased systematically by increasing P, 
 * and finally it becomes C &asymp; 0.5 .
 * \param initial_num is a number of nodes of the initial connected network.
 * \param final_num is a number of finally generated network.
 * \param P is the probability to perform a triangle formation step. 
 */
unweighted_edge_list generate_hk(V final_num, V initial_num, double P) {
  CHECK(initial_num > 2 && final_num > initial_num);
  unweighted_edge_list es;
  vector<vector<V>> adj(final_num);
  for (int v = 0; v <= initial_num; ++v) {
    for (int u = 0; u < v; ++u) {
      es.emplace_back(u, v);
      adj[u].emplace_back(v);
      adj[v].emplace_back(u);
    }
  }

  std::uniform_real_distribution<> p_rng(0.0, 1.0);
  for (int v = initial_num + 1; v < final_num; ++v) {
    set<V> next;
    V u = -1;
    while (next.size() < (size_t)initial_num) {
      if (next.size() > 0 && p_rng(agl::random) < P) {
        std::uniform_int_distribution<size_t> adj_rng(0, adj[u].size() - 1);
        u = adj[u][adj_rng(agl::random)];
      } else {
        std::uniform_int_distribution<size_t> rng(0, es.size() - 1);
        size_t e = rng(agl::random);
        u = rng(agl::random) % 2 ? es[e].first : es[e].second;
      }
      next.insert(u);
    }
    for (auto u : next) {
      es.emplace_back(u, v);
      adj[u].emplace_back(v);
      adj[v].emplace_back(u);
    }
  }

  return es;
}

/**
 * Generate a random small-world network by the Watts-Strogatz (WS) model.
 * \param num_vertices is a number of nodes of the generated network.
 * \param avg_deg is the average degree of the vertices, which must be even.
 * \param P is the probability of reconnecting each edge.
 */
unweighted_edge_list generate_ws(V num_vertices, V avg_deg, double P) {
  CHECK(num_vertices > 1);
  CHECK(0 < avg_deg && avg_deg < num_vertices);
  CHECK(avg_deg % 2 == 0);
  CHECK(0.0 <= P && P <= 1.0);

  unweighted_edge_list out;
  set<pair<V, V> > es;

  for (V i : make_irange(num_vertices)) {
    for (int j = 1; j <= avg_deg / 2; j++) {
      V u = i, v = (i + j) % num_vertices;
      out.emplace_back(u, v);
      es.insert({u, v});
      es.insert({v, u});
    }
  }

  std::uniform_real_distribution<> p_rng(0.0, 1.0);
  std::uniform_int_distribution<V> v_rng(0, num_vertices - 1);

  for (size_t i = 0; i < out.size(); i++) {
    V u = out[i].first, v = out[i].second;
    if (p_rng(agl::random) > P) continue;

    bool coin = p_rng(agl::random) < 0.5;
    if (coin) {
      swap(u, v);
    }
    es.erase({u, v});
    es.erase({v, u});
    do {
      v = v_rng(agl::random);
    } while (u == v || es.find({u, v}) != es.end());

    if (coin) {
      swap(u, v);
    }
    es.insert({u, v});
    es.insert({v, u});

    out[i].first = u;
    out[i].second = v;
  }

  return out;
}

/**
 * Generate a random network by the configuration model with a given degree sequence.
 * This function discards self-loops and multi-edges although the configuration model
 * generates a multi-graph in general.
 * Therefore the original degree sequence may not be reproduced.
 * \param num_vertices is a number of nodes of the generated network.
 * \param deg_seq is a degree sequence. The degree of the i-th node will be deg_seq[i].
 * Note that the sum of the sequence must be even.
 */
unweighted_edge_list generate_config(V num_vertices, const vector<size_t> &deg_seq) {
  CHECK(num_vertices > 0);
  CHECK(deg_seq.size() == (size_t)num_vertices);
  
  size_t deg_sum = 0;
  for (V i : make_irange(num_vertices)) {
    CHECK(deg_seq[i] < (size_t)num_vertices);
    deg_sum += deg_seq[i];
  }

  CHECK(deg_sum % 2 == 0);

  unweighted_edge_list out;
  set<pair<V, V> > es;
  vector<V> half_edges;
  for (V i : make_irange(num_vertices)) {
    for (size_t j : make_irange(deg_seq[i])) {
      half_edges.push_back(i);
    }
  }

  shuffle(half_edges.begin(), half_edges.end(), agl::random);
  

  for (size_t j : make_irange(deg_sum / 2)) {
    V u = half_edges[j * 2], v = half_edges[j * 2 + 1];
    if (u != v && es.find({u, v}) == es.end()) {
      out.emplace_back(u, v);
      es.insert({u, v});
      es.insert({v, u});
    }
  }

  return out;
}

unweighted_edge_list generate_kronecker_(int scale, size_t num_edges, size_t N, const std::vector<std::vector<double>> &matrix) {
  vector<double> prob_sum(N * N + 1);
  for (int i : make_irange(N)) {
    for (int j : make_irange(N)) {
      prob_sum[i * N + j + 1] = prob_sum[i * N + j] + matrix[i][j];
    }
  }
  CHECK_MSG(abs(prob_sum[N * N] - 1.0) < 1e-3, "the sum of the elements in the matrix must be 1.0");

  unweighted_edge_list out;
  set<pair<V, V>> es;

  size_t inserted_edge = 0;
  uniform_real_distribution<double> rng(0.0, 1.0);
  for (size_t i = 0; i < 2 * num_edges && inserted_edge != num_edges; i++) {
    size_t row = 0, col = 0;
    size_t base = 1;
    for (int j : make_irange(scale)) {
      auto it = upper_bound(prob_sum.begin(), prob_sum.end(), rng(agl::random));
      size_t index = distance(prob_sum.begin(), it) - 1;
      row += base * (index / N);
      col += base * (index % N);
      base *= N;
    }
    
    if (row == col || es.find({row, col}) != es.end()) continue;
    out.emplace_back(row, col);
    es.insert({row, col});
    ++inserted_edge;
  }

  return out;
}

/**
 * Generate a (directed) Kronecker graph by the given square probability matrix.
 * The number of vertices in the resulting graph will be (scale)-th power of N,
 * where N is the size of the row (or column) of the matrix.
 * This function avoids self-loops and multi-edges.
 * Note that the function aborts edge insertion after
 * 2 * (expected #edges) trials, which may result in fewer edges
 * depending on the value of the matrix.
 * \param scale is the number of hierarchy.
 * \param matrix is a probability matrix. It must be a square matrix
 * whose elements are between 0.0 and 1.0.
 */
unweighted_edge_list generate_kronecker(int scale, const std::vector<std::vector<double>> &matrix) {
  CHECK(scale > 0);
  CHECK(matrix.size() > 0);
  size_t N = matrix.size();
  double sum = 0;
  vector<vector<double>> normalized_matrix(N, vector<double>(N));
  for (size_t i : make_irange(N)) {
    CHECK(N == matrix[i].size());
    for (size_t j : make_irange(N)) {
      CHECK(0 <= matrix[i][j] && matrix[i][j] <= 1.0);
      sum += matrix[i][j];
    }
  }

  for (size_t i: make_irange(N)) {
    for (size_t j: make_irange(N)) {
      normalized_matrix[i][j] = matrix[i][j] / sum;
    }
  }

  size_t num_edges = pow(sum, scale);
  return generate_kronecker_(scale, num_edges, N, normalized_matrix);
}

/**
 * Generate a (directed) Kronecker graph by the given normalized square matrix and average degree.
 * The number of vertices in the resulting graph will be (scale)-th power of N,
 * where N is the size of the row (or column) of the matrix.
 * This function avoids self-loops and multi-edges.
 * Note that the function aborts edge insertion after
 * 2 * (expected #edges) trials, which may result in fewer edges
 * depending on the value of the matrix.
 * \param scale is the number of hierarchy.
 * \param avg_deg is the average degree of a generated graph.
 * \param matrix is a normalized matrix. It must be a square matrix
 * whose elements are between 0.0 and 1.0, and the sum of the elements is 1.0.
 */
unweighted_edge_list generate_kronecker(int scale, size_t avg_deg, const std::vector<std::vector<double>> &matrix) {
  CHECK(scale > 0);
  CHECK(matrix.size() > 0);
  size_t N = matrix.size();
  double sum = 0;
  vector<vector<double>> normalized_matrix(N, vector<double>(N));
  for (size_t i : make_irange(N)) {
    CHECK(N == matrix[i].size());
    for (size_t j : make_irange(N)) {
      CHECK(0 <= matrix[i][j] && matrix[i][j] <= 1.0);
      sum += matrix[i][j];
    }
  }

  size_t num_edges = avg_deg;
  for (int i : make_irange(scale)) {
    num_edges *= N;
  }

  return generate_kronecker_(scale, num_edges, N, matrix);
}

/**
 * Generate a (u, v)-flower graph by the given parameters.
 * Starting from (u + v)-vertices cycle graph, the number of vertices in the
 * resulting graph will be equal to or larger than the given number. The fractal
 * dimension D of this graph will be D = log(u + v) / log(u) (u <= v)
 * \param required_num is the minimum number of the vertices. The resulting
 * graph the number of vertices in the resulting graph will be equal to or
 * larger than this number.
 * \param u is the smaller parameter of (u, v)-flower.
 * \param v is the larger parameter of (u, v)-flower.
 */
unweighted_edge_list generate_uv_flower(V required_num, V u, V v) {
  assert(u <= v && u >= 1 && u + v >= 3);
  unweighted_edge_list es;
  es.emplace_back(u + v - 1, 0);
  for (int i = 1; i < u + v; ++i) {
    es.emplace_back(i - 1, i);
  }
  V current_vertices = u + v;
  while (current_vertices < required_num) {
    unweighted_edge_list next;
    for (auto e : es) {
      V s = e.first, t = e.second;
      for (int i = 0; i < u; ++i) {
        V left = current_vertices - 1;
        V right = current_vertices;
        if (i == 0) left = s;
        if (i == u - 1) right = t;
        next.emplace_back(left, right);
        if (right != t) current_vertices++;
      }
      for (int i = 0; i < v; ++i) {
        V left = current_vertices - 1;
        V right = current_vertices;
        if (i == 0) left = s;
        if (i == v - 1) right = t;
        next.emplace_back(left, right);
        if (right != t) current_vertices++;
      }
    }
    es.swap(next);
  }

  return es;
}

/**
 * Generate a Song-Havlin-Makse model (SHM-model) graph.
 * Starting from a star tree, the resulting graph will be a tree.
 * The fractal dimension D of this graph will be D = log(2t + 1) / log3
 * \param required_num is the minimum number of the vertices. The resulting
 * graph the number of vertices in the resulting graph will be equal to or
 * larger than this number.
 * \param initial_num is the number of vertices of the star tree of the first
 * generation.
 * \param t decides the fractal dimension of this graph. The fractal dimension D
 * of this graph will be D = log(2t + 1) / log3
 */
unweighted_edge_list generate_shm(V required_num, V initial_num, int t, double P) {
  assert(P >= 0.0 && P <= 1.0);
  std::uniform_real_distribution<> p_rng(0.0, 1.0);
  assert(t >= 2 && initial_num >= 3);
  unweighted_edge_list es;
  vector<vector<V>> adj(initial_num);
  for (int i = 1; i < initial_num; ++i) {
    es.emplace_back(0, i);
    adj[0].push_back(i);
    adj[i].push_back(0);
  }

  while ((V)adj.size() < required_num) {
    V current_num = adj.size();
    V next_num = current_num;
    for (int i = 0; i < current_num; ++i) next_num += adj[i].size() * t;

    vector<vector<V>> next(next_num);
    unweighted_edge_list next_es;
    V new_comer = current_num;
    for (V s = 0; s < current_num; ++s)
      for (int i = 0; i < (int)adj[s].size() * t; ++i) {
        next_es.emplace_back(s, new_comer);
        next[s].push_back(new_comer);
        next[new_comer].push_back(s);
        new_comer++;
      }
    for (auto e : es) {
      V s = e.first;
      V ns = -1;
      for (V n : next[s])
        if (next[n].size() == 1) {
          ns = n;
          break;
        }
      assert(ns >= 0);
      V t = e.second;
      V nt = -1;
      for (V n : next[t])
        if (next[n].size() == 1) {
          nt = n;
          break;
        }
      assert(nt >= 0);
      next_es.emplace_back(ns, nt);
      next[ns].push_back(nt);
      next[nt].push_back(ns);
    }
    if (p_rng(agl::random) <= P)
      for (auto e : es) {
        next_es.emplace_back(e);
        next[e.second].push_back(e.first);
        next[e.first].push_back(e.second);
      }
    es.swap(next_es);
    adj.swap(next);
  }
  return es;
}

unweighted_edge_list generate_random_planar(V num_vertices, size_t num_edges) {
  using namespace agl::geometry2d;

  uniform_real_distribution<double> urd(0.0, 1.0);
  vector<point_type> points(num_vertices);
  for (auto &p : points) p = point_type(urd(random), urd(random));

  vector<tuple<double, V, V>> es_sorted;
  for (V v : make_irange(num_vertices)) {
    for (V u : make_irange(v)) {
      es_sorted.emplace_back(abs(points[u] - points[v]), u, v);
    }
  }
  sort(es_sorted.begin(), es_sorted.end());

  set<pair<V, V>> es;
  for (const auto &t : es_sorted) {
    if (es.size() >= num_edges) break;
    V u = get<1>(t);
    V v = get<2>(t);
    if (u == v) continue;

    bool ins = false;
    for (const auto &e : es) {
      if (u == e.first || u == e.second ||
          v == e.first || v == e.second) continue;
      ins = does_intersect(segment_type(points[u], points[v]),
                            segment_type(points[e.first], points[e.second]));
      if (ins) break;
    }

    if (!ins) es.emplace(u, v);
    cout << es.size() << endl;
  }

  return unweighted_edge_list(es.begin(), es.end());
}

unweighted_edge_list gen_random_planar_(V num_vertices, size_t num_edges) {
  using namespace agl::geometry2d;

  uniform_real_distribution<double> urd(0.0, 1.0);
  vector<point_type> points(num_vertices);
  for (auto &p : points) p = point_type(urd(random), urd(random));

  set<pair<V, V>> es;
  while (es.size() < num_edges) {
    V u = random(num_vertices);
    V v = random(num_vertices);
    if (u == v) continue;

    bool ins = false;
    for (const auto &e : es) {
      if (u == e.first || u == e.second ||
          v == e.first || v == e.second) continue;
      ins = does_intersect(segment_type(points[u], points[v]),
                            segment_type(points[e.first], points[e.second]));
      if (ins) break;
    }

    if (!ins) es.emplace(u, v);
    cout << es.size() << endl;
  }

  return unweighted_edge_list(es.begin(), es.end());
}

unweighted_edge_list generate_random_spanning_tree(V num_vertices) {
  union_find uf(num_vertices);
  unweighted_edge_list es;
  std::uniform_int_distribution<V> rng(0, num_vertices - 1);

  while (es.size() + 1 < (size_t)num_vertices) {
    V u = rng(agl::random), v = rng(agl::random);
    if (uf.is_same(u, v)) continue;
    es.emplace_back(u, v);
    uf.unite(u, v);
  }

  return es;
}

unweighted_edge_list make_undirected(const unweighted_edge_list& es) {
  unweighted_edge_list out(es.size() * 2);
  for (auto i : make_irange(es.size())) {
    V u = es[i].first, v = es[i].second;
    out[i * 2 + 0] = make_pair(u, v);
    out[i * 2 + 1] = make_pair(v, u);
  }
  sort(out.begin(), out.end());
  out.erase(unique(out.begin(), out.end()), out.end());
  return out;
}

template<>
unweighted_edge_list add_random_weight<unweighted_graph>(const unweighted_edge_list &es) {
  return es;
}

template<>
unweighted_edge_list add_unit_weight<unweighted_graph>(const unweighted_edge_list &es) {
  return es;
}
}  // namespace agl
