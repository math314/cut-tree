#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <functional>

namespace agl {
// Cygwin/MinGW では to_string が定義されていないバグがある．
// http://stackoverflow.com/questions/12975341/to-string-is-not-a-member-of-std-says-so-g
// http://gcc.gnu.org/bugzilla/show_bug.cgi?id=52015
// これはみさわさんに教えてもらったテク：
// template<typename ...> static inline int getchar_unlocked(void){ return getchar(); }
template<typename T, typename ...>
std::string to_string(const T& n) {
  std::ostringstream stm;
  stm << n;
  return stm.str();
}

// Cygwin/MinGW では stoi もない...
// http://stackoverflow.com/questions/20145488/cygwin-g-stdstoi-error-stoi-is-not-a-member-of-std
template<typename ...>
int stoi(const std::string& str) {
  return atoi(str.c_str());
}

template<typename RangeType>
auto range_to_vector(RangeType range) -> std::vector<decltype(*range.begin())> {
  decltype(range_to_vector(range)) v;
  for (auto x : range) v.emplace_back(x);
  return v;
}

double get_current_time_sec();

}  // namespace agl
