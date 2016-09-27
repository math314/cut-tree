#include "macros.h"
#include "functions.h"

#include <thread>
#include <future>
#include <chrono>
#include <sys/time.h>
#include <sys/utsname.h>
#include <unistd.h>
using namespace std;

namespace agl {

double get_current_time_sec() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec * 1e-6;
}

}  // namespace agl
