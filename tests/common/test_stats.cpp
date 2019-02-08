#include "../catch.hpp"
#include "common/stats.hpp"
#include <unistd.h>
#include <string>
#include <stdlib.h>
using namespace Ebi;

TEST_CASE("Test becnhmarking and results", "[stats][common]"){
    SECTION("test timer"){

        std::string timer_name("TEST: timer");
        stats.start_timer(timer_name);
        int sleep_time = 1e3; //microseconds

        usleep(sleep_time);
        stats.stop_timer(timer_name);

        CHECK(fabs(stats._results[timer_name+" time (s)"] - sleep_time*1e-6) <1e-3);
    stats.print_results();
    }
    SECTION("test parallel timer"){

        std::string timer_name("TEST: parallel timer");
        stats.start_timer(timer_name);

        int sleep_time = 1e3; //microseconds
#pragma omp parallel for
        for(int i = 0 ; i < 4; i++) // maybe fail if computer  doesn't have 4 cores
            usleep(sleep_time);

        stats.stop_timer(timer_name);
        CHECK(fabs(stats._results[timer_name+" time (s)"] - sleep_time*1e-6) <1e-3);
    stats.print_results();
    CHECK(stats._results.count(string("TEST: parallel timer time (s)")) == 1);
    CHECK(stats._results.count(string("TEST: timer time (s)")) == 1);
    }
}
