#ifndef _STATS_HPP_
#define _STATS_HPP_

#define BENCHMARK true
#include "ebi.hpp"
#include "ebiobject.hpp"
#include <stack>

BEGIN_EBI_NAMESPACE
typedef pair<string, double> Timer;

class SystemStats: EbiObject{
    
public:
    string _file_prefix;
    map<string, double> _results;
    stack<Timer> _timer_stack;
    bool append;

        // push a Timer with name timer_name onto the stack with the 
        // appropriate start time
    void start_timer(string timer_name);
        
    // compute the end-time first
        // pop a Timer with name timer_name from the top of the stack,
        // and add into _results the pair (timer_name, end_time - start_time)
        //
        //use omp_get_wtime()
    void stop_timer(string timer_name);

        //  add into _results the pair (result_name, result)
        //  overwrite existing result with the same name
    void add_result(string result_name, double result);
    void result_plus_equals(string result_name, double result);

        // print contents of _results to stdout
    void print_results();
        // store the contents of _results as a table with one row or column with
        // desired results, in the file whose name has the relavent parameters 
        // in the title
    void dump();

    void dump_key_values(string test_name, string filename, string ext=".pydict");
    void clear(){
        _results.clear();
    }

    void store_relative_error(Vec true_potential, Vec computed_potential, 
        NormType norm_type, string result_name);


};

extern SystemStats stats;

END_EBI_NAMESPACE
#endif
