#include "stats.hpp"
#include <string.h>
#include <ctime>
#include <sstream>
BEGIN_EBI_NAMESPACE



// push a Timer with name timer_name onto the stack with the 
// appropriate start time
void SystemStats::start_timer(string timer_name){
    Timer timer(timer_name, omp_get_wtime());
    _timer_stack.push(timer);
}


// compute the end-time first
// pop a Timer with name timer_name from the top of the stack,
// and add the pair (timer_name, end_time - start_time) into _results 
void SystemStats::stop_timer(string timer_name){
    // compute the end time
    Timer most_recent_timer = _timer_stack.top();
    double execution_time = omp_get_wtime() - most_recent_timer.second;

    // make sure you placed timers in a reasonable fashion (ending the most
    // recently started timer first)
    assert(most_recent_timer.first == timer_name.c_str());

    // remove that timer
    _timer_stack.pop();
    string result_name =most_recent_timer.first+" time (s)";
    _results[result_name] = execution_time;

}

//  add into _results the pair (result_name, result)
//  overwrite existing result with the same name
void SystemStats::add_result(string result_name, double result){
#pragma omp critical 
    {
    _results[result_name] = result;
}
}

void SystemStats::result_plus_equals(string result_name, double result){
#pragma omp critical 
    {
    _results[result_name] += result;
}
    //cerr << result_name << ": " << _results[result_name] << endl;
}
// print contents of _results to stdout
void SystemStats::print_results(){
    for(map<string, double>::iterator it = _results.begin();
            it != _results.end();
            it++){
        cout << it->first << ": " << it->second << endl;
    }

}

void SystemStats::store_relative_error(Vec true_potential, Vec computed_potential, 
        NormType norm_type, string result_name){

    assert(norm_type == NORM_1 ||
            norm_type == NORM_2 ||
            norm_type == NORM_INFINITY);
    double true_minus_computed_norm = 0.;
    double true_norm= 0.;
    
    // Compute \| true_potential - computed_potential\|_inf
    Vec difference;
    VecDuplicate(true_potential, &difference);
    VecCopy(true_potential, difference); //difference = true_potential
    
    const double minus_one = -1.; // TODO move to EbiGlobals this is annoying
    VecAXPY(difference, minus_one, computed_potential); //difference = true_potential - computed_potential

    VecNorm(difference, norm_type, &true_minus_computed_norm);
    VecNorm(true_potential, norm_type, &true_norm);
    double relative_error = true_minus_computed_norm/true_norm;
    string to_store(result_name+" relative error (");
    switch(norm_type){
        case NORM_1:
            to_store += "1-norm)";
            break;
        case NORM_2:
            to_store += "2-norm)";
            break;
        case NORM_INFINITY:
            to_store += "inf-norm)";
            break;
        default:
            cerr << "only use 1-, 2-, or inf-norms" << endl;
            assert(0);
    }
    add_result(to_store, relative_error);
}

// store the contents of _results as a table with one row or column with
// desired results, in the file whose name has the relavent parameters 
// in the title
void SystemStats::dump(){
    //system("mkdir ${MOBO_DIR}/src/ebi/results");
    // Generate the timestamp of the current experiment
    // This will be appended to the filename
    ostringstream s;
    time_t t = time(0);
    struct tm* datetime = localtime(&t);
    s  << datetime->tm_hour <<
            ":" << datetime->tm_min <<
            ":" << datetime->tm_sec  <<
            "_" << datetime->tm_mon + 1 <<
            "-" << datetime->tm_mday <<
            "-" <<datetime->tm_year + 1900;
    string timestamp = s.str();
    s.str("");

    // Get all relevant parameters from options files
    PetscBool err;
    // Geometry related business
    double spacing;
    double refined_spacing;
    int64_t patch_order;
    int64_t refinement_factor;
    int64_t boundary_type;
    char surface[100];

    // Quadrature/PDE related business
    int64_t multipole_order;
    int64_t qbkix_order;
    int64_t kernel;
    int64_t is_bounded;

    // Load geometry business
    PetscOptionsGetReal(NULL, "", "-bis3d_spacing", &spacing, &err);
    ebiAssert(err);
    PetscOptionsGetReal(NULL, "", "-bis3d_rfdspacing", &refined_spacing, &err);
    ebiAssert(err);
    PetscOptionsGetInt(NULL, "", "-bd3d_facemap_patch_order", &patch_order, &err);
    ebiAssert(err);
    PetscOptionsGetInt(NULL, "", "-bd3d_facemap_refinement_factor", &refinement_factor, &err);
    ebiAssert(err);
    PetscOptionsGetInt(NULL, "", "-bdtype", &boundary_type, &err);
    ebiAssert(err);
    PetscOptionsGetString(NULL, "", "-bd3d_filename", surface, 100, &err);
    ebiAssert(err);

    // Load quadrature business
    PetscOptionsGetInt(NULL, "", "-bis3d_np", &multipole_order, &err);
    ebiAssert(err);
    PetscOptionsGetInt(NULL, "", "-near_interpolation_num_samples", &qbkix_order, &err);
    ebiAssert(err);
    PetscOptionsGetInt(NULL, "", "-kt", &kernel, &err);
    ebiAssert(err);
    PetscOptionsGetInt(NULL, "", "-dom", &is_bounded, &err); //0 = bounded, 1 = unbounded
    ebiAssert(err);

    // build row to add to results table
    s << kernel << ", " 
      << surface << ", " 
      << is_bounded << ", " 
      << spacing << ", " 
      << refined_spacing << ", " 
      << patch_order<< ", " 
      << refinement_factor<< ", " 
      << qbkix_order<< ", " 
      << boundary_type<< ", " 
      << multipole_order<< ", " 
      << timestamp << ", "
      << _results["Number of patches"] << ", "
      << _results["Far Evaluation time (s)"] << ", " 
      << _results["QBKIX Near time (s)"] << ", " 
      << _results["QBKIX On time (s)"] << ", " 
      << _results["Solve time (s)"] << ", " 
      << _results["QBKIX Near evaluation relative error (inf-norm)"] << ", " 
      << _results["QBKIX On-surface evaluation relative error (inf-norm)"] << ", " 
      << _results["Smooth quadrature (far evaluation) relative error (inf-norm)"]<< ", " ;
    if(_results.count("boundary_distance_ratio")){
        s << _results["boundary_distance_ratio"] << ", ";
    }
    if(_results.count("interpolation_spacing_ratio")){
        s << _results["interpolation_spacing_ratio"] << ", ";
    }
    s << endl;

    ofstream results;
    ofstream lock;

    // Try to open the lock file
    // Got nervous about concurrent writes, so made each job write out resutls
    // separately.
    /*
    lock.open("../results/results.lock",ios::out);
    while(lock.fail()){
        lock.open("../results/results.lock",ios::out);
    }
    cout << "opened lock" << endl;

    lock << "open";
    */
    char outfile[100];
    PetscBool flg = PETSC_FALSE;
    PetscOptionsGetString(NULL, "", "-results", outfile, 100, &flg);
    ebiAssert(flg);

    string output_file(outfile);
    output_file = "../results/"+output_file;
    //results.open("../results/results.csv", std::ios_base::app);
    if(this->append){ // add a line to the results file
        results.open(output_file.c_str(), std::ios_base::app);
    } else {  // overwrite final contents
        results.open(output_file.c_str(), fstream::out);
    }
    cout << "wrote result" << endl;
    results << s.str();
    results.close();
    //lock.close();
    //cout << "deleting lock" << endl;
    //remove("../results/results.lock");
}

void SystemStats::dump_key_values(string test_name, string filename, string ext){
    //string make_test_dir = string("mkdir ${MOBO_DIR}/output/")+test_name;
    //cout << make_test_dir << endl;
    //system(make_test_dir.c_str());
    ostringstream s;
    time_t t = time(0);
    struct tm* datetime = localtime(&t);
    s  << "_" << datetime->tm_mon + 1 <<
            "-" << datetime->tm_mday <<
            "-" <<datetime->tm_year + 1900 << 
            "_" << datetime->tm_hour <<
            "-" << datetime->tm_min <<
            "-" << datetime->tm_sec;
    string timestamp = s.str();
    //filename += timestamp;
    filename += ext;
    s.str("");
    s << "{ " << endl;
    for(const auto& key_value : _results){
        s << "'" << key_value.first << "' : "  << key_value.second <<", " << endl;

    }
    s << "}";
    ofstream results; 
    results.open(filename.c_str(), fstream::out);
    results << s.str();
    results.close();
}

SystemStats stats;
END_EBI_NAMESPACE
