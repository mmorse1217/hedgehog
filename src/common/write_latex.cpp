#include "common/write_latex.hpp"
BEGIN_EBI_NAMESPACE



string mat_to_latex_table(DblNumMat data, vector<string> column_names){
    assert(data.n() == int(column_names.size()));
    /*
     * \begin{table}[]
     * \centering
     * \caption{My caption}
     * \label{my-label}
     * \begin{tabular}{lllll}
     *  &  &  &  &  \\
     *  &  &  &  &  \\
     *  &  &  &  &  \\
     *  &  &  &  & 
     * \end{tabular}
     * \end{table}
     */
    ostringstream table_stream;
    // formatting
    table_stream << "\\begin{table}[]\n";
    table_stream << "\\centering\n";
    table_stream << "\\caption{}\n";
    table_stream << "\\label{}\n";
    table_stream << "\\begin{tabular}{lllll}\n";
    
    int num_rows = data.m();
    int num_columns = data.n();

    // set up column names
    for(int i = 0; i < num_columns-1; i++){
        table_stream << column_names[i] << " & ";
    }
    table_stream << column_names[num_columns-1] << " \\\\\n ";

    for(int i = 0; i < num_rows; i++){
        for(int j = 0; j < num_columns-1; j++){
            table_stream << data(i,j) << " & ";
        }
        table_stream << data(i,num_columns-1) << " \\\\\n";
    }

    table_stream << "\\end{tabular}\n";
    table_stream << "\\end{table}\n";

    string table = table_stream.str();
    return table;
}

END_EBI_NAMESPACE
