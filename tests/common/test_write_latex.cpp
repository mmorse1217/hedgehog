#include "../catch.hpp"
#include "common/write_latex.hpp"
#include <string>
#include <stdlib.h>
using namespace Ebi;

TEST_CASE("Test latex table compiler", "[write_latex][common]"){
    SECTION("test constant data against tablesgenerator.com"){
        DblNumMat data(15, 4);
        setvalue(data, 50.);
        
        vector<string> column_names;
        for (int i = 0; i < data.n(); i++) {
            column_names.push_back(string("Column ")+to_string(i));
        }

        string table = mat_to_latex_table(data, column_names);
    }
    SECTION("test random data against tablesgenerator.com"){
    }
}
