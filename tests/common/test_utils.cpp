#include "../catch.hpp"
#include "common/utils.hpp"

TEST_CASE("Test spherical harmonic function", "[sph-harm]"){

    vector<double> u;
    vector<double> v;
    vector<double> data;
    ifstream in;
    in.open("sph-harm-uv.txt");
    double number =0;
    while(in >> number){
        data.push_back(number);
    }
    int num_data_pts = data.size()/2;
    u.reserve(num_data_pts);
    v.reserve(num_data_pts);
    
    for(int i =0; i < num_data_pts; i++){
        u.push_back(data[i]);
    }

    for(int i =0; i < num_data_pts; i++){
        v.push_back(data[num_data_pts + i]);
    }



    ofstream f;
    f.precision(16);
    f.open("sphharm_values.txt");
    ofstream f_leg;
    f_leg.precision(16);
    f_leg.open("associated_legendre.txt");
    assert(v.size() == u.size());
    int num_grid_pts = v.size();
    int m =1;
    int n = 1;

    double factorial = 1.;
    for(int i = 0; i < 2*m; i++){
        factorial *= double(n + m - i);
    }
    cout << factorial << endl;
    //factorial = sqrt(1./factorial);
    cout << num_grid_pts << ", " << num_data_pts << endl;
    for(int i =0; i < num_grid_pts; i++){
        double theta = u[i];
        double phi = v[i];
        double associated_legendre_value = 
            Test::associated_legendre_function(m, n, cos(theta));
        f_leg << associated_legendre_value << " ";
        double coeff = sqrt((2.*n + 1.)/(4.*M_PI)/factorial);
        if(i < 20)
        cout << (coeff * associated_legendre_value * cos(m*phi)) << endl;
        double sphmn = (coeff * Test::associated_legendre_function(m,n,cos(theta)) * cos(m*phi));
        double sphmmn = (coeff * Test::associated_legendre_function(-m,n,cos(theta))* cos(-m*phi));
        //f << (coeff * associated_legendre_value * cos(m*phi)) << " ";
        f << 1./sqrt(2.)*(sphmn + sphmmn) << " ";
    }
    f.close();
    f_leg.close();

}
