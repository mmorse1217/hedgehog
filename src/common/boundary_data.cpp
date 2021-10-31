#include "kernel3d.hpp"
#include "vec3t.hpp"
#include "utils.hpp"
BEGIN_EBI_NAMESPACE

void Kernel3d::laplace_dirichlet(DblNumMat source_positions,
        DblNumMat source_strengths,
        DblNumMat target_positions, 
        DblNumMat& boundary_data){
    assert(boundary_data.n() == target_positions.n());
    assert(boundary_data.m() == _tdof);
    
    assert(target_positions.m() == DIM);
    assert(source_strengths.n() == source_positions.n());
    assert(source_positions.m() == DIM);
    assert(source_strengths.m() == _sdof);
    
    int num_targets = target_positions.n();
    int num_sources = source_positions.n();

    for(int ti = 0; ti < num_targets; ti++){
        Point3 x(target_positions.clmdata(ti));
        for(int si = 0; si < num_sources; si++){
            Point3 y(source_positions.clmdata(si));
            Point3 r = x - y;
            double norm_r = r.length();

            // if x is close to y, return 1.
            if(norm_r <= 1e-14)
                norm_r = 1.;

            double coeff = 1./(4*M_PI);
            boundary_data(0, ti) += coeff * 1./norm_r * source_strengths(0,si);
        }
    }

}

void Kernel3d::laplace_neumann(DblNumMat source_positions,
        DblNumMat source_strengths,
        DblNumMat target_positions, 
        DblNumMat target_normals, 
        DblNumMat& boundary_data){
    assert(boundary_data.n() == target_positions.n());
    assert(boundary_data.m() == _tdof);

    assert(target_normals.n() == target_positions.n());
    assert(target_normals.m() == DIM);
    assert(target_positions.m() == DIM);

    assert(source_strengths.n() == source_positions.n());
    assert(source_positions.m() == DIM);
    assert(source_strengths.m() == _sdof);
    
    int num_targets = target_positions.n();
    int num_sources = source_positions.n();

    for(int ti = 0; ti < num_targets; ti++){
        Point3 x(target_positions.clmdata(ti));
        Point3 n_x(target_normals.clmdata(ti));
        for(int si = 0; si < num_sources; si++){
            Point3 y(source_positions.clmdata(si));
            Point3 r = x - y;
            double norm_r = r.length();

            // if x is close to y, return 0.
            if(norm_r <= 1e-14)
                norm_r = 1.;

            double r3 = norm_r*norm_r*norm_r;
            double r_dot_nx = dot(r, n_x);
            double coeff = 1./(4*M_PI);
            boundary_data(0, ti) += coeff * 1./r3 *r_dot_nx * source_strengths(0,si);

        }
    }
}

int Tidx(int i, int j, int k){
    // tensor indexint
    return 9*i + 3*j +k;
}
int Midx(int i, int j){
    // tensor indexint
    return 3*i + j;
}
double delta(int i, int j){
    // kronecker delta
    return i == j ? 1. : 0.;
}
vector<double> T(Point3 r, double mu, double nu){
    vector<double> E(27,0.);
    double R = r.length();
    double R3 = R*R*R;
    double R5 = R3*R*R;
    
    double C = 1./(16.*M_PI*(1-nu));
    double C1 = mu*(4*nu-2.);
    double C2 = 6*mu;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                double ri = r(i);
                double rj = r(j);
                double rk = r(k);

                double term1 = -1./R3*(delta(i,k)*rj + delta(i,j)*rk - delta(j,k)*ri);
                double term2 = 1./R5*ri*rj*rk;

                E[Tidx(i,j,k)] = C*(C1*term1 + C2*term2);
            }
        }
    }
    return E;
}

vector<double> stress_tensor(Point3 r, Point3 g, double mu, double nu){
    vector<double> sigma(9,0.);
    auto TT = T(r,mu,nu);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                sigma[Midx(i,k)] += TT[Tidx(j,i,k)]*g(j); 
            }
        }
    }
    return sigma;
}

void Kernel3d::navier_dirichlet(DblNumMat source_positions,
        DblNumMat source_strengths,
        DblNumMat target_positions, 
        DblNumMat& boundary_data){
    assert(boundary_data.n() == target_positions.n());
    assert(boundary_data.m() == _tdof);

    assert(target_positions.m() == DIM);
    assert(source_strengths.n() == source_positions.n());
    assert(source_positions.m() == DIM);
    assert(source_strengths.m() == _sdof);

    int num_targets = target_positions.n();
    int num_sources = source_positions.n();

    for(int ti = 0; ti < num_targets; ti++){
        Point3 x(target_positions.clmdata(ti));
        for(int si = 0; si < num_sources; si++){
            Point3 y(source_positions.clmdata(si));
            Point3 r = x - y;
            double norm_r = r.length();

            // if x is close to y, return 1.
            if(norm_r >= 1e-14){

                double r3 = norm_r*norm_r*norm_r;

                double mu = coefs(0);
                double nu = coefs(1);
                double coeff = 1./(16.*M_PI*mu*(1.-nu));

                Point3 f(source_strengths.clmdata(si));
                double r_dot_f = dot(r, f);

                boundary_data(0, ti) += coeff *((3.-4.*nu)/norm_r*f(0) + r(0)*r_dot_f/r3);
                boundary_data(1, ti) += coeff *((3.-4.*nu)/norm_r*f(1) + r(1)*r_dot_f/r3);
                boundary_data(2, ti) += coeff *((3.-4.*nu)/norm_r*f(2) + r(2)*r_dot_f/r3);
            }
        }

    }
}
void Kernel3d::navier_neumann(DblNumMat source_positions,
        DblNumMat source_strengths,
        DblNumMat target_positions, 
        DblNumMat target_normals, 
        DblNumMat& boundary_data){
    assert(boundary_data.n() == target_positions.n());
    assert(target_normals.n() == target_positions.n());
    assert(boundary_data.m() == _tdof);

    assert(target_positions.m() == DIM);
    assert(target_normals.m() == DIM);
    assert(source_strengths.n() == source_positions.n());
    assert(source_positions.m() == DIM);
    assert(source_strengths.m() == _sdof);

    setvalue(boundary_data, 0.);
    int num_targets = target_positions.n();
    int num_sources = source_positions.n();

    for(int ti = 0; ti < num_targets; ti++){
        Point3 x(target_positions.clmdata(ti));
        Point3 n_x(target_normals.clmdata(ti));
        for(int si = 0; si < num_sources; si++){
            Point3 y(source_positions.clmdata(si));
            Point3 r = x - y;
            double norm_r = r.length();

            // if x is close to y, return 1.
            if(norm_r >= 1e-14){

                double mu = coefs(0);
                double nu = coefs(1);

                Point3 f(source_strengths.clmdata(si));
                auto sigma = stress_tensor(r, f, mu, nu);
                for(int t = 0; t < _tdof; t++){
                    for(int s = 0; s < _sdof; s++){
                        boundary_data(t, ti) += sigma[Midx(t,s)]*n_x(s);
                    }
                }
            }
        }

    }
}


void Kernel3d::dirichlet_bc_from_singularities(
        DblNumMat source_positions,
        DblNumMat source_strengths,
        DblNumMat target_positions, 
        DblNumMat& boundary_data){
    double temp;
    switch(_equation_type){
        case LAPLACE:
            laplace_dirichlet(source_positions, source_strengths, 
                    target_positions, boundary_data);
            break;
        case MOD_HELMHOLTZ:
            assert(0);
        case STOKES:
            temp = _coefs[1];
            _coefs[1] = .5;
            navier_dirichlet(source_positions, source_strengths, 
                    target_positions, boundary_data);
            _coefs[1] = temp;
            break;
        case NAVIER:
            navier_dirichlet(source_positions, source_strengths, 
                    target_positions, boundary_data);
            break;
        default:
            assert(0);

    }

}

void Kernel3d::neumann_bc_from_singularities(
        DblNumMat source_positions,
        DblNumMat source_strengths,
        DblNumMat target_positions, 
        DblNumMat target_normals, 
        DblNumMat& boundary_data){
    double temp;
    switch(_equation_type){
        case LAPLACE:
            laplace_neumann(source_positions, source_strengths, 
                    target_positions, target_normals, boundary_data);
            break;
        case MOD_HELMHOLTZ:
            assert(0);
        case STOKES:
            temp = _coefs[1];
            _coefs[1] = .5;
            navier_neumann(source_positions, source_strengths, 
                    target_positions, target_normals, boundary_data);
            _coefs[1] = temp;
            break;
        case NAVIER:
            navier_neumann(source_positions, source_strengths, 
                    target_positions, target_normals, boundary_data);
            break;
        default:
            assert(0);
    }
}
END_EBI_NAMESPACE
