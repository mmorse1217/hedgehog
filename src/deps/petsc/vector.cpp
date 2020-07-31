#include "deps/petsc/petsc.hpp"
#include "deps/petsc/vector.hpp"

namespace Petsc {

    VecType Vector::to_VecType(VectorType t ){
        switch ( t ) {
            case VectorType::seq:
                return VECSEQ;
            case VectorType::mpi:
                return VECMPI;
            case VectorType::standard:
                return VECSTANDARD;
        }
    }

    Vector::VectorType Vector::to_VectorType( VecType t ) {
        if ( std::strcmp( t, VECSEQ ) ) return VectorType::seq;
        if ( std::strcmp( t, VECMPI ) ) return VectorType::mpi;
        if ( std::strcmp( t, VECSTANDARD ) ) return VectorType::standard;
        throw std::out_of_range( std::string( "VectorType not supported |" ) +
                static_cast<const char*>( t ) + "|" );
    }

    Vector::Vector( size_t size,const VectorType t, MPI_Comm comm ) {
        if ( t == VectorType::seq ) comm = PETSC_COMM_SELF;
        VecCreate( comm, &v_ );
        VecSetType( v_, to_VecType( t ) );
        VecSetSizes( v_, PETSC_DECIDE, static_cast<int>( size ) );
    }

    Vector::Vector( const PetscScalar* input, const size_t size, const VectorType t){
        assert( t != VectorType::standard );
        if ( t == VectorType::seq )
            VecCreateSeqWithArray( PETSC_COMM_SELF, 1, static_cast<int>( size ),
                    input, &v_ );
        if ( t == VectorType::mpi )
            VecCreateMPIWithArray( PETSC_COMM_WORLD, 1,
                    static_cast<int>( size ), PETSC_DECIDE,
                    input, &v_ );
    }
    /*************
    //rule of 4.5:
    *************/
    // copy constructor:
    Vector::Vector( const Vector& other ){
        VecDuplicate( other.v_, &v_ );
        VecCopy( other.v_, v_ );
    }

    // move constructor:
    Vector::Vector( Vector&& other ) : v_( other.v_ ) {
        other.v_ = PETSC_NULL;
    }

    // this should only be used in the implementation, but might be used if I
    // forgot a function;
    Vector::Vector( Vec& in, owner o) : v_( in ) {
        switch ( o ) {
            case owner::self:
                deleter = []( Vec v ) { VecDestroy( &v ); };
                break;
            case owner::other:
                deleter = []( Vec ) {};
                break;
            case owner::bv:
                throw std::invalid_argument(
                    "bv can't be made with this constructor" );
        }
    }

    // assignment operator:
    Vector& Vector::operator=( Vector other )
    {
        swap( *this, other );
        return *this;
    }



    void swap( Vector& first, Vector& second ) // nothrow
    {
        using std::swap;

        swap( first.v_, second.v_ );
    }

    /*************
    //modifiers:
     *************/

    // set value:
    void Vector::set_value( const int n, PetscScalar v )
    {
        VecSetValue( v_, n, v, INSERT_VALUES );
    }
    Vector& Vector::set_all( const PetscScalar v )
    {
        VecSet( v_, v );
        return *this;
    }

    // assemble:
    void Vector::assemble()
    {
        VecAssemblyBegin( v_ );
        VecAssemblyEnd( v_ );
    }

    Vector& Vector::conjugate()
    {
        VecConjugate( v_ );
        return *this;
    }
    /*************
    // calculations:
     *************/
    double Vector::norm( NormType nt ) const
    {
        double norm;
        VecNorm( v_, nt, &norm );
        return norm;
    }

    double Vector::normalize()
    {
        double norm;
        VecNormalize( v_, &norm );
        return norm;
    }

    PetscScalar Vector::normalize_sign()
    {
        PetscScalar sign, *xx;
        if ( !this->rank() ) {
            auto r = this->get_ownership_rows();
            VecGetArray( v_, &xx );
            auto x = *xx;
            for ( int i = 1; i < r[1] - r[0]; ++i )
                x = std::abs( x ) > std::abs( xx[i] ) ? x : xx[i];
            sign = x / std::abs( x );
            VecRestoreArray( v_, &xx );
        }
        MPI_Bcast( &sign, 1, MPIU_SCALAR, 0, this->comm() );
        VecScale( v_, 1.0 / sign );
        return sign;
    }

    Vector Vector::operator-( Vector other ) const
    {
        if ( this->rank() == 0 )
            std::cerr << "using operator-: this is inefficient!!!" << std::endl;
        VecAYPX( other.v_, -1.0, this->v_ );
        return other;
    }
    Vector Vector::operator+( Vector other ) const
    {
        if ( this->rank() == 0 )
            std::cerr << "using operator+: this is inefficient!!!" << std::endl;
        VecAYPX( other.v_, 1.0, this->v_ );
        return other;
    }
    Vector Vector::operator*( const PetscScalar& other ) const
    {
        if ( this->rank() == 0 )
            std::cerr << "using operator*: this is inefficient!!!" << std::endl;
        Vector v{*this};
        VecScale( v.v_, other );
        return v;
    }

    Vector& Vector::operator/=( const PetscScalar& other )
    {
        VecScale( v_, 1. / other );
        return *this;
    }
    Vector& Vector::operator*=( const PetscScalar& other )
    {
        VecScale( v_, other );
        return *this;
    }
    Vector& Vector::operator*=( const Vector& other )
    {
        VecPointwiseMult( v_, v_, other.v_ );
        return *this;
    }
    Vector Vector::operator/( const PetscScalar& other ) const
    {
        if ( this->rank() == 0 )
            std::cerr << "using operator-: this is inefficient!!!" << std::endl;
        Vector v{*this};
        VecScale( v.v_, 1. / other );
        return v;
    }
    Vector& Vector::operator+=( const Vector& other )
    {
        VecAXPY( v_, 1, other.v_ );
        return *this;
    }
    Vector& Vector::operator-=( const Vector& other )
    {
        VecAXPY( v_, -1, other.v_ );
        return *this;
    }
    Vector& Vector::axpy( const PetscScalar& alpha, const Vector& x )
    {
        VecAXPY( v_, alpha, x.v_ );
        return *this;
    }
    /*************
    //getters:
     *************/
    int Vector::rank() const
    {
        int r;
        MPI_Comm_rank( comm(), &r );
        return r;
    }
    MPI_Comm Vector::comm() const
    {
        MPI_Comm c;
        PetscObjectGetComm( (PetscObject)v_, &c );
        return c;
    }

    size_t Vector::size() const
    {
        PetscInt i;
        VecGetSize( v_, &i );
        return static_cast<size_t>( i );
    }

    std::array<PetscInt, 2> Vector::get_ownership_rows() const
    {
        std::array<PetscInt, 2> m;
        VecGetOwnershipRange( v_, &m[0], &m[1] );
        return m;
    }

    Vector Vector::duplicate() const
    {
        Vec n;
        VecDuplicate( this->v_, &n );
        return Vector{n};
    }

    void Vector::print() const { VecView( v_, PETSC_VIEWER_STDOUT_WORLD ); }


    void Vector::to_file( const std::string& filename ) const
    {
        PetscViewer view;
        PetscViewerBinaryOpen( this->comm(), filename.c_str(), FILE_MODE_WRITE,
                &view );
        VecView( v_, view );
        PetscViewerDestroy( &view );
    }


    Vector operator*( const PetscScalar& alpha, Vector b )
    {
        if ( b.rank() == 0 )
            std::cerr << "using operator*: this is inefficient!!!" << std::endl;
        VecScale( b.v_, alpha );
        return b;
    }

    Vector conjugate( Vector v )
    {
        VecConjugate( v.v_ );
        return v;
    }

    Vector abs_square( Vector v )
    {
        VecAbs( v.v_ );
        VecPow( v.v_, 2 );
        return v;
    }

    void Vector::get_local(){
        _has_local = true;
        _local_v = hedgehog::get_local_vector(1, size(), v_);
    }
    void Vector::get_local(int m, int n){
        assert(m*n == int(size()));
        _has_local = true;
        _local_v = hedgehog::get_local_vector(m, n, v_);
    }
    void Vector::restore_local(){
                assert(_local_v._v == v_);
        _local_v.restore_local_vector();
    }

    PetscScalar& Vector::operator()(int i, int j){
        assert(_has_local);
        return _local_v(i,j);
        //return _has_local ? _local_v(i,j) : &std::nan("");
    }


}
