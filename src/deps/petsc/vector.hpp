#ifndef __PETSC_VECTOR_HPP__
#define __PETSC_VECTOR_HPP__

#include <memory>
#include <iostream>
#include <array>
#include <mutex>
#include <vector>
#include <cassert>
#include <cstring>
#include "common/nummat.hpp"
using Ebi::DblNumMat;
namespace Petsc
{
    class Vector
    {

        public:
            enum class VectorType : char { seq, mpi, standard };

        private:
            static VecType to_VecType(VectorType t );
            static VectorType to_VectorType( VecType t );
            bool _has_local;
            DblNumMat _local_v;


        public:
            Vector( size_t size,
                    const VectorType t = VectorType::standard,
                    MPI_Comm comm = PETSC_COMM_WORLD );
            Vector( const size_t m,
                    const size_t n,
                    const VectorType t = VectorType::seq,
                    MPI_Comm comm = PETSC_COMM_WORLD ):
                Vector(m*n, t, comm){}
            Vector( const PetscScalar* input,
                    const size_t size,
                    const VectorType t = VectorType::seq );
            Vector( const PetscScalar* input,
                    const size_t m,
                    const size_t n,
                    const VectorType t = VectorType::seq ):
                Vector(input, m*n, t){}



            /*************
            //rule of 4.5:
             *************/
            // copy constructor:
            Vector( const Vector& other );

            // move constructor:
            Vector( Vector&& other );

            // this should only be used in the implementation, but might be used if I
            // forgot a function;
            enum class owner : char { self, other, bv };
            Vector( Vec& in, owner o = owner::self );

            Vector( Vec& in, std::function<void( Vec )> deleter_ )
                : v_( in ), deleter( deleter_ ){}

            // assignment operator:
            Vector& operator=( Vector other );

            // destructor:
            ~Vector() { deleter( v_ ); }

            friend void swap( Vector& first, Vector& second ); // nothrow

            /*************
            //modifiers:
             *************/

            // set value:
            void set_value( const int n, PetscScalar v );

            Vector& set_all( const PetscScalar v );

            // assemble:
            void assemble();

            Vector& conjugate();

            /*************
            // calculations:
             *************/
            double norm( NormType nt = NORM_2 ) const;
            double normalize();
            PetscScalar normalize_sign();

            Vector operator-( Vector other ) const;
            Vector operator+( Vector other ) const;
            Vector operator*( const PetscScalar& other ) const;
            Vector operator/( const PetscScalar& other ) const;
            Vector& operator*=( const PetscScalar& other );
            Vector& operator*=( const Vector& other );
            Vector& operator/=( const PetscScalar& other );
            Vector& axpy(const PetscScalar& alpha, const Vector& x);
            Vector& operator-=( const Vector& other );
            Vector& operator+=( const Vector& other );
            /*************
            //getters:
             *************/
            MPI_Comm comm() const;
            int rank() const;

            size_t size() const;

            std::array<PetscInt, 2> get_ownership_rows() const;

            Vector duplicate() const;
            void print() const;

            void to_file( const std::string& filename ) const;


            // this is public incase people want to use methods that aren't defined
            // yet;
            Vec v_;

            // TODO:  make iterator for vector?

            void get_local();
            void get_local(int m, int n);
            void restore_local();
            PetscScalar& operator()(int i, int j);
            DblNumMat local(){
                return _local_v;
            }


        private:
            // state:
            std::function<void( Vec )> deleter{[]( Vec v ) { VecDestroy( &v ); }};
            friend class Matrix;
    };

    Vector operator*( const PetscScalar& alpha, Vector b );
    Vector conjugate( Vector v );
    Vector abs_square( Vector v );
}

#endif
