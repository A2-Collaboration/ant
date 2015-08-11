#pragma once

#include "base/printable.h"

#include "TVectorT.h"
#include "TMatrixT.h"
#include "TVector2.h"

#include <stack>


namespace ant {

class matrixstack {
public:
    struct Vector : TVectorT<double>, printable_traits {

        Vector();
        Vector(const double x, const double y);
        Vector(const TVector2& v);
        Vector(const TVectorT<double> &v);

        virtual ~Vector() {}

        Vector& operator/= ( double x );
        Vector  operator-( const Vector& v);
        Vector  operator+( const Vector& v);
        Vector  operator* (const double x);
        Vector  operator/ (const double x);

        const double& X() const { return (*this)(0); }
        double& X() { return (*this)(0); }

        const double& Y() const { return (*this)(1); }
        double& Y() { return (*this)(1); }

        virtual std::ostream& Print(std::ostream& stream) const override;
    };

    struct Matrix : TMatrixT<double>, printable_traits {
        using TMatrixT<double>::TMatrixT;
        virtual std::ostream& Print(std::ostream& stream) const override;
        virtual ~Matrix() {}
    };

private:
    typedef std::stack<Matrix*> Matrix_Stack;
    Matrix_Stack stack;
    void clear();

public:
    matrixstack();

    virtual ~matrixstack();

    void PushMatrix();
    void ApplyMatrix( const Matrix& m );
    void PopMatrix();
    void LoadIdentity();

    void Rotate( const double phi );
    void Translate( const Vector& d );
    void Scale( const double x, const double y);

    const Matrix& matrix() const { return *(stack.top()); }
    Matrix& matrix()       { return *(stack.top()); }

    Vector Transform(const Vector &v );


};

}
