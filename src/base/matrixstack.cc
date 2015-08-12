#include "matrixstack.h"
#include <cmath>
#include "TVector2.h"

using namespace ant;

void matrixstack::clear()
{
    while(!stack.empty()) {
        delete stack.top();
        stack.pop();
    }
}

matrixstack::matrixstack()
{
    Matrix* initial = new Matrix(3,3);
    initial->UnitMatrix();
    stack.push(initial);
}

matrixstack::~matrixstack()
{
    clear();
}

void matrixstack::PushMatrix()
{
    Matrix* new_matrix = new Matrix( matrix() );
    stack.push(new_matrix);
}

void matrixstack::ApplyMatrix(const matrixstack::Matrix &m)
{
    Matrix tmp(matrix());
    matrix().Mult(tmp,m);
}

void matrixstack::PopMatrix()
{
    if(stack.size()>1) {
        delete stack.top();
        stack.pop();
    }
}

void matrixstack::LoadIdentity()
{
    matrix().UnitMatrix();
}

void matrixstack::Rotate(const double phi)
{
    Matrix m(3,3);
    m.UnitMatrix();
    const double s = sin(phi);
    const double c = cos(phi);
    m(0,0) =  c;
    m(0,1) = -s;
    m(1,0) =  s;
    m(1,1) =  c;

    ApplyMatrix( m );

}

void matrixstack::Translate(const matrixstack::Vector &d)
{
    Matrix m(3,3);
    m.UnitMatrix();
    m(0,2) = d(0);
    m(1,2) = d(1);
    ApplyMatrix( m );
}

void matrixstack::Scale(const double x, const double y)
{
    Matrix m(3,3);
    m.UnitMatrix();
    m(0,0) = x;
    m(1,1) = y;
    ApplyMatrix(m);
}

matrixstack::Vector matrixstack::Transform(const matrixstack::Vector &v)
{
    TVectorT<double> a(3);
    a(0) = v.X();
    a(1) = v.Y();
    a(2) = 1;

    a *= matrix();

    Vector b(a);
    b /= a(2);

    return b;
}

matrixstack::Vector::Vector(): TVectorT<double>(2)
{}

matrixstack::Vector::Vector(const double x, const double y): TVectorT<double>(2)
{
    X() = x;
    Y() = y;
}

matrixstack::Vector::Vector(const TVector2 &v): TVectorT<double>(2)
{
    X() = v.X();
    Y() = v.Y();
}

matrixstack::Vector::Vector(const TVectorT<double> &v): TVectorT<double>(2)
{
    X() = v(0);
    Y() = v(1);
}

matrixstack::Vector &matrixstack::Vector::operator/=(double x)
{
    X() /= x;
    Y() /= x;
    return *this;
}

matrixstack::Vector matrixstack::Vector::operator-(const matrixstack::Vector &v)
{
    Vector tmp(*this);
    tmp-=v;
    return tmp;
}

matrixstack::Vector matrixstack::Vector::operator+(const matrixstack::Vector &v)
{
    Vector tmp(*this);
    tmp+=v;
    return tmp;
}

matrixstack::Vector matrixstack::Vector::operator*(const double x)
{
    Vector tmp(*this);
    tmp*=x;
    return tmp;
}

matrixstack::Vector matrixstack::Vector::operator/(const double x)
{
    Vector tmp(*this);
    tmp/=x;
    return tmp;
}
