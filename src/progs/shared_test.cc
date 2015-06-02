#include <iostream>

#include <memory>
#include <list>
using namespace std;

struct myclass {
    myclass(int x): a(x) { cout << __PRETTY_FUNCTION__ << endl; }
    ~myclass() { cout << __PRETTY_FUNCTION__ << endl; }
    int a;

};

void readit( const shared_ptr<myclass> x) {
    cout << "x=" <<x->a << endl;
}

void readall( const std::list< shared_ptr<myclass> > li ) {
    for(auto& a : li)
        readit(a);

    const shared_ptr<myclass> q = li.front();
    readit(q);
}

class base {
public:
    virtual void action() { cout << "base" << endl; }
    virtual base* Copy() { return new base(); }
    virtual unique_ptr<base> CopyU() { return unique_ptr<base>(new base()); }
};

template <typename T>
class derived: public base {
public:
    T data;
    virtual void action() { cout << __PRETTY_FUNCTION__ << endl; }
    virtual base* Copy() { return new derived(); }
    virtual unique_ptr<base> CopyU() { return unique_ptr<base>(new derived()); }
};

int main() {

    shared_ptr<myclass> a = make_shared<myclass>( 5 );

    cout << a->a << endl;

    auto b=a;

    b->a=4;
    cout << b->a << endl;
    cout << a->a << endl;

    readit(b);

    //std::shared_ptr<const myclass> c = make_shared<const myclass>(2);
    //readit(c); //ERROR

    std::list< shared_ptr<myclass> > l;
    l.push_back(b);

    readall(l);
    cout << "=============" << endl;

    unique_ptr<base> up1 = unique_ptr<base>( new base() );
    up1->action();

    unique_ptr<base> up2 = unique_ptr<base>( new derived<int>() );
    up2->action();

    unique_ptr<base> up1cp = up1->CopyU();
    up1cp->action();

    unique_ptr<base> up2cp = up2->CopyU();
    up2cp->action();



    return 0;
}
