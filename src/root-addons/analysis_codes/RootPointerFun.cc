#include "RootPointerFun.h"
#include <iostream>
#include <iomanip>
#include "TFile.h"
#include "TTree.h"

using namespace std;

TestClass::TestClass(const std::string& Name):
    ptr(nullptr),
    name(Name) {}

void TestClass::Print() const
{
    cout << name << " @" << std::hex << this << " ptr => " << std::hex << ptr;
    if(ptr)
        cout << " " << ptr->name;

    cout << endl;
}

TestClassContainer::TestClassContainer()
{

}

TestClassContainer::~TestClassContainer()
{
    for(auto& e:data1)
        delete e;
    for(auto& e:data2)
        delete e;
}

void TestClassContainer::Print() const
{
    cout << "1====" << endl;
    for(const auto& e : data1) {
       e->Print();
    }
    cout << "2====" << endl;
    for(const auto& e : data2) {
       e->Print();
    }
}

void PtrTest::Write()
{
    TFile of("Ptrrest.root","recreate");

    TestClassContainer* c = new TestClassContainer();

    TTree* t = new TTree("tree","");
    t->Branch("c", c);

    c->data1.push_back(new TestClass("a"));
    c->data1.push_back(new TestClass("b"));
    c->data1.push_back(new TestClass("c"));

    c->data2.push_back(new TestClass("a2"));
    c->data2.push_back(new TestClass("b2"));
    c->data2.push_back(new TestClass("c2"));
    c->data2.push_back(new TestClass("d2"));

    c->data2.at(0)->ptr=c->data1.at(0); // cross-vector ref, single
    c->data2.at(1)->ptr=c->data1.at(2); // cross-vector ref, double
    c->data2.at(2)->ptr=c->data1.at(2); // cross-vector ref, double

    c->data2.at(3)->ptr=c->data2.at(0); // in-vector ref

    cout << "Generated structure:" << endl;
    c->Print();

    t->Fill();

    of.Write();
    of.Close();

}

void PtrTest::Read()
{
    TFile infile("Ptrrest.root","read");

    TestClassContainer* c = new TestClassContainer();

    TTree* t = nullptr;
    infile.GetObject("tree", t);

    t->SetBranchAddress("c", &c);

    t->GetEntry(0);

    cout << "Read structure:" << endl;
    c->Print();

    infile.Close();
}
