#include <iostream>
#include "TH1D.h"
#include "Particle.h"
#include <memory>
#include <list>
#include <algorithm>
#include "plot/Histogram.h"
#include "plot/SmartHist.h"
#include "plot/HistogramFactories.h"
#include "TH1D.h"
#include "TFile.h"
#include "TRint.h"

using namespace std;
using namespace ant;


int main(int argc, char** argv) {
    TRint app("A", &argc, argv);

    TFile output("histo-dev.root", "recreate");

    SmartHist1<int> a = SmartHist1<int>::makeHist([] (int a) { return a * 1;},"ta","a","b",BinSettings(10),"t_a");
    a.Fill(5);

    SmartHist1<int> b(SmartHist1<int>::makeHist([] (int a) { return a * 1;},"tb","a","b",BinSettings(10),"t_b"));
    b.Fill(2);



    SmartHist1<double> f = SmartHist1<double>::makeHist("Direct double","x","y",BinSettings(10));

    std::list<SmartHist1<double>> liste;
    liste.emplace_back(SmartHist1<double>::makeHist("Direct double list1","x","y",BinSettings(10)));
    liste.push_back(   SmartHist1<double>::makeHist("Direct double list2","x","y",BinSettings(10)));

    for( auto& h : liste ) {
        h.Fill(1);
    }

    SmartHistFactory sf("factory");

    auto x = sf.makeHist<int>([] (int a) { return 2*a; },"form factory1","","",BinSettings(20),"from factory1");

    SmartHistFactory sf2("another_factory");
    auto x2 = sf2.makeHist<int>([] (int a) { return 3*a; },"form factory2","","",BinSettings(20),"from another factory2");

    TH1D* roothist = HistogramFactory::Default().Make1D("root","x","y");
    roothist->Fill(2);

    SmartHistFactory sf3("sub_factory", sf2);
    auto x4 = sf3.makeHist<int>([] (int a) { return 3*a; },"form sub factory3","","",BinSettings(20),"from another f3");

    auto x5 = sf2.makeHist<int>([] (int a) { return 3*a; },"form factory2 again","","",BinSettings(20),"from another factory2b");


    SmartHist1<int> c = sf2.Copy(b, "c_of_b");

    c.Fill(3);
    c += b;

    // ERROR!
    // SmartHist1<int> d = c;
    //SmartHist1<int> d(c,"d");

    SmartHist1<int> e;
    // will throw exception
    //e.Fill(20);

//  Working on this:
//    SmartHist1Base addition = c + b;

//    sf2.Add(addition, "addition");

    //********************************************************

    output.Write();

    app.Run();
}
