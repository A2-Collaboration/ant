#ifndef HIST_H
#define HIST_H

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

#include "plot/Histogram.h"
#include "root_draw.h"

#include <memory>
#include <string>

#include <iostream>

#include <exception>

using namespace std;
namespace ant {

template <typename T>
class BaseFillFunction {
public:
    virtual ~BaseFillFunction() {}
    virtual void Fill(TH1D& hist, const T& data, const double weight=1.0)=0;
    virtual std::unique_ptr<BaseFillFunction<T>> Copy() const =0;
};


template <typename T, typename FunctionType>
class FillFunction: public BaseFillFunction<T> {
protected:
    FunctionType function;

public:
    FillFunction(FunctionType func): function(func) {}
    virtual ~FillFunction() {}
    virtual void Fill(TH1D& hist, const T& data, const double weight=1.0) {
        hist.Fill(function(data),weight);
    }

    virtual std::unique_ptr<BaseFillFunction<T>> Copy() const {
        return std::unique_ptr<BaseFillFunction<T>>( new FillFunction<T,FunctionType>(function) );
    }
};

template<typename T, typename FunctionType>
std::unique_ptr<BaseFillFunction<T>> makeFunc(FunctionType f) {
    return std::unique_ptr<BaseFillFunction<T>>(new FillFunction<T,FunctionType>(f));
}

class uninitialized_histogram : public std::exception {
public:
    const char *what() const throw() { return "Unitialized SmartHist used!"; }
};

class SmartHist1Base: public ant::root_drawable_traits {
public:
    TH1D* histogram;
    bool cleanup;

    void Normalize(double a=1.0) {
        if(histogram && histogram->GetEntries() !=0 )
            histogram->Scale(a / histogram->GetEntries());
    }

    void Scale(double f) {
        if(histogram)
            histogram->Scale(f);
    }

    void Draw(const string &option) const {
        if(histogram)
            histogram->Draw(option.c_str());
    }

    SmartHist1Base& operator+= (const SmartHist1Base& rhs) {
        if(histogram)
            histogram->Add(rhs.histogram);
        return *this;
    }

    SmartHist1Base& operator-= (const SmartHist1Base& rhs) {
        if(histogram)
            histogram->Add(rhs.histogram, -1.0);
        return *this;
    }

    SmartHist1Base& operator*= (const SmartHist1Base& rhs) {
        histogram->Multiply(rhs.histogram);
        return *this;
    }

    SmartHist1Base& operator/= (const SmartHist1Base& rhs) {
        if(histogram)
            histogram->Divide(rhs.histogram);
        return *this;
    }

    SmartHist1Base& operator*= (const double& rhs) {
        if(histogram)
            histogram->Scale(rhs);
        return *this;
    }

    SmartHist1Base& operator/= (const double& rhs) {
        if(histogram)
            histogram->Scale(1.0/rhs);
        return *this;
    }

    virtual void Fill(const double& data, const double weight=1.0) {
        if(histogram)
            histogram->Fill(data, weight);
    }

    TH1D* GetRootHistogram() { return histogram; }

    virtual ~SmartHist1Base() {
        if(cleanup)
            delete histogram;
    }

//    SmartHist1Base( const SmartHist1Base& other ): histogram(new TH1D(*other.histogram)), cleanup(true) {histogram->SetName("");}
//    SmartHist1Base operator+ (const SmartHist1Base& rhs ) {
//        SmartHist1Base result(*this);
//        result += rhs;
//        return result;
//    }

//    SmartHist1Base& operator= (const SmartHist1Base& rhs ) {
//        if(cleanup)
//            delete histogram;
//        histogram = (TH1D*)(rhs.histogram->Clone(""));
//        cleanup = true;
//        return *this;
//    }

protected:
    SmartHist1Base(TH1D* hist): histogram(hist), cleanup(false) {}

};

template <typename T>
class SmartHist1: public SmartHist1Base {
public:
    using FillFuncPtr = std::unique_ptr< BaseFillFunction<T> >;

public:
    FillFuncPtr fillfunction;

public:
    SmartHist1(TH1D& _histogram, FillFuncPtr _function):
        SmartHist1Base(&_histogram),
        fillfunction(std::move(_function)) {}

    SmartHist1( SmartHist1&& ) = default;

    SmartHist1(): SmartHist1Base(nullptr), fillfunction(makeFunc<T>([] (const T& x) { throw uninitialized_histogram(); return 0;})) {}

    void Fill(const T& data, const double weight=1.0) {
        fillfunction->Fill(*histogram, data, weight);
    }

    SmartHist1& operator+= (const SmartHist1Base& rhs) {
        if(histogram)
            histogram->Add(rhs.histogram);
        return *this;
    }

    SmartHist1& operator-= (const SmartHist1Base& rhs) {
        if(histogram)
            histogram->Add(rhs.histogram, -1.0);
        return *this;
    }

    SmartHist1& operator*= (const SmartHist1Base& rhs) {
        histogram->Multiply(rhs.histogram);
        return *this;
    }

    SmartHist1& operator/= (const SmartHist1Base& rhs) {
        if(histogram)
            histogram->Divide(rhs.histogram);
        return *this;
    }

    SmartHist1& operator*= (const double& rhs) {
        if(histogram)
            histogram->Scale(rhs);
        return *this;
    }

    SmartHist1& operator/= (const double& rhs) {
        if(histogram)
            histogram->Scale(1.0/rhs);
        return *this;
    }

    template<typename FunctionType>
    static SmartHist1<T> makeHist(FunctionType func,
        const std::string& title,
        const std::string& xlabel,
        const std::string& ylabel,
        const BinSettings& bins,
        const std::string& name="",
        HistogramFactory& factory=ant::HistogramFactory::Default())
    {
        TH1D* hist = factory.Make1D(
            title,
            xlabel,
            ylabel,
            bins,
            name
            );
        return move(SmartHist1<T>(*hist, makeFunc<T>(func)));
    }

    static SmartHist1<T> makeHist(
        const std::string& title,
        const std::string& xlabel,
        const std::string& ylabel,
        const BinSettings& bins,
        const std::string& name="",
        HistogramFactory& factory=ant::HistogramFactory::Default())
    {
        TH1D* hist = factory.Make1D(
            title,
            xlabel,
            ylabel,
            bins,
            name
            );
        return move(SmartHist1<T>(*hist, makeFunc<T>([] (const T& data) { return data;})));
    }

    SmartHist1& operator= (SmartHist1&& rhs) {
        histogram = rhs.histogram;
        fillfunction = std::move(rhs.fillfunction);
        return *this;
    }
};


// specialization for strings
template<>
SmartHist1<std::string> SmartHist1<std::string>::makeHist(
    const std::string& title,
    const std::string& xlabel,
    const std::string& ylabel,
    const BinSettings& bins,
    const std::string& name,
    HistogramFactory& factory);

}

#endif
