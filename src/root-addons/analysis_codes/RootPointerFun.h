#ifndef ROOTPOINTERFUN_H
#define ROOTPOINTERFUN_H

#include <vector>
#include <string>

class TestClass {
public:
    TestClass* ptr;  //||

    const std::string name;

    TestClass(const std::string& Name="");

    void Print() const;
};

class TestClassContainer {
public:
    std::vector<TestClass*> data1;
    std::vector<TestClass*> data2;

    TestClassContainer();

    ~TestClassContainer();

    void Print() const;

};

struct PtrTest {
    static void Write();
    static void Read();
};

#endif
