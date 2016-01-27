#ifndef ROOTPOINTERFUN_H
#define ROOTPOINTERFUN_H

#include <vector>
#include <string>

/**
 * @brief A test class that has a pointer to another object
 */
class TestClass {
public:
    TestClass* ptr;

    const std::string name;

    TestClass(const std::string& Name="");

    void Print() const;
};

/**
 * @brief A container of TestClass objects. gets written to a TTree.
 */
class TestClassContainer {
public:
    std::vector<TestClass*> data1;
    std::vector<TestClass*> data2;

    TestClassContainer();

    ~TestClassContainer();

    void Print() const;

};

/**
 * @brief The PtrTest struct - functions to be called from ROOT
 */
struct PtrTest {

    /**
     * @brief Generate a test structure out of TestClass and TestClassContainers and write it to a TTree
     */
    static void Write();

    /**
     * @brief Read back what Write() produced
     */
    static void Read();
};

#endif
