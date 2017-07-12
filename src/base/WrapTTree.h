#pragma once

#include "TTree.h"

#include <vector>
#include <string>
#include <algorithm>
#include <memory>

namespace ant {

/**
 * @brief The WrapTTree struct simplifies TTree handling in physics classes
 *
 * Example usage:
 *
 *     struct treeTest_t : WrapTTree {
 *       ADD_BRANCH_T(bool,   IsSignal)
 *       ADD_BRANCH_T(double, KinFitChi2)
 *       // ... more branches using macro ADD_BRANCH_T
 *       // ... if reading, ADD_BRANCH_OPT_T specifies branches which may not be present
 *     };
 *
 *     treeTest_t treeTest;
 *
 * ... in constructor, setup for "writing" ...
 *
 *     treeTest.CreateBranches(HistFac.makeTTree("test"));
 *
 * ... get/set branches via operator() ...
 *
 *     treeTest.IsSignal() = b_IsSignal;
 *
 * ... fill (or get entry) ...
 *
 *     treeTest.Tree->Fill();
 *
 * In order to read the tree created above,
 * `LinkBranches()` is used (using `ant::WrapTFileInput`):
 *
 *     WrapTFileInput inputfile("/some/path/to/filename");
 *     treeTest_t treeTest;
 *     if(inputfile.GetObject("test", treeTest.Tree))
 *       treeTest.LinkBranches(); // uses already set Tree
 *
 * Note that WrapTTree even supports branches created with "branchname[sizebranch]" via
 * `WrapTTree::ROOTArray<T>`, which wraps it into an conviniently usable `std::vector<T>`
 */
class WrapTTree {
public:
    /**
     * @brief Tree to be used as usual TTree
     */
    TTree* Tree = nullptr;

    /**
     * @brief CreateBranches prepares the instance for filling the TTree
     * @param tree the tree to be filled
     * @param skipOptional if true, do not create ADD_BRANCH_OPT_T branches
     */
    void CreateBranches(TTree* tree, bool skipOptional = false);

    /**
     * @brief LinkBranches prepares the instance for reading the TTree
     * @param tree the tree to read from, or use already set Tree class member
     * @param requireOptional if true, do not ignore ADD_BRANCH_OPT_T branches silently
     */
    void LinkBranches(TTree* tree = nullptr, bool requireOptional = false);
    // to avoid bogus LinkBranchs(nullptr, true) call
    void LinkBranches(bool requireOptional);

    /**
     * @brief Matches checks if the branch names are all available
     * @param tree the tree to check
     * @param exact if false, the TTree may have additional branches
     * @return true if successful
     */
    bool Matches(TTree* tree = nullptr, bool exact = true, bool nowarn = false) const;

    /**
     * @brief CopyFrom copies contents in branches by name
     * @param src the source of the contents to be copied
     * @return true if successful, false on mismatch
     */
    bool CopyFrom(const WrapTTree& src);

    /**
     * @brief operator bool returns true if Tree is not null
     */
    explicit operator bool() const {
        return Tree != nullptr;
    }

    template<typename T>
    struct Branch_t {
        template<typename... Args>
        Branch_t(WrapTTree& wraptree, const std::string& name, bool* optionalIsPresent,
                 Args&&... args) :
            Name(name),
            // can't use unique_ptr because of std::addressof below
            Value(new T(std::forward<Args>(args)...))
        {
            static_assert(std::is_same<T, TClonesArray>::value ? sizeof... (Args) > 0 : true,
                          "TClonesArray cannot be default constructed (provide contained class as string!)");
            if(Name.empty())
                throw Exception("Branch name empty");
            auto& branches = wraptree.branches;
            // check if branch with this name already exists
            if(std::find(branches.begin(), branches.end(), Name) != branches.end())
                throw Exception("Branch with name '"+Name+"' already exists");

            branches.emplace_back(name,
                                  TClass::GetClass(typeid(T)),
                                  TDataType::GetType(typeid(T)),
                                  reinterpret_cast<void**>(std::addressof(Value.Ptr)),
                                  std::is_base_of<ROOTArray_traits, T>::value,
                                  optionalIsPresent);
        }
        ~Branch_t() = default;
        Branch_t(const Branch_t&) = delete;
        Branch_t& operator= (const Branch_t& other) {
            *Value = *(other.Value);
            return *this;
        }
        Branch_t(Branch_t&&) = delete;
        Branch_t& operator= (Branch_t&&) = delete;

        const std::string Name;

        // implicit conversion
        operator T& () { return *Value; }
        operator const T& () const { return *Value; }
        // assignment/move
        T& operator= (const T& v) { *Value = v; return *Value; }
        T& operator= (T&& v) { *Value = std::move(v); return *Value; }
        // if you need to call methods of T, sometimes operator() is handy
        T& operator() () { return *Value; }
        const T& operator() () const { return *Value; }
        // subscript access for more convenient access
        // templated to use SFINAE for typedefs T::reference, T::const_reference
        template<typename U = T>
        typename U::reference operator[](std::size_t n) { return (*Value)[n]; }
        template<typename U = T>
        typename U::const_reference operator[](std::size_t n) const { return (*Value)[n]; }
    private:
        struct Value_t {
            explicit Value_t(T* ptr) : Ptr(ptr) {}
            T& operator* () { return *Ptr; }
            const T& operator* () const { return *Ptr; }
            ~Value_t() {
                delete Ptr;
            }
            // forbid copy/move completely
            Value_t(const Value_t&) = delete;
            Value_t& operator=(const Value_t&) = delete;
            Value_t(Value_t&&);
            Value_t& operator=(Value_t&&) = delete;
            // access by Branch_t ctor
            T* Ptr;
        };
        Value_t Value;
    };

    template<typename T>
    struct Branch_Opt_t : Branch_t<T> {
        template<typename... Args>
        Branch_Opt_t(WrapTTree& wraptree, const std::string& name, Args&&... args) :
            Branch_t<T>(wraptree, name, std::addressof(IsPresent), std::forward<Args>(args)...),
            IsPresent(false)
        {
        }

        bool IsPresent;

        // assignment/move (not inherited from base class)
        T& operator= (const T& v) { Branch_t<T>::operator=(v); return *this; }
        T& operator= (T&& v) { Branch_t<T>::operator=(v); return *this; }
    };


private:
    // this interface is used only internally in WrapTTree
    struct ROOTArray_traits {
        virtual void  ROOTArray_setSize(int n) =0;
        virtual void* ROOTArray_getPtr() =0;
        virtual int   ROOTArray_getInternalSize() =0;
    protected:
        virtual ~ROOTArray_traits() = default;
    };

    // check at compile-time if T is std::array with size>1
    // and arithmetic underlying type
    // size==1 is excluded as std::array is superfluous then
    template<typename T>
    struct is_arith_std_array {
        static constexpr auto value = false;
        static constexpr auto N = 1;

    };
    template<typename T, std::size_t N_>
    struct is_arith_std_array<std::array<T,N_>> {
        static constexpr auto value = std::is_arithmetic<T>::value && N_>1;
        static constexpr auto N = N_;
    };

public:

    template<typename T>
    struct ROOTArray : std::vector<T>, ROOTArray_traits {
        virtual ~ROOTArray() = default;
    protected:
        static_assert(std::is_arithmetic<T>::value || is_arith_std_array<T>::value,
                      "ROOTArray can only be used with arithmetic types or "
                      "std::array of arithmetic type (size>1)");

        friend class WrapTTree;
        virtual void ROOTArray_setSize(int n) override {
            this->resize(n);
        }
        virtual void* ROOTArray_getPtr() override {
            return this->data();
        }
        virtual int ROOTArray_getInternalSize() override {
            return is_arith_std_array<T>::N;
        }
    };

    // provide some typedefs to avoid comma in macro ADD_BRANCH_T
    // Note: Only use those for two-dimensional arrays, N>1 (not 1 or 0)
    template<std::size_t N>
    using ROOTArray_Float = ROOTArray<std::array<float,N>>;
    template<std::size_t N>
    using ROOTArray_Double = ROOTArray<std::array<double,N>>;
    template<std::size_t N>
    using ROOTArray_Int = ROOTArray<std::array<int,N>>;
    template<std::size_t N>
    using ROOTArray_Long = ROOTArray<std::array<long,N>>;


    // some public exceptions
    struct Exception : std::runtime_error {
        // use ctors
        using std::runtime_error::runtime_error;
    };

    struct ROOTArrayException : Exception {
        using Exception::Exception;
    };

protected:
    // force user to inherit from this class (by making ctor protected)
    // use ADD_BRANCH_T to define branches (see comments above as well)
    explicit WrapTTree(const std::string& branchNamePrefix_ = "");
    ~WrapTTree();

private:
    struct ROOT_branchinfo_t {
        const std::string Name;
        TClass * const ROOTClass;
        const EDataType ROOTType;
        explicit ROOT_branchinfo_t(const std::string& name, TClass* rootClass, EDataType rootType) :
            Name(name), ROOTClass(rootClass), ROOTType(rootType)
        {}
        bool operator==(const std::string& other_name) const {
            return Name == other_name;
        }
    };

    struct ROOT_branch_t : ROOT_branchinfo_t {
        void** const ValuePtr;
        const bool IsROOTArray;
        bool* const OptionalIsPresent; // is nullptr if branch non-optional

        ROOT_branch_t(const std::string& name,
                      TClass* rootClass,
                      EDataType rootType,
                      void** valuePtr,
                      bool isROOTArray,
                      bool* optionalIsPresent) :
            ROOT_branchinfo_t(name, rootClass, rootType),
            ValuePtr(valuePtr),
            IsROOTArray(isROOTArray),
            OptionalIsPresent(optionalIsPresent)
        {
            if(ROOTClass==0 && ROOTType == kOther_t && !IsROOTArray)
                throw Exception("Cannot use type of branch "+Name+" as ROOT branch, as its unknown to ROOT");
        }

        ROOT_branch_t(ROOT_branch_t&&) = default;
    };

    const std::string branchNamePrefix;
    std::vector<ROOT_branch_t> branches;

    struct ROOTArrayNotifier_t;
    const std::unique_ptr<ROOTArrayNotifier_t> ROOTArrayNotifier;
    void HandleROOTArray(const std::string& branchname, void** valuePtr);
};

}

// macro to define branches consistently
#define ADD_BRANCH_T(type, name, args...) WrapTTree::Branch_t<type> name{*this, #name, nullptr, args};
#define ADD_BRANCH_OPT_T(type, name, args...) WrapTTree::Branch_Opt_t<type> name{*this, #name, args};
