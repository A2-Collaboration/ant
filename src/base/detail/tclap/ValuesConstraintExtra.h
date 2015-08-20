#pragma once

#include <string>
#include "base/detail/tclap/Constraint.h"


namespace TCLAP {


/**
 *@brief ValuesConstraint class for templated container type
 */
template<class Cont>
class ValuesConstraintExtra : public Constraint<typename Cont::value_type>
{
protected:
    void build_descr();
public:

    using T = typename Cont::value_type;

    /**
     * Constructor.
     * \param allowed - container of allowed values.
     */
    ValuesConstraintExtra(const Cont& allowed);

    /**
     * Constructor.
     * \param allowed - container of allowed values.
     */
    ValuesConstraintExtra(Cont&& allowed);

    /**
     * Virtual destructor.
     */
    virtual ~ValuesConstraintExtra() {}

    /**
     * Returns a description of the Constraint.
     */
    virtual std::string description() const;

    /**
     * Returns the short ID for the Constraint.
     */
    virtual std::string shortID() const;

    /**
     * The method used to verify that the value parsed from the command
     * line meets the constraint.
     * \param value - The value that will be checked.
     */
    virtual bool check(const T& value) const;

protected:

    /**
     * The list of valid values.
     */
    const Cont _allowed;

    /**
     * The string used to describe the allowed values of this constraint.
     */
    std::string _typeDesc;

};

const std::string& to_string(const std::string& s) { return s; }

template<class Cont>
void ValuesConstraintExtra<Cont>::build_descr()
{
    for (auto i = _allowed.cbegin(); i!=_allowed.end(); ++i)
    {

        if ( i !=_allowed.begin() ) {
            _typeDesc += "|";
        }
        _typeDesc += to_string(*i);
    }
}

template <class Cont>
ValuesConstraintExtra<Cont>::ValuesConstraintExtra(const Cont& allowed)
    : _allowed(allowed),
      _typeDesc("")
{
    build_descr();
}

template<class Cont>
ValuesConstraintExtra<Cont>::ValuesConstraintExtra(Cont&& allowed)
    : _allowed(allowed),
      _typeDesc("")
{
    build_descr();
}

template<class Cont>
bool ValuesConstraintExtra<Cont>::check( const T& val ) const
{
    if ( std::find(_allowed.begin(),_allowed.end(),val) == _allowed.end() )
        return false;
    else
        return true;
}

template<class Cont>
std::string ValuesConstraintExtra<Cont>::shortID() const
{
    return _typeDesc;
}

template<class Cont>
std::string ValuesConstraintExtra<Cont>::description() const
{
    return _typeDesc;
}

}
