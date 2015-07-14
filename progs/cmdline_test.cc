

#include <string>
#include <iostream>
#include <algorithm>
#include "base/detail/CmdLine.h"
#include <memory>

using namespace TCLAP;
using namespace std;

int main(int argc, char** argv)
{
    // Wrap everything in a try block.  Do this every time,
    // because exceptions will be thrown for problems.
    try {

    // Define the command line object.
    CmdLine cmd("Command description message", ' ', "0.9");

    // Define a value argument and add it to the command line.
    ValueArg<string> nameArg("n","name","Name to print",true,"homer","string");
    cmd.add( nameArg );

    auto a = cmd.add<ValueArg<string>>("v","vname","VName to print",false,"homer","string");

    // Define a switch and add it to the command line.
    SwitchArg reverseSwitch("r","reverse","Print name backwards", false);
    cmd.add( reverseSwitch );

    // Define a multiarg. This can appear more than one... like ''prog -i 5 -i 10 -i 2''
    MultiArg<int> itest("i", "intTest", "multi int test", false, "int" );
    cmd.add( itest );

    // Parse the args.
    cmd.parse( argc, argv );

    // Get the value parsed by each arg.
    string name = nameArg.getValue();
    bool reverseName = reverseSwitch.getValue();

    // Do what you intend too...
    if ( reverseName )
    {
        reverse(name.begin(),name.end());
        cout << "My name (spelled backwards) is: " << name << endl;
    }
    else
        cout << "My name is: " << name << endl;

    for(auto& i : itest.getValue()) {
        cout << i << " ";
    }
    cout << endl;

    if(a->isSet())
        cout << a->getValue() << endl;


    } catch (ArgException &e)  // catch any exceptions
    { cerr << "error: " << e.error() << " for arg " << e.argId() << endl; }
}
