#include "Target.h"

#include "base/root_printable.h"

using namespace ant::analysis::data;
using namespace std;


ostream& Target::Print(ostream& stream) const
{
    stream << "Target"
           << "(Vertex=" << vertex
           << ")";
    return stream;
}
