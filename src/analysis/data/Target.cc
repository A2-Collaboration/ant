#include "Target.h"

#include "base/root_printable.h"

using namespace ant::analysis::data;
using namespace std;


ostream& Target_t::Print(ostream& stream) const
{
    stream << "Target"
           << "(Vertex=" << Vertex
           << ")";
    return stream;
}
