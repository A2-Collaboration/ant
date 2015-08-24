#include "Slowcontrol.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::slowcontrol;



Variable::Variable(Receiver* receiver, const TSlowControl::Key& Key): key(Key)
{
    receiver->RequestSlowcontrol(this);
}


void Receiver::RequestSlowcontrol(Variable* var)
{
    requested_slowcontrols.emplace_back(var);
}


Livetime::Livetime(Receiver* receiver):
    Variable(receiver, TSlowControl::Key(TSlowControl::Type_t::EpicsScaler, "TRIG:TotalLivetime"))
{
}
