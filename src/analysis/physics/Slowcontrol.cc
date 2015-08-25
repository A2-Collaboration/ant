#include "Slowcontrol.h"

using namespace std;
using namespace ant;
using namespace ant::analysis::slowcontrol;



Variable::Variable(Receiver* receiver)
{
    receiver->RequestSlowcontrol(this);
}


void Receiver::RequestSlowcontrol(Variable* var)
{
    requested_slowcontrols.emplace_back(var);
}



void Distributor::Register(Receiver& rec)
{
    for(auto& var : rec.GetRequestedSlowcontrols()) {

        std::shared_ptr<DataGetter> getter = var->Getter();

        for(const auto& key : getter->GetRequiredKeys() ) {
            requestedKeys.insert(key);
        }

        getters.emplace_back(move(getter));
    }

}

void Distributor::Process(double d)
{
    for(auto& getter : getters) {
        getter->Process(d);
    }
}

