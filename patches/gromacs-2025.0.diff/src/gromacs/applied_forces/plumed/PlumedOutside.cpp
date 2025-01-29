#include "PlumedOutside.h"

/*PLUMED*/
namespace PLMD
{
PlumedOutside::PlumedOutside() {}
PlumedOutside::~PlumedOutside() {}
void PlumedOutside::setPlumed(Plumed* plumed)
{
    plumed_ = plumed;
    active_ = true;
}

bool PlumedOutside::plumedHREX()
{
    return false;
    // not yet integrated
    // return plumedHREX_;
}

bool PlumedOutside::active()
{
    return active_;
}

PlumedOutside& plumedOutside()
{
    static PlumedOutside outside;
    return outside;
} /* */
} // namespace PLMD
/*ENDPLUMED*/