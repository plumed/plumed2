#ifndef PLMD_OUTSIDE_H
#define PLMD_OUTSIDE_H
#include "external/plumed/PlumedInclude.h"

#include "gromacs/utility/exceptions.h"
namespace PLMD
{
class Plumed;
class PlumedOutside
{
    Plumed* plumed_{ nullptr };
    bool    plumedHREX_{ false };
    bool    active_{ false };

public:
    PlumedOutside();
    ~PlumedOutside();
    void setPlumed(Plumed* plumed);
    template<typename Key, typename T>
    void cmd(Key&& key, T&& val)
    try
    {
        plumed_->cmd(std::forward<Key>(key), std::forward<T>(val));
    }
    catch (const std::exception& ex)
    {
        GMX_THROW(gmx::InternalError(
                std::string("An error occurred while running a command from the PLUMED patch:\n")
                + ex.what()));
    }

    bool plumedHREX();
    bool active();
};

PlumedOutside& plumedOutside();
} // namespace PLMD

#endif // PLMD_OUTSIDE_H