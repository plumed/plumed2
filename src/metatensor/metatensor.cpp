/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Copyright (c) 2024 Guillaume Fraux

    This module is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License as published by the Free
    Software Foundation, either version 3 of the License, or (at your option)
    any later version.
 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#if !defined(__PLUMED_HAS_LIBTORCH) || !defined(__PLUMED_HAS_METATENSOR)

// give a nice error message if the user tries to enable
// metatensor without enabling the corresponding libraries
#error "can not compile the metatensor module without the corresponding libraries, either the disable metatensor module or configure with `--enable-metatensor --enable-libtorch` and make sure the libraries can be found"

#else

#include "core/ActionAtomistic.h"
#include "core/ActionWithValue.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"


#include <torch/script.h>
#include <metatensor/torch.hpp>


namespace PLMD {

class MetatensorPlumedAction: public ActionAtomistic, public ActionWithValue {
public:
    static void registerKeywords(Keywords& keys);
    explicit MetatensorPlumedAction(const ActionOptions&);

    void calculate() override;
    void apply() override;
    unsigned getNumberOfDerivatives() override;

private:

    metatensor_torch::TorchTensorMap output_;
};

PLUMED_REGISTER_ACTION(MetatensorPlumedAction, "METATENSOR")

void MetatensorPlumedAction::registerKeywords(Keywords& keys) {
    Action::registerKeywords(keys);
    ActionAtomistic::registerKeywords(keys);
    ActionWithValue::registerKeywords(keys);

    throw std::runtime_error("unimplemented");
}

MetatensorPlumedAction::MetatensorPlumedAction(const ActionOptions& options):
    Action(options),
    ActionAtomistic(options),
    ActionWithValue(options)
{
    throw std::runtime_error("unimplemented");
}

unsigned MetatensorPlumedAction::getNumberOfDerivatives() {
    // gradients w.r.t. positions (3 x N values) + gradients w.r.t. strain (9 values)
    return 3 * this->getNumberOfAtoms() + 9;
}


void MetatensorPlumedAction::calculate() {
    throw std::runtime_error("unimplemented");
}


void MetatensorPlumedAction::apply() {
    throw std::runtime_error("unimplemented");
}

} // namespace PLMD


#endif
