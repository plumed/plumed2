/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2022-2023 of Luigi Bonati and Enrico Trizio.

The pytorch module is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The pytorch module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#ifdef __PLUMED_HAS_LIBTORCH
// convert LibTorch version to string
//#define STRINGIFY(x) #x
//#define TOSTR(x) STRINGIFY(x)
//#define LIBTORCH_VERSION TO_STR(TORCH_VERSION_MAJOR) "." TO_STR(TORCH_VERSION_MINOR) "." TO_STR(TORCH_VERSION_PATCH)

#include "core/PlumedMain.h"
#include "function/Function.h"
#include "core/ActionRegister.h"

#include <torch/torch.h>
#include <torch/script.h>

#include <fstream>
#include <cmath>

// Note: Freezing a ScriptModule (torch::jit::freeze) works only in >=1.11
// For 1.8 <= versions <=1.10 we need a hack
// (see https://discuss.pytorch.org/t/how-to-check-libtorch-version/77709/4 and also
// https://github.com/pytorch/pytorch/blob/dfbd030854359207cb3040b864614affeace11ce/torch/csrc/jit/api/module.cpp#L479)
// adapted from NequIP https://github.com/mir-group/nequip
#if ( TORCH_VERSION_MAJOR == 2 || TORCH_VERSION_MAJOR == 1 && TORCH_VERSION_MINOR <= 10 )
#define DO_TORCH_FREEZE_HACK
// For the hack, need more headers:
#include <torch/csrc/jit/passes/freeze_module.h>
#include <torch/csrc/jit/passes/frozen_graph_optimizations.h>
#endif

using namespace std;

namespace PLMD {
namespace function {
namespace pytorch {

//+PLUMEDOC FUNCTION PYTORCH_MODEL
/*
Load a PyTorch model compiled with TorchScript.

This can be a function defined in Python or a more complex model, such as a neural network optimized on a set of data. In both cases the derivatives of the outputs with respect to the inputs are computed using the automatic differentiation (autograd) feature of Pytorch.

By default it is assumed that the model is saved as: `model.ptc`, unless otherwise indicated by the `FILE` keyword. The function automatically checks for the number of output dimensions and creates a component for each of them. The outputs are called node-i with i between 0 and N-1 for N outputs.

Note that this function requires an installation of the LibTorch C++ library. Check the instructions on [the module page](module_pytorch.md) for instructions on how to enable the module.

## Examples

Load a model called `torch_model.ptc` that takes as input two dihedral angles and returns two outputs.

```plumed
phi: TORSION ATOMS=5,7,9,15
psi: TORSION ATOMS=7,9,15,17
model: PYTORCH_MODEL ...
  FILE=regtest/pytorch/rt-pytorch_model_2d/torch_model.ptc
  ARG=phi,psi
...
PRINT FILE=COLVAR ARG=model.node-0,model.node-1
```

*/
//+ENDPLUMEDOC


class PytorchModel :
  public Function {
  unsigned _n_in;
  unsigned _n_out;
  torch::jit::script::Module _model;
  torch::Device device = torch::kCPU;

public:
  explicit PytorchModel(const ActionOptions&);
  void calculate();
  static void registerKeywords(Keywords& keys);

  std::vector<float> tensor_to_vector(const torch::Tensor& x);
};

PLUMED_REGISTER_ACTION(PytorchModel,"PYTORCH_MODEL")

void PytorchModel::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys);
  keys.add("optional","FILE","Filename of the PyTorch compiled model");
  keys.addOutputComponent("node", "default", "Model outputs");
}

// Auxiliary function to transform torch tensors in std vectors
std::vector<float> PytorchModel::tensor_to_vector(const torch::Tensor& x) {
  return std::vector<float>(x.data_ptr<float>(), x.data_ptr<float>() + x.numel());
}

PytorchModel::PytorchModel(const ActionOptions&ao):
  Action(ao),
  Function(ao) {
  // print libtorch version
  std::stringstream ss;
  ss << TORCH_VERSION_MAJOR << "." << TORCH_VERSION_MINOR << "." << TORCH_VERSION_PATCH;
  std::string version;
  ss >> version; // extract into the string.
  log.printf(("  LibTorch version: "+version+"\n").data());

  //number of inputs of the model
  _n_in=getNumberOfArguments();

  //parse model name
  std::string fname="model.ptc";
  parse("FILE",fname);

  //deserialize the model from file
  try {
    _model = torch::jit::load(fname, device);
  }

  //if an error is thrown check if the file exists or not
  catch (const c10::Error& e) {
    std::ifstream infile(fname);
    bool exist = infile.good();
    infile.close();
    if (exist) {
      plumed_merror("Cannot load FILE: '"+fname+"'. Please check that it is a Pytorch compiled model (exported with 'torch.jit.trace' or 'torch.jit.script').");
    } else {
      plumed_merror("The FILE: '"+fname+"' does not exist.");
    }
  }
  checkRead();

// Optimize model
  _model.eval();
#ifdef DO_TORCH_FREEZE_HACK
  // Do the hack
  // Copied from the implementation of torch::jit::freeze,
  // except without the broken check
  // See https://github.com/pytorch/pytorch/blob/dfbd030854359207cb3040b864614affeace11ce/torch/csrc/jit/api/module.cpp
  bool optimize_numerics = true;  // the default
  // the {} is preserved_attrs
  auto out_mod = torch::jit::freeze_module(
                   _model, {}
                 );
  // See 1.11 bugfix in https://github.com/pytorch/pytorch/pull/71436
  auto graph = out_mod.get_method("forward").graph();
  OptimizeFrozenGraph(graph, optimize_numerics);
  _model = out_mod;
#else
  // Do it normally
  _model = torch::jit::freeze(_model);
#endif

// Optimize model for inference
#if (TORCH_VERSION_MAJOR == 2 || TORCH_VERSION_MAJOR == 1 && TORCH_VERSION_MINOR >= 10)
  _model = torch::jit::optimize_for_inference(_model);
#endif

  //check the dimension of the output
  log.printf("  Checking output dimension:\n");
  std::vector<float> input_test (_n_in);
  torch::Tensor single_input = torch::tensor(input_test).view({1,_n_in});
  single_input = single_input.to(device);
  std::vector<torch::jit::IValue> inputs;
  inputs.push_back( single_input );
  torch::Tensor output = _model.forward( inputs ).toTensor();
  vector<float> cvs = this->tensor_to_vector (output);
  _n_out=cvs.size();

  //create components
  for(unsigned j=0; j<_n_out; j++) {
    string name_comp = "node-"+std::to_string(j);
    addComponentWithDerivatives( name_comp );
    componentIsNotPeriodic( name_comp );
  }

  //print log
  log.printf("  Number of input: %d \n",_n_in);
  log.printf("  Number of outputs: %d \n",_n_out);
  log.printf("  Bibliography: ");
  log<<plumed.cite("Bonati, Trizio, Rizzi and Parrinello, J. Chem. Phys. 159, 014801 (2023)");
  log<<plumed.cite("Bonati, Rizzi and Parrinello, J. Phys. Chem. Lett. 11, 2998-3004 (2020)");
  log.printf("\n");

}


void PytorchModel::calculate() {

  // retrieve arguments
  vector<float> current_S(_n_in);
  for(unsigned i=0; i<_n_in; i++) {
    current_S[i]=getArgument(i);
  }
  //convert to tensor
  torch::Tensor input_S = torch::tensor(current_S).view({1,_n_in}).to(device);
  input_S.set_requires_grad(true);
  //convert to Ivalue
  std::vector<torch::jit::IValue> inputs;
  inputs.push_back( input_S );
  //calculate output
  torch::Tensor output = _model.forward( inputs ).toTensor();


  for(unsigned j=0; j<_n_out; j++) {
    auto grad_output = torch::ones({1}).expand({1, 1}).to(device);
    auto gradient = torch::autograd::grad({output.slice(/*dim=*/1, /*start=*/j, /*end=*/j+1)},
    {input_S},
    /*grad_outputs=*/ {grad_output},
    /*retain_graph=*/true,
    /*create_graph=*/false)[0]; // the [0] is to get a tensor and not a vector<at::tensor>

    vector<float> der = this->tensor_to_vector ( gradient );
    string name_comp = "node-"+std::to_string(j);
    //set derivatives of component j
    for(unsigned i=0; i<_n_in; i++) {
      setDerivative( getPntrToComponent(name_comp),i, der[i] );
    }
  }

  //set CV values
  vector<float> cvs = this->tensor_to_vector (output);
  for(unsigned j=0; j<_n_out; j++) {
    string name_comp = "node-"+std::to_string(j);
    getPntrToComponent(name_comp)->set(cvs[j]);
  }

}


} //PLMD
} //function
} //pytorch

#endif //PLUMED_HAS_LIBTORCH
