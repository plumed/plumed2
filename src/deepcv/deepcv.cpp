/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
MIT License

Copyright (c) 2021 Rangsiman Ketkaew

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*
Deep Learning for Collective Variables (DeepCV)
https://lubergroup.pages.uzh.ch/deepcv/

Info:
28/11/2020 : Create project
02/12/2021 : Add TensorFlow C++ API
*/

#include "function/Function.h"
#include "core/ActionRegister.h"
#include "cassert"
#include <tensorflow/cc/saved_model/loader.h>
#include <tensorflow/cc/saved_model/tag_constants.h>

#define DEBUG
#define DEBUG_2

namespace PLMD {
namespace function {
namespace deepcv {
// namespace tensorflow {

//+PLUMEDOC DEEPCV_FUNCTION DEEPCV
/*
Load a pretrained DAENN (TensorFlow) model to train a model for real-time learning collective variables on-the-fly.

This module loads external TensorFlow's saved model that is trained using a deep autoencoder neural network (DAENN) \cite Ketkaew2022_deepcv 
in the DeepCV engine \cite Ketkaew2022_daenn. The loaded model will be used to train a model and generate collective variables during enhanced sampling simulation. 

Tensor Session is dynamically created and optimized training parameters including weights and biases are 
called via TensorFlow C++ API (see https://www.tensorflow.org/api_docs/cc).

\par Examples
Define a model that takes as inputs internal coordinates. In this case, we use bond distances, bond angles, 
and torsion angles as input for primary input layer of the DAENN model, and use only bond distances as input 
for secondary input layer. The primary and secondary input leyers are treated separately by primary and 
secondary loss function. respectively.

\plumedfile
# Define distances 
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=2,3
d3: DISTANCE ATOMS=3,4
d4: DISTANCE ATOMS=4,5
d5: DISTANCE ATOMS=5,6

a1: ANGLE ATOMS=1,2,3
a2: ANGLE ATOMS=2,3,4
a3: ANGLE ATOMS=3,4,5
a4: ANGLE ATOMS=4,5,6

t1: TORSION ATOMS=1,2,3,4
t2: TORSION ATOMS=2,3,4,5
t3: TORSION ATOMS=3,4,5,6

# define CV in a separate block
DEEPCV ...
    PRIMARY_INPUT=d1,d2,d3,d4,d5,a1,a2,a3,a4,t1,t2,t3
    SECONDARY_INPUT=d1,d2,d3,d4,d5
    MODEL=output_model/
    NUM_OUTPUT=2
    LABEL=CV
... DEEPCV

PRINT ARG=CV.node-0,CV.node-1 FILE=colvar
\endplumedfile

You can also define signature key, input layer and output layer names, like this

\plumedfile
DEEPCV ...
    PRIMARY_INPUT=d1,d2,d3,d4,d5,a1,a2,a3,a4,t1,t2,t3
    SECONDARY_INPUT=d1,d2,d3,d4,d5
    MODEL=output_model/
    SIGNATURE_DEF=serving_default
    INPUT_TENSOR_1=serving_default_args_0:0
    INPUT_TENSOR_2=serving_default_args_0_1:0
    OUTPUT_TENSOR=StatefulPartitionedCall:0
    NUM_OUTPUT=2
    LABEL=CV
... DEEPCV
\endplumedfile

\par Tips
To figure out the name of input and output tensors in a TensorFlow saved model (SavedModel format), 
you can use TensorFlow's saved_model_cli tool which can be found in Python environments folder
<yourenv>/bin/saved_model_cli. Run the following command to check the model saved by DAENN after training 
in, e.g., output_model/ folder:

\verbatim
$ saved_model_cli show --tag_set serve --signature_def serving_default --all --dir output_model/
\endverbatim

The output should look like this:

\verbatim
...
signature_def['serving_default']:
  The given SavedModel SignatureDef contains the following input(s):
    inputs['args_0'] tensor_info:
        dtype: DT_FLOAT
        shape: (-1, 12)
        name: serving_default_args_0:0
    inputs['args_0_1'] tensor_info:
        dtype: DT_FLOAT
        shape: (-1, 5)
        name: serving_default_args_0_1:0
  The given SavedModel SignatureDef contains the following output(s):
    outputs['dense_2'] tensor_info:
        dtype: DT_FLOAT
        shape: (-1, 2)
        name: StatefulPartitionedCall:0
  Method name is: tensorflow/serving/predict
...
\endverbatim

where in this case serving_default_args_0:0 and serving_default_args_0_1:0 are the names of input tensors 1 
and 2 (corresponding to primary and secondary inputs), respectively, and StatefulPartitionedCall:0 is the names 
of the output tensor (corresponding to the middle hidden layer, aka the output layer of encoder), respectively.

You can also notice that the shape of the first and second layers are (-1, 12) and (-1, 5). It means that 
the model was trained in DeepCV with number of inputs of 12 (d1, d2, d3, d4, d5, a1, a2, a3, a4, t1, t2, t3) 
for the primary loss, and 5 (d1, d2, d3, d4, d5) for the secondary loss functions, respectively.
And for the shape of output layer is (-1, 2), which "2" corresponds to the two neurons in the output layer of 
encoder.

*/
//+ENDPLUMEDOC

class DeepCV: public Function {
    private:
    /// Path to model folder
    std::string model_folder;
    /// Primary layer input
    std::vector<Value*> primary_input;
    std::vector<double> primary_input_vec;
    /// Secondary layer input
    std::vector<Value*> secondary_input;
    std::vector<double> secondary_input_vec;
    /// Signature def key
    /// Default "serving_default" signature name set by tf.saved_model.save method
    std::string signature_def_key;
    /// The name of input tensor
    std::string input_tensor_1_name;
    std::string input_tensor_2_name;
    /// The name of output tensor
    std::string output_tensor_name;
    /// Number of output nodes (CVs)
    int num_output;
    /// Model bundle
    tensorflow::SavedModelBundle model_bundle;
    public:
    /// Create manual
    static void registerKeywords( Keywords& keys );
    /// Constructor
    explicit DeepCV(const ActionOptions&);
    /// Check status of a session
    void checkStatus(tensorflow::Status&);
    /// Get shape of the tensor
    void get_tensor_shape(tensorflow::Tensor&, std::vector<int>&);
    /// Load model --> tensor
    void load_model();
    /// Calculate CVs
    void calculate() override;
};

// Register the name of the module code block in plumed
PLUMED_REGISTER_ACTION(DeepCV,"DEEPCV")

void DeepCV::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys);
  keys.remove("ARG");
  keys.add("compulsory","MODEL","Folder name of pretrained TensorFlow model (SavedModel)");
  keys.add("compulsory","PRIMARY_INPUT","Inputs for primary layer");
  keys.add("compulsory","SECONDARY_INPUT","Inputs for secondary layer");
  keys.add("optional","SIGNATURE_DEF","Name that identifies an entry point of the graph");
  keys.add("optional","INPUT_TENSOR_1","The name of input tensor 1 (primary layer of encoder)");
  keys.add("optional","INPUT_TENSOR_2","The name of input tensor 2 (secondary layer of encoder)");
  keys.add("optional","OUTPUT_TENSOR","The name of output tensor (output layer of encoder");
  keys.add("compulsory","NUM_OUTPUT","2","Number of output nodes corresponding to number of CVs. Default to 2.");
  keys.addOutputComponent("node","default","All \\f$n\\f$ CVs are calculated via two unit nodes. "
                          "The first and second CVs will be labeled <em>label</em>.node-0 and "
                          "<em>label</em>.node-1, respectively.");
}

DeepCV::DeepCV(const ActionOptions&ao):
    Action(ao),
    Function(ao)
    {
    // Parse path to model folder
    parse("MODEL", model_folder);
    // Get primary inputs
    parseArgumentList("PRIMARY_INPUT", primary_input);
    requestArguments(primary_input);
    if (primary_input.size() > 0) log.printf("  More than 1 input will be used for primary loss.\n");
    primary_input_vec.resize(primary_input.size());
    for (unsigned i=0; i<primary_input.size(); ++i) {
        log.printf("  %s : %f\n", primary_input[i]->getName().c_str(), primary_input[i]->get());
        primary_input_vec[i] = primary_input[i]->get();
    }
    // Get Secondary inputs
    parseArgumentList("SECONDARY_INPUT", secondary_input);
    requestArguments(secondary_input);
    if (secondary_input.size() > 0) log.printf("  More than 1 input will be used for secondary loss.\n");
    secondary_input_vec.resize(secondary_input.size());
    for (unsigned i=0; i<secondary_input.size(); ++i) {
        log.printf("  %s : %f\n", secondary_input[i]->getName().c_str(), secondary_input[i]->get());
        secondary_input_vec[i] = secondary_input[i]->get();
    }

    parse("SIGNATURE_DEF", signature_def_key);
    parse("INPUT_TENSOR_1", input_tensor_1_name);
    parse("INPUT_TENSOR_2", input_tensor_2_name);
    parse("OUTPUT_TENSOR", output_tensor_name);
    parse("NUM_OUTPUT", num_output);
    // Check if PLUMED input is loaded properly
    checkRead();
    // Create all components
    log.printf("  Output model folder %s\n", model_folder.c_str());
    if (signature_def_key.empty()) signature_def_key = "serving_default"; 
    if (input_tensor_1_name.empty()); input_tensor_1_name = "serving_default_args_0:0";
    if (input_tensor_2_name.empty()); input_tensor_2_name = "serving_default_args_0_1:0";
    if (output_tensor_name.empty()); output_tensor_name = "StatefulPartitionedCall:0";
    if (num_output > 2) warning("Specified number of outputs is greater than 2. Suggested number of CVs is 2.");
    for (unsigned i=0; i<num_output; ++i) {
    std::string num; Tools::convert(i,num);
    addComponentWithDerivatives("node-"+num);
    componentIsNotPeriodic("node-"+num);
    getPntrToComponent(i)->resizeDerivatives( getNumberOfDerivatives() );
  }
}

void DeepCV::checkStatus(tensorflow::Status& status) {
    if (!status.ok()) {
        // error("Status failed: " + status.ToString());
        error("Status failed >> Failed to load model: " + status.error_message());
    }
    return;
}

void DeepCV::get_tensor_shape(tensorflow::Tensor& tensor, std::vector<int>& shape){
    int num_dimensions = tensor.shape().dims();
    for (unsigned ii=0; ii<num_dimensions; ++ii) {
        shape.push_back(tensor.shape().dim_size(ii));
    }
    return;
}

void DeepCV::load_model() {
    // Create sessions
    tensorflow::SessionOptions session_options = tensorflow::SessionOptions();
    tensorflow::RunOptions run_options = tensorflow::RunOptions();
    tensorflow::Status status = tensorflow::LoadSavedModel(session_options, 
                                                           run_options, 
                                                           model_folder, 
                                                           {tensorflow::kSavedModelTagServe},
                                                           &model_bundle);
    // Check if session is successfully loaded
    checkStatus(status);

#ifdef DEBUG
    // const tensorflow::SignatureDef& signature_def = model_bundle.GetSignatures().at(signature_def_key);
    auto sig_map = model_bundle.GetSignatures();
    auto model_def = sig_map.at(signature_def_key);

    log.printf("Action DEEPCV\n");
    log.printf("  Model Signature:\n");
    for (auto const& p : sig_map) {
        log.printf("    key: %s\n", p.first.c_str());
    }
    log.printf("  Model Input Nodes:\n");
    for (auto const& p : model_def.inputs()) {
        log.printf("    key: %s value: %s\n", p.first.c_str(), p.second.name().c_str());
    }
    log.printf("  Model Output Nodes:\n");
    for (auto const& p : model_def.outputs()) {
        log.printf("    key: %s value: %s\n", p.first.c_str(), p.second.name().c_str());
    }
#endif
    return;
}

void DeepCV::calculate() {
    // Load a model
    load_model();
    // Get and store argument from PLUMED
    // int n_arg = getNumberOfArguments();
    int n_primary_arg = primary_input_vec.size();
    int n_secondary_arg = secondary_input_vec.size();
    //-------------------------------------------------
    // create input tensors to feed data
    //-------------------------------------------------
    // tensorflow::TensorShape(m) for 1D data
    // tensorflow::TensorShape({m,n}) for 2D data
    tensorflow::Tensor input_1(tensorflow::DT_FLOAT, tensorflow::TensorShape({1,n_primary_arg}));
    tensorflow::Tensor input_2(tensorflow::DT_FLOAT, tensorflow::TensorShape({1,n_secondary_arg}));
    // fill in tensor
    auto input_1_mat = input_1.matrix<float>();
    for (unsigned i=0; i<n_primary_arg; ++i) {
        // input_1_mat(0, i) = (double)getArgument(i);
        input_1_mat(0, i) = primary_input_vec[i];
    }
    auto input_2_mat = input_1.matrix<float>();
    for (unsigned i=0; i<n_secondary_arg; ++i) {
        input_2_mat(0, i) = secondary_input_vec[i];
    }
    // there might be more efficient way to copy vector to tensor, e.g., using std::copy or std:;memcpy
    // std::copy( getArgument.begin(), getArgument.end(), input.flat<float>().data() )
    // Make a dict for multiple input tensors
    typedef std::vector<std::pair<std::string, tensorflow::Tensor>> tensor_dict;
    tensor_dict feed_dict = {
        { input_tensor_1_name, input_1},
        { input_tensor_2_name, input_2}
    };
    //-------------------------------------------------
    // Create output placeholder tensors for results
    //-------------------------------------------------
    std::vector<std::string> output_names = {
        output_tensor_name
    };
    std::vector<tensorflow::Tensor> outputs;
    //-------------------------------------------------
    // Running inference
    //-------------------------------------------------
    tensorflow::Status status = model_bundle.session->Run(feed_dict, 
                                                     output_names, 
                                                     {}, 
                                                     &outputs);
    // Check if session is successfully loaded
    checkStatus(status);
    // Get result
    auto out = outputs[0].matrix<float>();
    // Get shape of the output tensor
    std::vector<int> shape;
    get_tensor_shape(outputs[0], shape);
    log.printf("  Shape of output tensor: (%d, %d)\n", shape[0], shape[1]);
    // Store output
    std::vector<double> deepcv_cv;
    for (unsigned i=0; i< shape[0] ; ++i) {
        for (unsigned j=0; j< shape[1]; ++j) {
            deepcv_cv.push_back( (double)out(i,j) );
        }
    }

#ifdef DEBUG_2
    log.printf("  Input 1 : %f\n", input_1.DebugString());
    log.printf("  Input 2 : %f\n", input_2.DebugString());
    log.printf("  Output  : %f\n", outputs[0].DebugString());
    log.printf("  Type of output object: %s\n", typeid(out).name());
#endif

#ifdef DEBUG
    // Print output
    log.printf("  Output components:\n");
    for (unsigned i=0; i< shape[0] ; ++i) {
        for (unsigned j=0; j< shape[1]; ++j) {
            log.printf("    %6.6f\n", out(i,j));
        }
        log.printf("\n");
    }
#endif

    // Get CV
    for (unsigned i=0; i<deepcv_cv.size(); ++i) {
    getPntrToComponent(i)->set( deepcv_cv[i] );
    }

    // Close session to free memory space
    model_bundle.session->Close();
}

// } // namespace tensorflow 
} // namespace deepcv
} // namespace function
} // namespace PLMD
