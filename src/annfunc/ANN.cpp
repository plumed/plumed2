/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2019 of Wei Chen and Andrew Ferguson

   This module is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This module is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "Function.h"
#include "ActionRegister.h"
#include "cassert"

#include <string>
#include <cmath>
#include <iostream>
#include <stdio.h>

using namespace std;

// #define DEBUG
// #define DEBUG_2
// #define DEBUG_3

namespace PLMD {
namespace function {
namespace annfunc {


class ANN : public Function 
{
private:
  int num_layers;
  vector<int> num_nodes;
  vector<string> layer_types;
  vector<vector<double> > coefficients_of_connections;
  vector<vector<double> > values_of_biased_nodes;
  vector<vector<double> > output_of_each_layer;
  vector<vector<double> > input_of_each_layer;
  vector<double** > coeff;  // each coeff of connection is a matrix

public:
  static void registerKeywords( Keywords& keys );
  explicit ANN(const ActionOptions&);
  virtual void calculate();
  void calculate_output_of_each_layer(const vector<double>& input);
  void back_prop(vector<vector<double> >& derivatives_of_each_layer, int index_of_output_component);
};

PLUMED_REGISTER_ACTION(ANN,"ANN")

void ANN::registerKeywords( Keywords& keys ) {
  Function::registerKeywords(keys);
  keys.use("ARG"); keys.use("PERIODIC");
  keys.add("compulsory", "NUM_LAYERS", "3");   // why default value for a compulsory argument?
  keys.add("compulsory", "NUM_NODES", "");
  keys.add("compulsory", "LAYER_TYPES", "");
  keys.add("numbered", "COEFFICIENTS_OF_CONNECTIONS", "");
  keys.add("numbered", "VALUES_OF_BIASED_NODES", "");
}

ANN::ANN(const ActionOptions&ao):
  Action(ao),
  Function(ao)
{
  parse("NUM_LAYERS", num_layers);
  num_nodes = vector<int>(num_layers);
  layer_types = vector<string>(num_layers - 1);
  output_of_each_layer = vector<vector<double> >(num_layers);
  input_of_each_layer = vector<vector<double> >(num_layers);
  coeff = vector<double** >(num_layers - 1);
  parseVector("NUM_NODES", num_nodes);
  parseVector("LAYER_TYPES", layer_types);
  log.printf("layer_types = %s, %s\n", layer_types[0].c_str(), layer_types[1].c_str());
  log.printf("num_nodes = %d, %d, %d\n", num_nodes[0], num_nodes[1], num_nodes[2]);
  vector<double> temp_single_coeff, temp_single_bias;
  for (int ii = 0; ; ii ++) {
    // parse coeff
    if( !parseNumberedVector("COEFFICIENTS_OF_CONNECTIONS", ii, temp_single_coeff) ) {
      temp_single_coeff=coefficients_of_connections[ii-1];
      break;
    }
    coefficients_of_connections.push_back(temp_single_coeff);
    log.printf("size of temp_single_coeff = %d\n", temp_single_coeff.size());
    log.printf("size of coefficients_of_connections = %d\n", coefficients_of_connections.size());
    // parse bias
    if( !parseNumberedVector("VALUES_OF_BIASED_NODES", ii, temp_single_bias) ) {
      temp_single_bias=values_of_biased_nodes[ii-1];
    }
    values_of_biased_nodes.push_back(temp_single_bias);
    log.printf("size of temp_single_bias = %d\n", temp_single_bias.size());
    log.printf("size of values_of_biased_nodes = %d\n", values_of_biased_nodes.size());
  }

  if(getNumberOfArguments() != num_nodes[0]) {
    error("Number of arguments is wrong");
  }

  auto temp_coeff = coefficients_of_connections; 
  for (int ii = 0; ii < num_layers - 1; ii ++) {
      int num_of_rows, num_of_cols; // num of rows/cols for the coeff matrix of this connection
      num_of_rows = num_nodes[ii + 1];
      num_of_cols = num_nodes[ii];
      assert (num_of_rows * num_of_cols == temp_coeff[ii].size()); // check whether the size matches
      // create a 2d array to hold coefficients
      coeff[ii] = new double*[num_of_rows];
      for (int kk = 0; kk < num_of_rows; kk ++) {
          coeff[ii][kk] = new double[num_of_cols];
      }
      for (int jj = 0; jj < temp_coeff[ii].size(); jj ++) {
          coeff[ii][jj / num_of_cols][jj % num_of_cols] = temp_coeff[ii][jj];
      }
  }
  // check coeff
  for (int ii = 0; ii < num_layers - 1; ii ++) {
    log.printf("coeff %d = \n", ii);
    for (int jj = 0; jj < num_nodes[ii + 1]; jj ++) {
      for (int kk = 0; kk < num_nodes[ii]; kk ++) {
        log.printf("%f ", coeff[ii][jj][kk]);
      }
      log.printf("\n");
    }
  }
  // check bias
  for (int ii = 0; ii < num_layers - 1; ii ++) {
    log.printf("bias %d = \n", ii);
    for (int jj = 0; jj < num_nodes[ii + 1]; jj ++) {
      log.printf("%f ", values_of_biased_nodes[ii][jj]);
    }
    log.printf("\n");
  } 
  log.printf("initialization ended\n");
  // create components
  for (int ii = 0; ii < num_nodes[num_layers - 1]; ii ++) {
    string name_of_this_component = to_string(ii);
    addComponentWithDerivatives(name_of_this_component);  
    componentIsNotPeriodic(name_of_this_component);  
  }
  checkRead();
}

void ANN::calculate_output_of_each_layer(const vector<double>& input) {
  // first layer
  output_of_each_layer[0] = input;
  // following layers
  for(int ii = 1; ii < num_nodes.size(); ii ++) {
    output_of_each_layer[ii].resize(num_nodes[ii]);
    input_of_each_layer[ii].resize(num_nodes[ii]);
    // first calculate input
    for (int jj = 0; jj < num_nodes[ii]; jj ++) {
      input_of_each_layer[ii][jj] = values_of_biased_nodes[ii - 1][jj];  // add bias term
      for (int kk = 0; kk < num_nodes[ii - 1]; kk ++) {
        input_of_each_layer[ii][jj] += coeff[ii - 1][jj][kk] * output_of_each_layer[ii - 1][kk];
      }
    }
    // then get output
    if (layer_types[ii - 1] == string("Linear")) {
      for(int jj = 0; jj < num_nodes[ii]; jj ++) {
          output_of_each_layer[ii][jj] = input_of_each_layer[ii][jj];
      }
    }
    else if (layer_types[ii - 1] == string("Tanh")) {
      for(int jj = 0; jj < num_nodes[ii]; jj ++) {
          output_of_each_layer[ii][jj] = tanh(input_of_each_layer[ii][jj]);
      }
    }
    else if (layer_types[ii - 1] == string("Circular")) {
    assert (num_nodes[ii] % 2 == 0);
    for(int jj = 0; jj < num_nodes[ii] / 2; jj ++) {
        double radius = sqrt(input_of_each_layer[ii][2 * jj] * input_of_each_layer[ii][2 * jj] 
                            +input_of_each_layer[ii][2 * jj + 1] * input_of_each_layer[ii][2 * jj + 1]);
        output_of_each_layer[ii][2 * jj] = input_of_each_layer[ii][2 * jj] / radius;
        output_of_each_layer[ii][2 * jj + 1] = input_of_each_layer[ii][2 * jj + 1] / radius;

        }
    }
    else {
        printf("layer type not found!\n\n");
        return;
    }
  }
#ifdef DEBUG_2
    // print out the result for debugging
    printf("output_of_each_layer = \n");
    // for (int ii = num_layers - 1; ii < num_layers; ii ++) {
    for (int ii = 0; ii < num_layers; ii ++) {
        printf("layer[%d]: ", ii);
        if (ii != 0) {
            cout << layer_types[ii - 1] << "\t";    
        }
        else {
            cout << "input \t" ;
        }
        for (int jj = 0; jj < num_nodes[ii]; jj ++) {
            printf("%lf\t", output_of_each_layer[ii][jj]);
        }
        printf("\n");
    }
    printf("\n");
#endif
  return;
}

void ANN::back_prop(vector<vector<double> >& derivatives_of_each_layer, int index_of_output_component) {
    derivatives_of_each_layer = output_of_each_layer;  // the data structure and size should be the same, so I simply deep copy it
    // first calculate derivatives for bottleneck layer
    for (int ii = 0; ii < num_nodes[num_nodes.size() - 1]; ii ++ ) {
      if (ii == index_of_output_component) {
        derivatives_of_each_layer[num_nodes.size() - 1][ii] = 1; 
      }
      else {
        derivatives_of_each_layer[num_nodes.size() - 1][ii] = 0;
      }
    }
    // the use back propagation to calculate derivatives for previous layers
    for (int jj = num_nodes.size() - 2; jj >= 0; jj --) {
        if (layer_types[jj] == string("Circular")) {
            vector<double> temp_derivative_of_input_for_this_layer;
            temp_derivative_of_input_for_this_layer.resize(num_nodes[jj + 1]);
#ifdef DEBUG
            assert (num_nodes[jj + 1] % 2 == 0);
#endif
            // first calculate the derivative of input from derivative of output of this circular layer
            for(int ii = 0; ii < num_nodes[jj + 1] / 2; ii ++) {
                // printf("size of input_of_each_layer[%d] = %d\n",jj,  input_of_each_layer[jj].size());
                double x_p = input_of_each_layer[jj + 1][2 * ii];
                double x_q = input_of_each_layer[jj + 1][2 * ii + 1];
                double radius = sqrt(x_p * x_p + x_q * x_q);
                temp_derivative_of_input_for_this_layer[2 * ii] = x_q / (radius * radius * radius) 
                                                        * (x_q * derivatives_of_each_layer[jj + 1][2 * ii] 
                                                        - x_p * derivatives_of_each_layer[jj + 1][2 * ii + 1]);
                temp_derivative_of_input_for_this_layer[2 * ii + 1] = x_p / (radius * radius * radius) 
                                                        * (x_p * derivatives_of_each_layer[jj + 1][2 * ii + 1] 
                                                        - x_q * derivatives_of_each_layer[jj + 1][2 * ii]);
            }
#ifdef DEBUG
            for (int mm = 0; mm < num_nodes[jj + 1]; mm ++) {
                printf("temp_derivative_of_input_for_this_layer[%d] = %lf\n", mm, temp_derivative_of_input_for_this_layer[mm]);
                printf("derivatives_of_each_layer[%d + 1][%d] = %lf\n", jj, mm, derivatives_of_each_layer[jj + 1][mm]);
            }
#endif
            // the calculate the derivative of output of layer jj, from derivative of input of layer (jj + 1)
            for (int mm = 0; mm < num_nodes[jj]; mm ++) {
                derivatives_of_each_layer[jj][mm] = 0;
                for (int kk = 0; kk < num_nodes[jj + 1]; kk ++) {
                        derivatives_of_each_layer[jj][mm] += coeff[jj][kk][mm] \
                                                            * temp_derivative_of_input_for_this_layer[kk];
#ifdef DEBUG
                        printf("derivatives_of_each_layer[%d][%d] = %lf\n", jj, mm, derivatives_of_each_layer[jj][mm]);
                        printf("coeff[%d][%d][%d] = %lf\n", jj, kk, mm, coeff[jj][kk][mm]);
#endif
                }
            } 
            // FIXME: some problem here      
        }
        else {
            for (int mm = 0; mm < num_nodes[jj]; mm ++) {
                derivatives_of_each_layer[jj][mm] = 0;
                for (int kk = 0; kk < num_nodes[jj + 1]; kk ++) {
                    if (layer_types[jj] == string("Tanh")) {
                        // printf("tanh\n");
                        derivatives_of_each_layer[jj][mm] += derivatives_of_each_layer[jj + 1][kk] \
                                    * coeff[jj][kk][mm] \
                                    * (1 - output_of_each_layer[jj + 1][kk] * output_of_each_layer[jj + 1][kk]);
                    }
                    else if (layer_types[jj] == string("Linear")) {
                        // printf("linear\n");
                        derivatives_of_each_layer[jj][mm] += derivatives_of_each_layer[jj + 1][kk] \
                                    * coeff[jj][kk][mm] \
                                    * 1;
                    }
                    else {
                        printf("layer type not found!\n\n");
                        return;
                    }
                }
            }
        }
    }
#ifdef DEBUG
    // print out the result for debugging
    printf("derivatives_of_each_layer = \n");
    for (int ii = 0; ii < num_layers; ii ++) {
        printf("layer[%d]: ", ii);
        for (int jj = 0; jj < num_nodes[ii]; jj ++) {
            printf("%lf\t", derivatives_of_each_layer[ii][jj]);
        }
        printf("\n");
    }
    printf("\n");
#endif
    return;
}

void ANN::calculate() {

  vector<double> input_layer_data(num_nodes[0]);
  for (int ii = 0; ii < num_nodes[0]; ii ++) {
    input_layer_data[ii] = getArgument(ii);
  }
  
  calculate_output_of_each_layer(input_layer_data);
  vector<vector<double> > derivatives_of_each_layer;
  
  for (int ii = 0; ii < num_nodes[num_layers - 1]; ii ++) {
    back_prop(derivatives_of_each_layer, ii);
    string name_of_this_component = to_string(ii);
    Value* value_new=getPntrToComponent(name_of_this_component);
    value_new -> set(output_of_each_layer[num_layers - 1][ii]);
    for (int jj = 0; jj < num_nodes[0]; jj ++) {
      value_new -> setDerivative(jj, derivatives_of_each_layer[0][jj]);  // TODO: setDerivative or addDerivative?
    }
#ifdef DEBUG_3
    printf("derivatives = ");
    for (int jj = 0; jj < num_nodes[0]; jj ++) {
      printf("%f ", value_new -> getDerivative(jj));
    }
    printf("\n");
#endif
  }

}

}
}
}