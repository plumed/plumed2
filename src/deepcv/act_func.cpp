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
05/01/2021 : Add common functions
*/

#include <iostream>
#include <vector>
#include <math.h>
#include <assert.h>

using namespace std;

// Sigmoid
double sigmoid(double z) {
	return 1 / (1 + exp(-z));
}

// Sigmoid derivative
double sigmoid_deri(double z) {
	return sigmoid(z) * (1 - sigmoid(z));
}

// ReLU
double relu(double x) {
	if (x > 0)
		return x;
	else
		return 0;
	}

// ReLU derivative
double relu_deri(double z) {
	if (z > 0)
		return 1;
	else
		return 0;
}

// Leaky ReLU
double leaky_relu(double z) {
	if (z > 0) return z;
	else return 0.01 * z;
}

// Leaky ReLU derivative
double leaky_relu_deri(double z) {
	if (z > 0) return 1;
	else return 0.01;
}

// Softmax
vector<double> softmax(double* x) {
	assert(x > 0);

	int i;
	int size = sizeof(x);
	double m = -INFINITY, sum = 0.0, constant;
	vector<double> y;

	for (i = 0; i < size; ++i) {
		if (m < x[i]) {
			m = x[i];
		}
	}

	for (i = 0; i < size; ++i) {
		sum += exp(x[i] - m);
	}

	constant = m + log(sum);
	for (i = 0; i < size; ++i) {
		y[i] = exp(x[i] - constant);
	}
	return y;
}

