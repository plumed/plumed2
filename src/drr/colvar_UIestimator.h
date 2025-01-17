/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Copyright (c) 2017 of Haochuan Chen (excluding colvar_UIestimator.h)
    Copyright (c) 2017 of Haohao Fu (colvar_UIestimator.h)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_drr_colvar_UIestimator_h
#define __PLUMED_drr_colvar_UIestimator_h
// The original code(https://github.com/fhh2626/colvars/blob/master/src/colvar_UIestimator.h) has been modified by Haochuan Chen.
// Modifications:
// 1. Disable colvars related code.
// 2. Change boltzmann constant.
// 3. Change output precision.
// I(Haochuan Chen) don't know how to maintain this code and how it runs. If you are interested in it, please contact Haohao Fu.

#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
#include <limits>

#include <typeinfo>

// only for colvar module!
// when integrated into other code, just remove this line and "...cvm::backup_file(...)"
// #include "colvarmodule.h"

namespace PLMD {
namespace drr {

namespace UIestimator
{
const int Y_SIZE = 21;
const int HALF_Y_SIZE = 10;
const double BOLTZMANN = 0.0083144621;
const int EXTENDED_X_SIZE = HALF_Y_SIZE;

class n_matrix    // spare matrix, stores the distribution matrix of n(x,y)
{
public:
  n_matrix() {}
  n_matrix(const std::vector<double> & lowerboundary_p,   // lowerboundary of x
           const std::vector<double> & upperboundary_p,   // upperboundary of
           const std::vector<double> & width_p,           // width of x
           const int y_size)           // size of y, for example, ysize=7, then when x=1, the distribution of y in [-2,4] is considered
    : lowerboundary(lowerboundary_p), upperboundary(upperboundary_p), width(width_p)
  {
    this->dimension = lowerboundary.size();
    this->y_size = y_size;     // keep in mind the internal (spare) matrix is stored in diagonal form
    this->y_total_size = int(pow(y_size, dimension) + 0.000001);

    // the range of the matrix is [lowerboundary, upperboundary]
    x_total_size = 1;
    for (int i = 0; i < dimension; i++)
    {
      x_size.push_back(int((upperboundary[i] - lowerboundary[i]) / width[i] + 0.000001));
      x_total_size *= x_size[i];
    }

    // initialize the internal matrix
    matrix.reserve(x_total_size);
    for (int i = 0; i < x_total_size; i++)
    {
      matrix.push_back(std::vector<int>(y_total_size, 0));
    }

    temp.resize(dimension);
  }

  int inline get_value(const std::vector<double> & x, const std::vector<double> & y)
  {
    //if (matrix[convert_x(x)][convert_y(x, y)]!=0)
    //{
    //std::cout<<convert_x(x)<<" "<<convert_y(x, y)<<" "<<x[0]<<" "<<x[1]<<" "<<y[0]<<" "<<y[1]<<" ";
    //std::cout<<matrix[convert_x(x)][convert_y(x, y)]<<"sadasfdasaaaaaaaa"<<std::endl;
    //}
    return matrix[convert_x(x)][convert_y(x, y)];
  }

  void inline set_value(const std::vector<double> & x, const std::vector<double> & y, const int value)
  {
    matrix[convert_x(x)][convert_y(x,y)] = value;
  }

  void inline increase_value(const std::vector<double> & x, const std::vector<double> & y, const int value)
  {
    matrix[convert_x(x)][convert_y(x,y)] += value;
  }

private:
  std::vector<double> lowerboundary;
  std::vector<double> upperboundary;
  std::vector<double> width;
  int dimension;
  std::vector<int> x_size;       // the size of x in each dimension
  int x_total_size;              // the size of x of the internal matrix
  int y_size;                    // the size of y in each dimension
  int y_total_size;              // the size of y of the internal matrix

  std::vector<std::vector<int> > matrix;  // the internal matrix

  std::vector<int> temp;         // this vector is used in convert_x and convert_y to save computational resource

  int convert_x(const std::vector<double> & x)        // convert real x value to its interal index
  {
    for (int i = 0; i < dimension; i++)
    {
      temp[i] = int((x[i] - lowerboundary[i]) / width[i] + 0.000001);
    }

    int index = 0;
    for (int i = 0; i < dimension; i++)
    {
      if (i + 1 < dimension)
      {
        int x_temp = 1;
        for (int j = i + 1; j < dimension; j++)
          x_temp *= x_size[j];
        index += temp[i] * x_temp;
      }
      else
        index += temp[i];
    }
    return index;
  }

  int convert_y(const std::vector<double> & x, const std::vector<double> & y)        // convert real y value to its interal index
  {
    for (int i = 0; i < dimension; i++)
    {
      temp[i] = round((round(y[i] / width[i] + 0.000001) - round(x[i] / width[i] + 0.000001)) + (y_size - 1) / 2 + 0.000001);
    }

    int index = 0;
    for (int i = 0; i < dimension; i++)
    {
      if (i + 1 < dimension)
        index += temp[i] * int(pow(y_size, dimension - i - 1) + 0.000001);
      else
        index += temp[i];
    }
    return index;
  }

  double round(double r)
  {
    return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
  }
};

// vector, store the sum_x, sum_x_square, count_y
template <typename T>
class n_vector
{
public:
  n_vector() {}
  n_vector(const std::vector<double> & lowerboundary,   // lowerboundary of x
           const std::vector<double> & upperboundary,   // upperboundary of
           const std::vector<double> & width_p,                // width of x
           const int y_size,           // size of y, for example, ysize=7, then when x=1, the distribution of y in [-2,4] is considered
           const T & default_value)          //   the default value of T
    :width(width_p)
  {
    this->dimension = lowerboundary.size();

    x_total_size = 1;
    for (int i = 0; i < dimension; i++)
    {
      this->lowerboundary.push_back(lowerboundary[i] - (y_size - 1) / 2 * width[i] - 0.000001);
      this->upperboundary.push_back(upperboundary[i] + (y_size - 1) / 2 * width[i] + 0.000001);

      x_size.push_back(int((this->upperboundary[i] - this->lowerboundary[i]) / this->width[i] + 0.000001));
      x_total_size *= x_size[i];
    }

    // initialize the internal vector
    vector.resize(x_total_size, default_value);

    temp.resize(dimension);
  }

  T inline get_value(const std::vector<double> & x)
  {
    return vector[convert_x(x)];
  }

  void inline set_value(const std::vector<double> & x, const T & value)
  {
    vector[convert_x(x)] = value;
  }

  void inline increase_value(const std::vector<double> & x, const T & value)
  {
    vector[convert_x(x)] += value;
  }
private:
  std::vector<double> lowerboundary;
  std::vector<double> upperboundary;
  std::vector<double> width;
  int dimension;
  std::vector<int> x_size;       // the size of x in each dimension
  int x_total_size;              // the size of x of the internal matrix

  std::vector<T> vector;  // the internal vector

  std::vector<int> temp;         // this vector is used in convert_x and convert_y to save computational resource

  int convert_x(const std::vector<double> & x)        // convert real x value to its interal index
  {
    for (int i = 0; i < dimension; i++)
    {
      temp[i] = int((x[i] - lowerboundary[i]) / width[i] + 0.000001);
    }

    int index = 0;
    for (int i = 0; i < dimension; i++)
    {
      if (i + 1 < dimension)
      {
        int x_temp = 1;
        for (int j = i + 1; j < dimension; j++)
          x_temp *= x_size[j];
        index += temp[i] * x_temp;
      }
      else
        index += temp[i];
    }
    return index;
  }
};

class UIestimator      // the implemension of UI estimator
{
public:
  UIestimator() {}

  //called when (re)start an eabf simulation
  UIestimator(const std::vector<double>& lowerboundary_p,
              const std::vector<double>& upperboundary_p,
              const std::vector<double>& width_p,
              const std::vector<double>& krestr_p,                // force constant in eABF
              const std::string& output_filename_p,              // the prefix of output files
              const int output_freq_p,
              const bool restart_p,                              // whether restart from a .count and a .grad file
              const std::vector<std::string>& input_filename_p,   // the prefixes of input files
              const double temperature_p)
    : lowerboundary(lowerboundary_p), upperboundary(upperboundary_p),
      width(width_p), krestr(krestr_p),
      output_filename(output_filename_p), output_freq(output_freq_p),
      restart(restart_p), input_filename(input_filename_p),
      temperature(temperature_p)
  {

    dimension = lowerboundary.size();

    for (int i = 0; i < dimension; i++)
    {
      sum_x.push_back(n_vector<double>(lowerboundary, upperboundary, width, Y_SIZE, 0.0));
      sum_x_square.push_back(n_vector<double>(lowerboundary, upperboundary, width, Y_SIZE, 0.0));

      x_av.push_back(n_vector<double>(lowerboundary, upperboundary, width, Y_SIZE, 0.0));
      sigma_square.push_back(n_vector<double>(lowerboundary, upperboundary, width, Y_SIZE, 0.0));
    }

    count_y = n_vector<int>(lowerboundary, upperboundary, width, Y_SIZE, 0);
    distribution_x_y = n_matrix(lowerboundary, upperboundary, width, Y_SIZE);

    grad = n_vector<std::vector<double> >(lowerboundary, upperboundary, width, 1, std::vector<double>(dimension, 0.0));
    count = n_vector<int>(lowerboundary, upperboundary, width, 1, 0);

    written = false;
    written_1D = false;

    if (dimension == 1)
    {
      std::vector<double> upperboundary_temp = upperboundary;
      upperboundary_temp[0] = upperboundary[0] + width[0];
      oneD_pmf = n_vector<double>(lowerboundary, upperboundary_temp, width, 1, 0.0);
    }

    if (restart == true)
    {
      input_grad = n_vector<std::vector<double> >(lowerboundary, upperboundary, width, 1, std::vector<double>(dimension, 0.0));
      input_count = n_vector<int>(lowerboundary, upperboundary, width, 1, 0);

      // initialize input_Grad and input_count
      std::vector<double> loop_flag(dimension, 0);
      for (int i = 0; i < dimension; i++)
      {
        loop_flag[i] = lowerboundary[i];
      }
      while (true)
      {
        for (int i = 0; i < dimension; i++)
        {
          input_grad.set_value(loop_flag, std::vector<double>(dimension,0));
        }
        input_count.set_value(loop_flag, 0);

        // iterate over any dimensions
        int i = dimension - 1;
        while (true)
        {
          loop_flag[i] += width[i];
          if (loop_flag[i] > upperboundary[i] - width[i] + 0.00001)
          {
            loop_flag[i] = lowerboundary[i];
            i--;
            if (i < 0)
              goto INITIAL_LOOPEND;
          }
          else
            break;
        }
      }
INITIAL_LOOPEND:
      read_inputfiles(input_filename);
    }
  }

  ~UIestimator() {}

  // called from MD engine every step
  bool update(const int step, std::vector<double> x, std::vector<double> y)
  {

    //std::cout<<"weeeee: :"<<std::endl;
    //for (int i = 0; i < dimension; i++)
    //{
    //	std::cout<<x[i]<<" "<<y[i]<<" ";
    //}
    //std::cout<<std::endl;

    if (step % output_freq == 0)
    {
      calc_pmf();
      write_files();
      //write_interal_data();
    }

    for (int i = 0; i < dimension; i++)
    {
      // for dihedral RC, it is possible that x = 179 and y = -179, should correct it
      // may have problem, need to fix
      if (x[i] > 150 && y[i] < -150)
      {
        //std::vector<double> x_temp(x);
        //x_temp[i] -= 360;
        //update(7, x_temp, y);
        y[i] += 360;
      }
      if (x[i] < -150 && y[i] > 150)
      {
        //std::vector<double> x_temp(x);
        //x_temp[i] += 360;
        //update(7, x_temp, y);
        y[i] -= 360;
      }

      if (x[i] < lowerboundary[i] - EXTENDED_X_SIZE * width[i] + 0.00001 || x[i] > upperboundary[i] + EXTENDED_X_SIZE * width[i] - 0.00001 \
          || y[i] - x[i] < -HALF_Y_SIZE * width[i] + 0.00001 || y[i] - x[i] > HALF_Y_SIZE * width[i] - 0.00001 \
          || y[i] - lowerboundary[i] < -HALF_Y_SIZE * width[i] + 0.00001 || y[i] - upperboundary[i] > HALF_Y_SIZE * width[i] - 0.00001)
        return false;
    }

    //for (int i = 0; i < dimension; i++)
    //{
    //	std::cout<<x[i]<<" "<<y[i]<<" ";
    //}
    //std::cout<<std::endl;

    for (int i = 0; i < dimension; i++)
    {
      sum_x[i].increase_value(y, x[i]);
      sum_x_square[i].increase_value(y, x[i] * x[i]);
    }
    count_y.increase_value(y, 1);

    for (int i = 0; i < dimension; i++)
    {
      //if (x[i] < lowerboundary[i] + 0.000001 || x[i] > upperboundary[i] - 0.000001)
      // adapt colvars precision
      if (x[i] < lowerboundary[i] + 0.00001 || x[i] > upperboundary[i] - 0.00001)
        return false;
    }
    distribution_x_y.increase_value(x, y, 1);

    return true;
  }

  // update the output_filename
  void update_output_filename(const std::string& filename)
  {
    output_filename = filename;
  }

private:
  std::vector<n_vector<double> > sum_x;                        // the sum of x in each y bin
  std::vector<n_vector<double> > sum_x_square;                 // the sum of x in each y bin
  n_vector<int> count_y;                              // the distribution of y
  n_matrix distribution_x_y;   // the distribution of <x, y> pair

  int dimension;

  std::vector<double> lowerboundary;
  std::vector<double> upperboundary;
  std::vector<double> width;
  std::vector<double> krestr;
  std::string output_filename;
  int output_freq;
  bool restart;
  std::vector<std::string> input_filename;
  double temperature;

  n_vector<std::vector<double> > grad;
  n_vector<int> count;

  n_vector<double> oneD_pmf;

  n_vector<std::vector<double> > input_grad;
  n_vector<int> input_count;

  // used in double integration
  std::vector<n_vector<double> > x_av;
  std::vector<n_vector<double> > sigma_square;

  bool written;
  bool written_1D;

  // calculate gradients from the internal variables
  void calc_pmf()
  {
    int norm;

    std::vector<double> loop_flag(dimension, 0);
    for (int i = 0; i < dimension; i++)
    {
      loop_flag[i] = lowerboundary[i] - HALF_Y_SIZE * width[i];
    }

    while (true)
    {
      norm = count_y.get_value(loop_flag) > 0 ? count_y.get_value(loop_flag) : 1;
      for (int i = 0; i < dimension; i++)
      {
        x_av[i].set_value(loop_flag, sum_x[i].get_value(loop_flag) / norm);
        sigma_square[i].set_value(loop_flag, sum_x_square[i].get_value(loop_flag) / norm - x_av[i].get_value(loop_flag) * x_av[i].get_value(loop_flag));
      }

      // iterate over any dimensions
      int i = dimension - 1;
      while (true)
      {
        loop_flag[i] += width[i];
        if (loop_flag[i] > upperboundary[i] + HALF_Y_SIZE * width[i] - width[i] + 0.00001)
        {
          loop_flag[i] = lowerboundary[i] - HALF_Y_SIZE * width[i];
          i--;
          if (i < 0)
            goto LOOPEND;
        }
        else
          break;
      }
    }
LOOPEND:

    // double integration
    std::vector<double> av(dimension, 0);
    std::vector<double> diff_av(dimension, 0);

    std::vector<double> loop_flag_x(dimension, 0);
    std::vector<double> loop_flag_y(dimension, 0);
    for (int i = 0; i < dimension; i++)
    {
      loop_flag_x[i] = lowerboundary[i];
      loop_flag_y[i] = loop_flag_x[i] - HALF_Y_SIZE * width[i];
    }

    while (true)
    {
      norm = 0;
      for (int i = 0; i < dimension; i++)
      {
        av[i] = 0;
        diff_av[i] = 0;
        loop_flag_y[i] = loop_flag_x[i] - HALF_Y_SIZE * width[i];
      }

      while (true)
      {
        //std::cout<<"pppppppppppppppppppppp "<<loop_flag_x[0]<<" "<<loop_flag_x[1]<<" "<<loop_flag_y[0]<<" "<<loop_flag_y[1]<<std::endl;
        norm += distribution_x_y.get_value(loop_flag_x, loop_flag_y);
        for (int i = 0; i < dimension; i++)
        {
          if (sigma_square[i].get_value(loop_flag_y) > 0.00001 || sigma_square[i].get_value(loop_flag_y) < -0.00001)
            av[i] += distribution_x_y.get_value(loop_flag_x, loop_flag_y) * ( (loop_flag_x[i] + 0.5 * width[i]) - x_av[i].get_value(loop_flag_y)) / sigma_square[i].get_value(loop_flag_y);

          diff_av[i] += distribution_x_y.get_value(loop_flag_x, loop_flag_y) * (loop_flag_x[i] - loop_flag_y[i]);
        }

        // iterate over any dimensions
        int i = dimension - 1;
        while (true)
        {
          loop_flag_y[i] += width[i];
          if (loop_flag_y[i] > loop_flag_x[i] + HALF_Y_SIZE * width[i] - width[i] + 0.00001)
          {
            loop_flag_y[i] = loop_flag_x[i] - HALF_Y_SIZE * width[i];
            i--;
            if (i < 0)
              goto LOOPEND2;
          }
          else
            break;
        }
      }
LOOPEND2:

      std::vector<double> grad_temp(dimension, 0);
      for (int i = 0; i < dimension; i++)
      {
        diff_av[i] /= (norm > 0 ? norm : 1);
        av[i] = BOLTZMANN * temperature * av[i] / (norm > 0 ? norm : 1);
        grad_temp[i] = av[i] - krestr[i] * diff_av[i];
      }
      grad.set_value(loop_flag_x, grad_temp);
      count.set_value(loop_flag_x, norm);

      // iterate over any dimensions
      int i = dimension - 1;
      while (true)
      {
        loop_flag_x[i] += width[i];
        if (loop_flag_x[i] > upperboundary[i] - width[i] + 0.00001)
        {
          loop_flag_x[i] = lowerboundary[i];
          i--;
          if (i < 0)
            goto LOOPEND3;
        }
        else
          break;
      }
    }
LOOPEND3:;
  }


  // calculate 1D pmf
  void calc_1D_pmf()
  {
    std::vector<double> last_position(1, 0);
    std::vector<double> position(1, 0);

    double min = 0;
    double dG = 0;

    oneD_pmf.set_value(lowerboundary, 0);
    last_position = lowerboundary;
    for (double i = lowerboundary[0] + width[0]; i < upperboundary[0] + 0.000001; i += width[0])
    {
      position[0] = i + 0.000001;
      if (restart == false || input_count.get_value(last_position) == 0)
      {
        dG = oneD_pmf.get_value(last_position) + grad.get_value(last_position)[0] * width[0];
      }
      else
      {
        dG = oneD_pmf.get_value(last_position) + ((grad.get_value(last_position)[0] * count.get_value(last_position) + input_grad.get_value(last_position)[0] * input_count.get_value(last_position)) / (count.get_value(last_position) + input_count.get_value(last_position))) * width[0];
      }
      if (dG < min)
        min = dG;
      oneD_pmf.set_value(position, dG);
      last_position[0] = i + 0.000001;
    }

    for (double i = lowerboundary[0]; i < upperboundary[0] + 0.000001; i += width[0])
    {
      position[0] = i + 0.000001;
      oneD_pmf.set_value(position, oneD_pmf.get_value(position) - min);
    }
  }

  // write 1D pmf
  void write_1D_pmf()
  {
    std::string pmf_filename = output_filename + ".UI.pmf";

    // only for colvars module!
// 			if (written_1D) cvm::backup_file(pmf_filename.c_str());

    std::ofstream ofile_pmf(pmf_filename.c_str());

    std::vector<double> position(1, 0);
    for (double i = lowerboundary[0]; i < upperboundary[0] + 0.000001; i += width[0])
    {
      ofile_pmf << i << " ";
      position[0] = i + 0.000001;
      ofile_pmf << oneD_pmf.get_value(position) << std::endl;
    }
    ofile_pmf.close();

    written_1D = true;
  }

  // write heads of the output files
  void writehead(std::ofstream& os) const
  {
    os << "# " << dimension << std::endl;
    for (int i = 0; i < dimension; i++)
    {
      os.precision(std::numeric_limits<double>::max_digits10);
      os << "# " << std::fixed << lowerboundary[i] << " " << width[i] << " " << int((upperboundary[i] - lowerboundary[i]) / width[i] + 0.000001) << " " << 0 << std::endl;
    }
    os << std::endl;
  }

  // write interal data, used for testing
//   void write_interal_data()
//   {
//     std::string internal_filaname = output_filename + ".UI.internal";
//
//     std::ofstream ofile_internal(internal_filaname.c_str());
//
//     std::vector<double> loop_flag(dimension, 0);
//     for (int i = 0; i < dimension; i++)
//     {
//       loop_flag[i] = lowerboundary[i];
//     }
//     while (true)
//     {
//       for (int i = 0; i < dimension; i++)
//       {
//         ofile_internal << loop_flag[i] + 0.5 * width[i] << " ";
//       }
//
//       for (int i = 0; i < dimension; i++)
//       {
//         ofile_internal << grad.get_value(loop_flag)[i] << " ";
//       }
//
//       std::vector<double> ii(dimension,0);
//       for (double i = loop_flag[0] - 10; i < loop_flag[0] + 10 + 0.00001; i+= width[0])
//       {
//         for (double j = loop_flag[1] - 10; j< loop_flag[1] + 10 +0.00001; j+=width[1])
//         {
//           ii[0] = i;
//           ii[1] = j;
//           ofile_internal << i <<" "<<j<<" "<< distribution_x_y.get_value(loop_flag,ii)<< " ";
//         }
//       }
//       ofile_internal << std::endl;
//
//       // iterate over any dimensions
//       int i = dimension - 1;
//       while (true)
//       {
//         loop_flag[i] += width[i];
//         if (loop_flag[i] > upperboundary[i] - width[i] + 0.00001)
//         {
//           loop_flag[i] = lowerboundary[i];
//           i--;
//           if (i < 0)
//             goto LOOPEND5;
//         }
//         else
//           break;
//       }
//     }
// LOOPEND5:
//     ofile_internal.close();
//   }

  // write output files
  void write_files()
  {
    std::string grad_filename = output_filename + ".UI.grad";
    std::string hist_filename = output_filename + ".UI.hist.grad";
    std::string count_filename = output_filename + ".UI.count";

    // only for colvars module!
// 			if (written) cvm::backup_file(grad_filename.c_str());
    //if (written) cvm::backup_file(hist_filename.c_str());
// 			if (written) cvm::backup_file(count_filename.c_str());

    std::ofstream ofile(grad_filename.c_str());
    std::ofstream ofile_hist(hist_filename.c_str(), std::ios::app);
    std::ofstream ofile_count(count_filename.c_str());

    writehead(ofile);
    writehead(ofile_hist);
    writehead(ofile_count);

    if (dimension == 1)
    {
      calc_1D_pmf();
      write_1D_pmf();
    }

    std::vector<double> loop_flag(dimension, 0);
    for (int i = 0; i < dimension; i++)
    {
      loop_flag[i] = lowerboundary[i];
    }
    while (true)
    {
      for (int i = 0; i < dimension; i++)
      {
        ofile << std::fixed << std::setprecision(9) << loop_flag[i] + 0.5 * width[i] << " ";
        ofile_hist << std::fixed << std::setprecision(9) << loop_flag[i] + 0.5 * width[i] << " ";
        ofile_count << std::fixed << std::setprecision(9) << loop_flag[i] + 0.5 * width[i] << " ";
      }

      if (restart == false)
      {
        for (int i = 0; i < dimension; i++)
        {
          ofile << std::fixed << std::setprecision(9) << grad.get_value(loop_flag)[i] << " ";
          ofile_hist << std::fixed << std::setprecision(9) << grad.get_value(loop_flag)[i] << " ";
        }
        ofile << std::endl;
        ofile_hist << std::endl;
        ofile_count << count.get_value(loop_flag) << " " <<std::endl;
      }
      else
      {
        double final_grad = 0;
        for (int i = 0; i < dimension; i++)
        {
          int total_count_temp = (count.get_value(loop_flag) + input_count.get_value(loop_flag));
          if (input_count.get_value(loop_flag) == 0)
            final_grad = grad.get_value(loop_flag)[i];
          else
            final_grad = ((grad.get_value(loop_flag)[i] * count.get_value(loop_flag) + input_grad.get_value(loop_flag)[i] * input_count.get_value(loop_flag)) / total_count_temp);
          ofile << std::fixed << std::setprecision(9) << final_grad << " ";
          ofile_hist << std::fixed << std::setprecision(9) << final_grad << " ";
        }
        ofile << std::endl;
        ofile_hist << std::endl;
        ofile_count << (count.get_value(loop_flag) + input_count.get_value(loop_flag)) << " " <<std::endl;
      }

      // iterate over any dimensions
      int i = dimension - 1;
      while (true)
      {
        loop_flag[i] += width[i];
        if (loop_flag[i] > upperboundary[i] - width[i] + 0.00001)
        {
          loop_flag[i] = lowerboundary[i];
          i--;
          ofile << std::endl;
          ofile_hist << std::endl;
          ofile_count << std::endl;
          if (i < 0)
            goto LOOPEND4;
        }
        else
          break;
      }
    }
LOOPEND4:
    ofile.close();
    ofile_count.close();
    ofile_hist.close();

    written = true;
  }

  // read input files
  void read_inputfiles(const std::vector<std::string>& input_filename)
  {
    char sharp;
    double nothing;
    int dimension_temp;

    std::vector<double> loop_bin_size(dimension, 0);
    std::vector<double> position_temp(dimension, 0);
    std::vector<double> grad_temp(dimension, 0);
    int count_temp = 0;
    for (unsigned int i = 0; i < input_filename.size(); i++)
    {
      int size = 1, size_temp = 0;

      std::string count_filename = input_filename[i] + ".UI.count";
      std::string grad_filename = input_filename[i] + ".UI.grad";

      std::ifstream count_file(count_filename.c_str(), std::ios::in);
      std::ifstream grad_file(grad_filename.c_str(), std::ios::in);

      count_file >> sharp >> dimension_temp;
      grad_file >> sharp >> dimension_temp;

      for (int j = 0; j < dimension; j++)
      {
        count_file >> sharp >> nothing >> nothing >> size_temp >> nothing;
        grad_file >> sharp >> nothing >> nothing >> nothing >> nothing;
        size *= size_temp;
      }

      for (int j = 0; j < size; j++)
      {
        do
        {
          for (int k = 0; k < dimension; k++)
          {
            count_file >> position_temp[k];
            grad_file >> nothing;
          }

          for (int l = 0; l < dimension; l++)
          {
            grad_file >> grad_temp[l];
          }
          count_file >> count_temp;
        }
        while (position_temp[i] < lowerboundary[i] - 0.000001 || position_temp[i] > upperboundary[i] + 0.000001);

        if (count_temp == 0)
        {
          continue;
        }

        for (int m = 0; m < dimension; m++)
        {
          grad_temp[m] = (grad_temp[m] * count_temp + input_grad.get_value(position_temp)[m] * input_count.get_value(position_temp)) / (count_temp + input_count.get_value(position_temp));
        }
        input_grad.set_value(position_temp, grad_temp);
        input_count.increase_value(position_temp, count_temp);
      }

      count_file.close();
      grad_file.close();
    }
  }
};
}

}
}

#endif
