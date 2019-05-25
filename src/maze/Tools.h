/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2019 Jakub Rydzewski (jr@fizyka.umk.pl). All rights reserved.

See http://www.maze-code.github.io for more information.

This file is part of maze.

maze is free software: you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option)
any later version.

maze is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.

See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with maze. If not, see <https://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_maze_Tools_h
#define __PLUMED_maze_Tools_h

/**
 * @file Tools.h
 *
 * @author J. Rydzewski
 */

#include "core/ActionSet.h"
#include "tools/Vector.h"
#include "Core.h"

namespace PLMD {
namespace maze {

/**
 * @class tls Tools.h "maze/Tools.h"
 *
 * @brief Helper functions.
 */
class tls {
public:
  template<typename T> static int sgn(T val);

  template<typename T>
  static std::vector<std::string> get_labels_actions(const ActionSet&);

  template<typename T>
  static T get_pointer_label(
    const std::string&,
    const ActionSet&, std::string&
  );

  template<typename T>
  static std::vector<T> get_pointers_labels(
    const std::vector<std::string>&,
    const ActionSet&, std::string&
  );

  template<typename T>
  static std::vector<T> Vector2vector(const Vector&);

  template<typename T>
  static Vector vector2Vector(const std::vector<T>&);

  template<typename T>
  static T vector_l(const std::vector<T>&);

  template<typename T>
  static std::vector<T> vector_n(const std::vector<T>&);

  template<typename T>
  static T mean(const std::vector<T>&);

  template<typename T>
  static T std(const std::vector<T>&);

  struct delete_ptr {
    template<typename T> void operator()(T t) {
      delete t;
    }
  };
};

template<typename T>
T tls::std(const std::vector<T>& v) {
  T m=mean(v);
  T sq_sum=std::inner_product(v.begin(), v.end(), v.begin(), 0.0);

  return std::sqrt(sq_sum/v.size()-m*m);
}

template<typename T>
T tls::mean(const std::vector<T>& v) {
  return std::accumulate(v.begin(), v.end(), 0.0)/v.size();
}

template<typename T>
T tls::vector_l(const std::vector<T>& v) {
  return std::sqrt(std::inner_product(v.begin(), v.end(), v.begin(), 0.0));
}

template<typename T>
std::vector<T> tls::vector_n(const std::vector<T>& v) {
  double l=vector_l(v);
  std::vector<double> n;
  for(std::size_t i=0; i<v.size(); ++i)
    n.push_back(v[i]/l);

  return n;
}

template<typename T>
std::vector<T> tls::Vector2vector(const Vector& v) {
  std::vector<T> t= {v[0], v[1], v[2]};

  return t;
}

template<typename T>
Vector tls::vector2Vector(const std::vector<T>& v) {
  Vector t(v[0], v[1], v[2]);

  return t;
}

template<typename T>
int tls::sgn(T val) {
  return (T(0)<val)-(val<T(0));
}

template<typename T>
std::vector<std::string> tls::get_labels_actions(const ActionSet& actionset) {
  std::vector<std::string> action_str(0);
  std::vector<T> action_pntrs=actionset.select<T>();

  for(unsigned int i=0; i<action_pntrs.size(); i++)
    action_str.push_back(action_pntrs[i]->getLabel());

  return action_str;
}

template<typename T> T
tls::get_pointer_label(
  const std::string& action_label,
  const ActionSet& actionset,
  std::string& error_msg) {

  std::vector<std::string> action_labels(1);
  action_labels[0]=action_label;
  std::vector<T> action_pntrs=get_pointers_labels<T>(action_labels, actionset, error_msg);

  return action_pntrs[0];
}

template<typename T>
std::vector<T> tls::get_pointers_labels(
  const std::vector<std::string>& action_labels,
  const ActionSet& actionset,
  std::string& error_msg) {

  std::vector<T> action_pntrs(action_labels.size(), NULL);
  error_msg="";
  std::vector<std::string> missing(0);

  for(unsigned int i=0; i<action_labels.size(); i++) {
    action_pntrs[i]=actionset.selectWithLabel<T>(action_labels[i]);
    if(action_pntrs[i]==NULL)
      missing.push_back(action_labels[i]);
  }

  return action_pntrs;
}

} // namespace maze
} // namespace PLMD

#endif // __PLUMED_maze_Tools_h
