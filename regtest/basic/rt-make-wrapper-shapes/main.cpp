#include <array>
#include <cstdlib>
#include <fstream>
#include <string>
#include <type_traits>
#include <vector>
struct WrapperVec3 {
  double data[3];
};
struct WrapperRows {
  double data[6];
};
struct WrapperDualArray {
  double data[6];
};
struct WrapperOverloadedAddress {
  double values[2];
  WrapperOverloadedAddress* operator&() {
    return 0;
  }
  const WrapperOverloadedAddress* operator&() const {
    return 0;
  }
};
struct WrapperConstFalse {
  double values[2];
  std::size_t size() const {
    return 1;
  }
  double* data() {
    return values;
  }
  const double* data() const {
    return values;
  }
};
struct WrapperMutableBuffer {
  double values[2];
  std::size_t size() & { return 2; }
  double* data() & { return values; }
};
namespace PLMD {
#if defined(__PLUMED_WRAPPER_CXX_ANONYMOUS_NAMESPACE) && __PLUMED_WRAPPER_CXX_ANONYMOUS_NAMESPACE
namespace {
#endif
namespace wrapper {
template<typename T> struct is_custom_array;
template<> struct is_custom_array<WrapperVec3> : std::true_type {
  typedef double value_type;
};
template<> struct is_custom_array<WrapperRows> : std::true_type {
  typedef double value_type[3];
};
template<> struct is_custom_array<WrapperOverloadedAddress> : std::true_type {
  typedef double value_type;
};
template<> struct is_custom_array<WrapperDualArray> : std::true_type {
  typedef double value_type;
};
template<> struct is_custom_array<const WrapperDualArray> : std::true_type {
  typedef double value_type[3];
};
template<> struct is_custom_array<const double[6]> : std::true_type {
  typedef double value_type[3];
};
template<> struct is_custom_array<std::array<double,6> > : std::true_type {
  typedef double value_type[3];
};
template<> struct is_custom_array<const std::array<double,6> > : std::true_type {
  typedef double value_type[2];
};
template<> struct is_custom_array<WrapperMutableBuffer> : std::false_type {};
template<> struct is_custom_array<WrapperConstFalse> : std::true_type {
  typedef double value_type;
};
template<> struct is_custom_array<const WrapperConstFalse> : std::false_type {};
}
#if defined(__PLUMED_WRAPPER_CXX_ANONYMOUS_NAMESPACE) && __PLUMED_WRAPPER_CXX_ANONYMOUS_NAMESPACE
}
#endif
}
#include "plumed/wrapper/Plumed.h"
static_assert(PLMD::wrapper::is_custom_array<const double[2][3]>::value,
              "const multidimensional C arrays must remain directly classifiable");
static_assert(std::is_same<PLMD::wrapper::is_custom_array<const double[3]>::value_type,
              const double>::value,
              "const C-array traits must retain their scalar const qualification");
static_assert(std::is_same<PLMD::wrapper::is_custom_array<const double[2][3]>::value_type,
              const double[3]>::value,
              "nested const C-array traits must retain the row const qualification");

static plumed_safeptr captured;
static std::size_t captured_shape[5];
extern "C" plumed plumed_create_reference(plumed p) {
  return p;
}
extern "C" void plumed_finalize(plumed) {}
extern "C" void* plumed_malloc(std::size_t n) {
  return std::malloc(n);
}
extern "C" void plumed_free(void* p) {
  std::free(p);
}
extern "C" void plumed_cmd_safe_nothrow(plumed, const char*, plumed_safeptr s, plumed_nothrow_handler) {
  captured=s;
  for(unsigned i=0; i<5; i++) {
    captured_shape[i]=0;
  }
  if(s.shape) {
    for(unsigned i=0; i<4 && s.shape[i]!=0; i++) {
      captured_shape[i]=s.shape[i];
    }
  }
}

static bool shape_is(std::size_t a, std::size_t b=0, std::size_t c=0) {
  return captured.shape && captured_shape[0]==a && captured_shape[1]==b && captured_shape[2]==c;
}

int main() {
  PLMD::Plumed p(plumed{0});
  const double matrix[2][3]= {{}};
  p.cmd("const-matrix",matrix);
  if(captured.ptr!=matrix || !shape_is(2,3) || (captured.flags&0x0e000000u)!=0x06000000u ||
      (captured.flags&0x10000000u)) {
    return 1;
  }
  const double customized_const_array[6]= {};
  p.cmd("const-array-hook",customized_const_array);
  if(captured.ptr!=customized_const_array || !shape_is(2,3)) {
    return 2;
  }

  const std::array<std::array<double,3>,2> nested= {};
  p.cmd("nested",nested);
  if(captured.ptr!=nested.data()->data() || !shape_is(2,3)) {
    return 3;
  }
  std::array<double,0> zero;
  p.cmd("zero",zero);
  if(captured.ptr!=zero.data() || !shape_is(0)) {
    return 4;
  }
  std::array<double,6> hooked_array= {{}};
  p.cmd("hooked-array",hooked_array);
  if(captured.ptr!=hooked_array.data() || !shape_is(2,3)) {
    return 5;
  }
  const std::array<double,6> const_hooked_array= {{}};
  p.cmd("const-hooked-array",const_hooked_array);
  if(captured.ptr!=const_hooked_array.data() || !shape_is(3,2)) {
    return 6;
  }

  const WrapperVec3 vectors[2]= {};
  p.cmd("custom-vector",vectors);
  if(!shape_is(2,3)) {
    return 7;
  }
  WrapperRows rows= {};
  p.cmd("custom-rows",rows);
  if(!shape_is(2,3)) {
    return 8;
  }
  WrapperOverloadedAddress overloaded_address= {{}};
  p.cmd("overloaded-address",overloaded_address);
  if(captured.ptr!=overloaded_address.values || !shape_is(2)) {
    return 9;
  }
  WrapperDualArray mutable_dual= {{}}, const_dual= {{}};
  p.cmd("mutable-hook",mutable_dual);
  if(!shape_is(6)) {
    return 10;
  }
  p.cmd("const-hook",static_cast<const WrapperDualArray&>(const_dual));
  if(!shape_is(2,3)) {
    return 11;
  }

  std::vector<std::array<double,3> > dynamic_fixed(2);
  p.cmd("dynamic-fixed",dynamic_fixed);
  if(captured.ptr!=dynamic_fixed.data()->data() || !shape_is(2,3)) {
    return 12;
  }
  const std::vector<double> empty;
  p.cmd("empty",empty);
  if(captured.ptr!=empty.data() || !shape_is(0)) {
    return 13;
  }
  WrapperMutableBuffer mutable_buffer= {{}};
  p.cmd("mutable-data-lvalue",mutable_buffer);
  if(captured.ptr!=mutable_buffer.values || !shape_is(2)) {
    return 14;
  }
  p.cmd("mutable-data-temporary",WrapperMutableBuffer());
  if(!shape_is(2) || !(captured.flags&0x10000000u)) {
    return 15;
  }
  const WrapperConstFalse const_false= {{}};
  p.cmd("const-false-hook",const_false);
  if(!shape_is(1) || (captured.flags&0x0e000000u)!=0x06000000u) {
    return 16;
  }

  std::string text("x");
  p.cmd("string",text);
  const std::size_t mutable_string_pointer_flag=
    std::is_same<decltype(std::declval<std::string&>().data()),char*>::value ? 0x04000000u : 0x06000000u;
  if(captured.ptr!=text.data() || !shape_is(2) ||
      (captured.flags&0x0e000000u)!=mutable_string_pointer_flag) {
    return 17;
  }
  const std::string const_text("y");
  p.cmd("const-string",const_text);
  if(captured.ptr!=const_text.data() || !shape_is(2) || (captured.flags&0x0e000000u)!=0x06000000u) {
    return 18;
  }
  p.cmd("temporary",std::array<std::array<double,2>,1>());
  if(!shape_is(1,2) || !(captured.flags&0x10000000u)) {
    return 19;
  }

  double raw=0;
  p.cmd("raw",&raw);
#if __PLUMED_WRAPPER_CXX_DETECT_SHAPES_STRICT
  if(!shape_is(1)) {
    return 26;
  }
  double* null_raw=0;
  p.cmd("null-raw",null_raw);
  if(!shape_is(0)) {
    return 27;
  }
  double fixed_storage[3]= {};
  double (*fixed_pointer)[3]=&fixed_storage;
  p.cmd("fixed-pointer",fixed_pointer);
  if(!shape_is(1,3)) {
    return 28;
  }
#else
  if(captured.shape) {
    return 24;
  }
  p.cmd("nelem",&raw,2);
  if(captured.nelem!=2 || captured.shape) {
    return 25;
  }
  double* null_raw=0;
  p.cmd("null-raw",null_raw);
  if(captured.shape) {
    return 28;
  }
  double fixed_storage[3]= {};
  double (*fixed_pointer)[3]=&fixed_storage;
  p.cmd("fixed-pointer",fixed_pointer);
  if(!shape_is(static_cast<std::size_t>(-1),3)) {
    return 29;
  }
#endif
  p.cmd("explicit",&raw, {2,3});
  if(!shape_is(2,3)) {
    return 30;
  }
  double fixed[2][3]= {{}};
  try {
    p.cmd("rank",fixed, {1,1,1,1});
    return 31;
  } catch(const PLMD::Plumed::ExceptionTypeError&) {}
  // Empty outer and inner dimensions cannot be represented as nonempty shapes.
  std::vector<std::array<double,3>> empty_rows;
  p.cmd("empty-rows",empty_rows);
  if(captured.ptr!=empty_rows.data() || !shape_is(0)) {
    return 40;
  }
  std::array<std::array<double,3>,0> zero_rows;
  p.cmd("zero-rows",zero_rows);
  if(captured.ptr!=zero_rows.data() || !shape_is(0)) {
    return 41;
  }
  std::array<std::array<double,0>,2> empty_columns;
  p.cmd("empty-columns",empty_columns);
  if(!shape_is(0)) {
    return 42;
  }
  const std::array<double,0> const_zero = {};
  p.cmd("const-zero",const_zero);
  if(captured.ptr!=const_zero.data() || !shape_is(0) ||
      (captured.flags&0x0e000000u)!=0x06000000u) {
    return 43;
  }
  const double (*null_rows)[3]=nullptr;
  p.cmd("null-rows",null_rows);
  if(captured.ptr || !shape_is(0)) {
    return 44;
  }
  p.cmd("explicit-zero",fixed, {0,3});
  if(!shape_is(0)) {
    return 45;
  }
#if !__PLUMED_WRAPPER_CXX_DETECT_SHAPES_STRICT
  p.cmd("zero-nelem",&fixed[0],0);
  if(!shape_is(0)) {
    return 46;
  }
#endif
  const WrapperRows const_rows = {};
  p.cmd("const-custom-rows",const_rows);
  if(captured.ptr!=const_rows.data || !shape_is(2,3) ||
      (captured.flags&0x0e000000u)!=0x06000000u) {
    return 47;
  }
  const WrapperOverloadedAddress const_address = {};
  p.cmd("const-address",const_address);
  if(captured.ptr!=const_address.values || !shape_is(2)) {
    return 48;
  }
  double rank4[1][2][3][4] = {};
  p.cmd("rank4",rank4);
  if(!shape_is(1,2,3) || captured_shape[3]!=4) {
    return 49;
  }
  std::vector<WrapperRows> custom_rows(2);
  p.cmd("custom-dynamic",custom_rows);
  if(captured.ptr!=custom_rows.data()->data || !shape_is(2,2,3)) {
    return 50;
  }
  std::string empty_string;
  p.cmd("empty-string",empty_string);
  if(captured.ptr!=empty_string.data() || !shape_is(1)) {
    return 51;
  }
  const std::vector<std::array<double,3>> const_dynamic_fixed(2);
  p.cmd("const-dynamic-fixed",const_dynamic_fixed);
  if(captured.ptr!=const_dynamic_fixed.data()->data() || !shape_is(2,3) ||
      (captured.flags&0x0e000000u)!=0x06000000u) {
    return 52;
  }
  p.cmd("explicit-rows",&matrix[0], {2});
  if(captured.ptr!=matrix || !shape_is(2,3)) {
    return 53;
  }
#if !__PLUMED_WRAPPER_CXX_DETECT_SHAPES_STRICT
  p.cmd("nelem-rows",&matrix[0],2);
  if(captured.ptr!=matrix || !shape_is(2,3)) {
    return 54;
  }
#endif
  std::ofstream output("output");
  output << "wrapper shape metadata checks passed\n";
  return 0;
}
