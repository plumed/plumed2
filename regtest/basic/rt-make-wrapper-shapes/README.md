# C++11 wrapper shape corner cases

The production change is confined to `src/wrapper/Plumed.h`.

| Case | Previous behavior | Corrected behavior |
| --- | --- | --- |
| Const user-defined fixed arrays, including grouped element types | The mutable customization was lost, or its cast removed const and failed to compile | Inherit the mutable trait, preserve const on the pointer, and honor exact const/false specializations |
| `std::array<T,0>` | Used `sizeof` as storage size and the object's address as element storage | Zero logical extent and the actual `data()` pointer |
| Nested empty arrays, empty vectors of fixed arrays, null pointers to fixed arrays | Recursive flattening replaced a zero extent with nonzero inner dimensions | Retain empty metadata throughout recursion |
| Custom arrays with overloaded `operator&` | Could use an unrelated or null pointer | Acquire the real object address with `std::addressof` |
| Contiguous containers with mutable-only or `&`-qualified methods | Probed an rvalue or called `size()` through const, rejecting valid inputs | Probe and access the named argument as an lvalue; temporary inputs still carry the no-retention flag |
| Nested dynamic containers, e.g. `vector<vector<double>>` | Accepted the addresses of container objects as shaped data | Reject at compile time; the elements do not form one contiguous scalar buffer |

The tests also protect const multidimensional C arrays, nested fixed/dynamic
combinations, exact trait overrides, strings (including the terminating NUL),
pointer/size/explicit-shape overloads, temporary lifetime flags, and the existing
four-dimension limit. The mock C entry point inspects metadata synchronously;
it does not require a kernel or dereference temporary/null element storage.

The one-parameter `is_custom_array<T>` customization contract is retained:
custom types must actually contain packed contiguous `value_type` storage.
No new layout restriction or PLUMED-native type adapter is introduced.
Non-contiguous containers, `vector<bool>`, nonpointer `data()` APIs and volatile
objects are not supported. Legacy raw-pointer behavior and strict-mode defaults
are unchanged. Shapes remain zero-terminated: an inferred zero axis is reported
as empty metadata, since this ABI cannot encode inner dimensions after zero.

Run the normal PLUMED regtest for runtime metadata assertions, or the independent
C++11/14/17, strict/legacy and compile-rejection matrix without configuring:

```sh
CXX=clang++ sh regtest/basic/rt-make-wrapper-shapes/check.sh
CXX=g++ sh regtest/basic/rt-make-wrapper-shapes/check.sh
```
