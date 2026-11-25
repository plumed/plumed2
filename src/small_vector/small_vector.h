/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * An implementation of `small_vector` (a vector with a small
 * buffer optimization). I would probably have preferred to
 * call this `inline_vector`, but I'll just go with the canonical
 * name for now.
 *
 * Copyright © 2020-2021 Gene Harvey
 *
 * This software may be modified and distributed under the terms
 * of the MIT license. See the LICENSE file for details.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_small_vector_small_vector_h
#define __PLUMED_small_vector_small_vector_h
/** small_vector.hpp
 * An implementation of `small_vector` (a vector with a small
 * buffer optimization). I would probably have preferred to
 * call this `inline_vector`, but I'll just go with the canonical
 * name for now.
 *
 * Copyright © 2020-2021 Gene Harvey
 *
 * This software may be modified and distributed under the terms
 * of the MIT license. See the LICENSE file for details.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef PLUMED_GCH_SMALL_VECTOR_HPP
#define PLUMED_GCH_SMALL_VECTOR_HPP

#ifdef __clang__
#  ifndef PLUMED_GCH_CLANG
#    define PLUMED_GCH_CLANG
#  endif
#  if defined (__cplusplus) && __cplusplus >= 201103L
#    ifndef PLUMED_GCH_CLANG_11
#      define PLUMED_GCH_CLANG_11
#    endif
#  endif
#  if defined (__cplusplus) && __cplusplus >= 201402L
#    ifndef PLUMED_GCH_CLANG_14
#      define PLUMED_GCH_CLANG_14
#    endif
#  endif
#  if defined (__cplusplus) && __cplusplus >= 201703L
#    ifndef PLUMED_GCH_CLANG_17
#      define PLUMED_GCH_CLANG_17
#    endif
#  endif
#  if defined (__cplusplus) && __cplusplus >= 202002L
#    ifndef PLUMED_GCH_CLANG_20
#      define PLUMED_GCH_CLANG_20
#    endif
#  endif
#endif

#ifndef PLUMED_GCH_CPP14_CONSTEXPR
#  if defined (__cpp_constexpr) && __cpp_constexpr >= 201304L
#    define PLUMED_GCH_CPP14_CONSTEXPR constexpr
#    ifndef PLUMED_GCH_HAS_CPP14_CONSTEXPR
#      define PLUMED_GCH_HAS_CPP14_CONSTEXPR
#    endif
#  else
#    define PLUMED_GCH_CPP14_CONSTEXPR
#  endif
#endif

#ifndef PLUMED_GCH_CPP17_CONSTEXPR
#  if defined (__cpp_constexpr) && __cpp_constexpr >= 201603L
#    define PLUMED_GCH_CPP17_CONSTEXPR constexpr
#    ifndef PLUMED_GCH_HAS_CPP17_CONSTEXPR
#      define PLUMED_GCH_HAS_CPP17_CONSTEXPR
#    endif
#  else
#    define PLUMED_GCH_CPP17_CONSTEXPR
#  endif
#endif

#ifndef PLUMED_GCH_CPP20_CONSTEXPR
#  if defined (__cpp_constexpr) && __cpp_constexpr >= 201907L
#    define PLUMED_GCH_CPP20_CONSTEXPR constexpr
#    ifndef PLUMED_GCH_HAS_CPP20_CONSTEXPR
#      define PLUMED_GCH_HAS_CPP20_CONSTEXPR
#    endif
#  else
#    define PLUMED_GCH_CPP20_CONSTEXPR
#  endif
#endif

#ifndef PLUMED_GCH_NORETURN
#  if defined (__has_cpp_attribute) && __has_cpp_attribute (noreturn) >= 200809L
#    define PLUMED_GCH_NORETURN [[noreturn]]
#  else
#    define PLUMED_GCH_NORETURN
#  endif
#endif

#ifndef PLUMED_GCH_NODISCARD
#  if defined (__has_cpp_attribute) && __has_cpp_attribute (nodiscard) >= 201603L
#    if ! defined (__clang__) || defined (PLUMED_GCH_CLANG_17)
#      define PLUMED_GCH_NODISCARD [[nodiscard]]
#    else
#      define PLUMED_GCH_NODISCARD
#    endif
#  else
#    define PLUMED_GCH_NODISCARD
#  endif
#endif

#ifndef PLUMED_GCH_INLINE_VARIABLE
#  if defined (__cpp_inline_variables) && __cpp_inline_variables >= 201606L
#    define PLUMED_GCH_INLINE_VARIABLE inline
#  else
#    define PLUMED_GCH_INLINE_VARIABLE
#  endif
#endif

#ifndef PLUMED_GCH_EMPTY_BASE
#  if defined (_MSC_FULL_VER) && _MSC_FULL_VER >= 190023918L
#    define PLUMED_GCH_EMPTY_BASE __declspec (empty_bases)
#  else
#    define PLUMED_GCH_EMPTY_BASE
#  endif
#endif

#ifndef PLUMED_GCH_IMPLICIT_CONVERSION
#  if defined (__cpp_conditional_explicit) && __cpp_conditional_explicit >= 201806L
#    define PLUMED_GCH_IMPLICIT_CONVERSION explicit (false)
#  else
#    define PLUMED_GCH_IMPLICIT_CONVERSION /* implicit */
#  endif
#endif

#if defined (__cpp_variable_templates) && __cpp_variable_templates >= 201304L
#  ifndef PLUMED_GCH_VARIABLE_TEMPLATES
#    define PLUMED_GCH_VARIABLE_TEMPLATES
#  endif
#endif

#if defined (__cpp_deduction_guides) && __cpp_deduction_guides >= 201703L
#  ifndef PLUMED_GCH_CTAD_SUPPORT
#    define PLUMED_GCH_CTAD_SUPPORT
#  endif
#endif

#if defined (__cpp_if_constexpr) && __cpp_if_constexpr >= 201606L
#  ifndef PLUMED_GCH_CONSTEXPR_IF
#    define PLUMED_GCH_CONSTEXPR_IF
#  endif
#endif

#if defined (__cpp_exceptions) && __cpp_exceptions >= 199711L
#  ifndef PLUMED_GCH_EXCEPTIONS
#    define PLUMED_GCH_EXCEPTIONS
#  endif
#endif

#ifndef PLUMED_GCH_TRY
#  ifdef PLUMED_GCH_EXCEPTIONS
#    define PLUMED_GCH_TRY try
#  else
#    ifdef PLUMED_GCH_CONSTEXPR_IF
#      define PLUMED_GCH_TRY if constexpr (true)
#    else
#      define PLUMED_GCH_TRY if (true)
#    endif
#  endif
#endif

#ifndef PLUMED_GCH_CATCH
#  ifdef PLUMED_GCH_EXCEPTIONS
#    define PLUMED_GCH_CATCH(...) catch (__VA_ARGS__)
#  else
#    ifdef PLUMED_GCH_CONSTEXPR_IF
#      define PLUMED_GCH_CATCH(...) else if constexpr (false)
#    else
#      define PLUMED_GCH_CATCH(...) else if (false)
#    endif
#  endif
#endif

#ifndef PLUMED_GCH_THROW
#  ifdef PLUMED_GCH_EXCEPTIONS
#    define PLUMED_GCH_THROW throw
#  else
#    define PLUMED_GCH_THROW
#  endif
#endif

#ifndef PLUMED_GCH_CONSTEVAL
#  if defined (__cpp_consteval) && __cpp_consteval >= 201811L
#    define PLUMED_GCH_CONSTEVAL consteval
#    ifndef PLUMED_GCH_HAS_CONSTEVAL
#      define PLUMED_GCH_HAS_CONSTEVAL
#    endif
#  else
#    define PLUMED_GCH_CONSTEVAL constexpr
#  endif
#endif

#if defined (__cpp_impl_three_way_comparison) && __cpp_impl_three_way_comparison >= 201907L
#  ifndef PLUMED_GCH_IMPL_THREE_WAY_COMPARISON
#    define PLUMED_GCH_IMPL_THREE_WAY_COMPARISON
#  endif
#endif

#if defined (__cpp_concepts) && __cpp_concepts >= 201907L
#  ifndef PLUMED_GCH_CONCEPTS
#    define PLUMED_GCH_CONCEPTS
#  endif
#endif

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <initializer_list>
#include <iterator>
#include <limits>
#include <memory>
#include <new>
#include <type_traits>
#include <utility>

#ifdef PLUMED_GCH_IMPL_THREE_WAY_COMPARISON
#  if defined (__has_include) && __has_include (<compare>)
#    include <compare>
#  endif
#endif

#ifdef PLUMED_GCH_CONCEPTS
#  if defined (__has_include) && __has_include (<concepts>)
#    include <concepts>
#  endif
#endif

#ifdef PLUMED_GCH_STDLIB_INTEROP
#  include <array>
#  include <valarray>
#  include <vector>
#endif

#ifdef PLUMED_GCH_EXCEPTIONS
#  include <stdexcept>
#else
#  include <cstdio>
#  include <cstdlib>
#endif

#if defined (__cpp_lib_three_way_comparison) && __cpp_lib_three_way_comparison >= 201907L
#  ifndef PLUMED_GCH_LIB_THREE_WAY_COMPARISON
#    define PLUMED_GCH_LIB_THREE_WAY_COMPARISON
#  endif
#endif

#if defined (__cpp_lib_concepts) && __cpp_lib_concepts >= 202002L
#  if ! defined (PLUMED_GCH_LIB_CONCEPTS) && ! defined (PLUMED_GCH_DISABLE_CONCEPTS)
#    define PLUMED_GCH_LIB_CONCEPTS
#  endif
#endif

#if defined (__cpp_lib_is_final) && __cpp_lib_is_final >= 201402L
#  ifndef PLUMED_GCH_LIB_IS_FINAL
#    define PLUMED_GCH_LIB_IS_FINAL
#  endif
#endif

#if defined (__cpp_lib_is_constant_evaluated) && __cpp_lib_is_constant_evaluated >= 201811L
#  ifndef PLUMED_GCH_LIB_IS_CONSTANT_EVALUATED
#    define PLUMED_GCH_LIB_IS_CONSTANT_EVALUATED
#  endif
#endif

#if defined (__cpp_lib_is_swappable) && __cpp_lib_is_swappable >= 201603L
#  ifndef PLUMED_GCH_LIB_IS_SWAPPABLE
#    define PLUMED_GCH_LIB_IS_SWAPPABLE
#  endif
#endif

#if defined (__cpp_lib_allocator_traits_is_always_equal)
#  if __cpp_lib_allocator_traits_is_always_equal >= 201411L
#    ifndef PLUMED_GCH_LIB_IS_ALWAYS_EQUAL
#      define PLUMED_GCH_LIB_IS_ALWAYS_EQUAL
#    endif
#  endif
#endif

#if defined (__cpp_lib_constexpr_memory) && __cpp_lib_constexpr_memory >= 201811L
#  ifndef PLUMED_GCH_LIB_CONSTEXPR_MEMORY
#    define PLUMED_GCH_LIB_CONSTEXPR_MEMORY
#  endif
#endif

// TODO:
//   Make sure we don't need any laundering in the internal class functions.
//   I also need some sort of test case to actually show where UB is occurring,
//   because it's still a bit unclear to me.
#if defined (__cpp_lib_launder) && __cpp_lib_launder >= 201606L
#  ifndef PLUMED_GCH_LIB_LAUNDER
#    define PLUMED_GCH_LIB_LAUNDER
#  endif
#endif

// defined if the entire thing is available for constexpr
#ifndef PLUMED_GCH_SMALL_VECTOR_CONSTEXPR
#  if defined (PLUMED_GCH_HAS_CPP20_CONSTEXPR) && defined (PLUMED_GCH_LIB_IS_CONSTANT_EVALUATED) \
                                        && defined (PLUMED_GCH_LIB_CONSTEXPR_MEMORY)
#    define PLUMED_GCH_SMALL_VECTOR_CONSTEXPR constexpr
#    ifndef PLUMED_GCH_HAS_CONSTEXPR_SMALL_VECTOR
#      define PLUMED_GCH_HAS_CONSTEXPR_SMALL_VECTOR
#    endif
#  else
#    define PLUMED_GCH_SMALL_VECTOR_CONSTEXPR
#  endif
#endif

#ifndef PLUMED_GCH_SMALL_VECTOR_DEFAULT_SIZE
#  define PLUMED_GCH_SMALL_VECTOR_DEFAULT_SIZE 64
#endif

namespace PLMD {
namespace gch {

#ifdef PLUMED_GCH_LIB_CONCEPTS

namespace concepts {

template <typename T>
concept Complete = requires { sizeof (T); };

// Note: this mirrors the named requirements, not the standard concepts, so we don't require
// the destructor to be noexcept for Destructible.
template <typename T>
concept Destructible = std::is_destructible<T>::value;

template <typename T>
concept TriviallyDestructible = std::is_trivially_destructible<T>::value;

template <typename T>
concept NoThrowDestructible = std::is_nothrow_destructible<T>::value;

// Note: this mirrors the named requirements, not the standard library concepts,
// so we don't require Destructible here.

template <typename T, typename... Args>
concept ConstructibleFrom = std::is_constructible<T, Args...>::value;

template <typename T, typename... Args>
concept NoThrowConstructibleFrom = std::is_nothrow_constructible<T, Args...>::value;

template <typename From, typename To>
concept ConvertibleTo =
  std::is_convertible<From, To>::value
&&  requires (typename std::add_rvalue_reference<From>::type (&f) (void)) {
  static_cast<To> (f ());
};

template <typename From, typename To>
concept NoThrowConvertibleTo =
  std::is_nothrow_convertible<From, To>::value
&&  requires (typename std::add_rvalue_reference<From>::type (&f) (void) noexcept) {
  {
    static_cast<To> (f ())
  }
  noexcept;
};

// Note: std::default_initializable requires std::destructible.
template <typename T>
concept DefaultConstructible =
  ConstructibleFrom<T>
  &&  requires { T { }; }
  &&  requires { ::new (static_cast<void *> (nullptr)) T; };

template <typename T>
concept MoveAssignable = std::assignable_from<T&, T>;

template <typename T>
concept CopyAssignable =
  MoveAssignable<T>
  &&  std::assignable_from<T&, T&>
  &&  std::assignable_from<T&, const T&>
  &&  std::assignable_from<T&, const T>;

template <typename T>
concept MoveConstructible = ConstructibleFrom<T, T> && ConvertibleTo<T, T>;

template <typename T>
concept NoThrowMoveConstructible =
  NoThrowConstructibleFrom<T, T>
  &&  NoThrowConvertibleTo<T, T>;

template <typename T>
concept CopyConstructible =
  MoveConstructible<T>
  &&  ConstructibleFrom<T,       T&> && ConvertibleTo<      T&, T>
  &&  ConstructibleFrom<T, const T&> && ConvertibleTo<const T&, T>
  &&  ConstructibleFrom<T, const T > && ConvertibleTo<const T, T>;

template <typename T>
concept NoThrowCopyConstructible =
  NoThrowMoveConstructible<T>
  &&  NoThrowConstructibleFrom<T,       T&> && NoThrowConvertibleTo<      T&, T>
  &&  NoThrowConstructibleFrom<T, const T&> && NoThrowConvertibleTo<const T&, T>
  &&  NoThrowConstructibleFrom<T, const T > && NoThrowConvertibleTo<const T, T>;

template <typename T>
concept Swappable = std::swappable<T>;

template <typename T>
concept EqualityComparable = std::equality_comparable<T>;

// T is a type
// X is a Container
// A is an Allocator
// if X::allocator_type then
//   std::same_as<typename X::allocator_type,
//                typename std::allocator_traits<A>::template rebind_alloc<T>>
// otherwise
//   no condition; we use std::allocator<T> regardless of A
//
// see [22.2.1].16
template <typename T, typename X, typename A, typename ...Args>
concept EmplaceConstructible =
  std::same_as<typename X::value_type, T>
  &&  (  (  requires { typename X::allocator_type; } // only perform this check if X is
            &&  std::same_as<typename X::allocator_type, // allocator-aware
            typename std::allocator_traits<A>::template rebind_alloc<T>>
&&  (  requires (A m, T *p, Args&&... args) {
  m.construct (p, std::forward<Args> (args)...);
}
||  requires (T *p, Args&&... args) {
#if __cplusplus >= 202002L // c++20 fully featured
  { std::construct_at (p, std::forward<Args> (args)...) } -> std::same_as<T *>;
#else
  ::new (std::declval<void *> ()) T (std::declval<Args> ()...);
#endif
}))
||  (! requires { typename X::allocator_type; }
&&  requires (T *p, Args&&... args) {
#if __cplusplus >= 202002L // c++20 fully featured
  { std::construct_at (p, std::forward<Args> (args)...) } -> std::same_as<T *>;
#else
  ::new (std::declval<void *> ()) T (std::declval<Args> ()...);
#endif
}));

template <typename T, typename X,
typename A = typename std::conditional<requires { typename X::allocator_type; },
         typename X::allocator_type,
         std::allocator<T>>::type>
         concept DefaultInsertable = EmplaceConstructible<T, X, A>;

template <typename T, typename X,
typename A = typename std::conditional<requires { typename X::allocator_type; },
         typename X::allocator_type,
         std::allocator<T>>::type>
         concept MoveInsertable = EmplaceConstructible<T, X, A, T>;

template <typename T, typename X,
typename A = typename std::conditional<requires { typename X::allocator_type; },
         typename X::allocator_type,
         std::allocator<T>>::type>
         concept CopyInsertable = MoveInsertable<T, X, A>
                                  &&  EmplaceConstructible<T, X, A,       T&>
                                  &&  EmplaceConstructible<T, X, A, const T&>;

// same method as with EmplaceConstructible
template <typename T, typename X,
typename A = typename std::conditional<requires { typename X::allocator_type; },
         typename X::allocator_type,
         std::allocator<T>>::type>
         concept Erasable =
           std::same_as<typename X::value_type, T>
           &&  (  (  requires { typename X::allocator_type; } // if X is allocator aware
                     &&  std::same_as<typename X::allocator_type,
                     typename std::allocator_traits<A>::template rebind_alloc<T>>
&&  (  requires (A m, T *p) {
  m.destroy (p);
}
||   std::is_destructible<T>::value))
||  (! requires { typename X::allocator_type; }
     &&  std::is_destructible<T>::value));

template <typename T>
concept ContextuallyConvertibleToBool = std::constructible_from<bool, T>;

template <typename T>
concept BoolConstant = std::derived_from<T, std::true_type>
|| std::derived_from<T, std::false_type>;

template <typename T>
concept NullablePointer =
EqualityComparable<T>
&&  DefaultConstructible<T>
&&  CopyConstructible<T>
&&  CopyAssignable<T>
&&  Destructible<T>
&&  ConstructibleFrom<T, std::nullptr_t>
&&  ConvertibleTo<std::nullptr_t, T>
&&  requires (T p, T q, std::nullptr_t np) {
  T (np);
  {
    p = np
  }
  -> std::same_as<T&>;
  {
    p  != q
  }
  -> ContextuallyConvertibleToBool;
  {
    p  == np
  }
  -> ContextuallyConvertibleToBool;
  {
    np == p
  }
  -> ContextuallyConvertibleToBool;
  {
    p  != np
  }
  -> ContextuallyConvertibleToBool;
  {
    np != p
  }
  -> ContextuallyConvertibleToBool;
};

static_assert (  NullablePointer<int *>);
static_assert (! NullablePointer<int>);

template <typename A, typename T, typename U = T *>
concept AllocatorFor =
NoThrowCopyConstructible<A>
&&  requires (A a,
              typename std::allocator_traits<A>::template rebind_alloc<U> b,
              U xp,
              typename std::allocator_traits<A>::pointer p,
              typename std::allocator_traits<A>::const_pointer cp,
              typename std::allocator_traits<A>::void_pointer vp,
              typename std::allocator_traits<A>::const_void_pointer cvp,
              typename std::allocator_traits<A>::value_type& r,
              typename std::allocator_traits<A>::size_type n) {
  /** Inner types **/
  // A::pointer
  requires NullablePointer<            typename std::allocator_traits<A>::pointer>;
  requires std::random_access_iterator<typename std::allocator_traits<A>::pointer>;
  requires std::contiguous_iterator<   typename std::allocator_traits<A>::pointer>;

  // A::const_pointer
  requires NullablePointer<            typename std::allocator_traits<A>::const_pointer>;
  requires std::random_access_iterator<typename std::allocator_traits<A>::const_pointer>;
  requires std::contiguous_iterator<   typename std::allocator_traits<A>::const_pointer>;

  requires std::convertible_to<typename std::allocator_traits<A>::pointer,
           typename std::allocator_traits<A>::const_pointer>;

  // A::void_pointer
  requires NullablePointer<typename std::allocator_traits<A>::void_pointer>;

  requires std::convertible_to<typename std::allocator_traits<A>::pointer,
           typename std::allocator_traits<A>::void_pointer>;

  requires std::same_as<typename std::allocator_traits<A>::void_pointer,
           typename std::allocator_traits<decltype (b)>::void_pointer>;

  // A::const_void_pointer
  requires NullablePointer<typename std::allocator_traits<A>::const_void_pointer>;

  requires std::convertible_to<typename std::allocator_traits<A>::pointer,
           typename std::allocator_traits<A>::const_void_pointer>;

  requires std::convertible_to<typename std::allocator_traits<A>::const_pointer,
           typename std::allocator_traits<A>::const_void_pointer>;

  requires std::convertible_to<typename std::allocator_traits<A>::void_pointer,
           typename std::allocator_traits<A>::const_void_pointer>;

  requires std::same_as<typename std::allocator_traits<A>::const_void_pointer,
           typename std::allocator_traits<decltype (b)>::const_void_pointer>;

  // A::value_type
  typename A::value_type;
  requires std::same_as<typename A::value_type, T>;
  requires std::same_as<typename A::value_type,
           typename std::allocator_traits<A>::value_type>;

  // A::size_type
  requires std::unsigned_integral<typename std::allocator_traits<A>::size_type>;

  // A::difference_type
  requires std::signed_integral<typename std::allocator_traits<A>::difference_type>;

  // A::template rebind<U>::other [optional]
  requires ! requires { typename A::template rebind<U>::other; }
  ||  requires {
    requires std::same_as<decltype (b), typename A::template rebind<U>::other>;
    requires std::same_as<A, typename decltype (b)::template rebind<T>::other>;
  };

  /** Operations on pointers **/
  {
    *p
  }
  -> std::same_as<typename A::value_type&>;
  {
    *cp
  }
  -> std::same_as<const typename A::value_type&>;

  // Language in the standard implies that `decltype (p)` must either
  // be a raw pointer or implement `operator->`. There is no mention
  // of `std::to_address` or `std::pointer_traits<Ptr>::to_address`.
  requires std::same_as<decltype (p), typename A::value_type *>
  ||  requires {
    { p.operator-> () } -> std::same_as<typename A::value_type *>;
  };

  requires std::same_as<decltype (cp), const typename A::value_type *>
  ||  requires {
    { cp.operator-> () } -> std::same_as<const typename A::value_type *>;
  };

  {
    static_cast<decltype (p)> (vp)
  }
  -> std::same_as<decltype (p)>;
  {
    static_cast<decltype (cp)> (cvp)
  }
  -> std::same_as<decltype (cp)>;

  {
    std::pointer_traits<decltype (p)>::pointer_to (r)
  }
  -> std::same_as<decltype (p)>;

  /** Storage and lifetime operations **/
  // a.allocate (n)
  {
    a.allocate (n)
  }
  -> std::same_as<decltype (p)>;

  // a.allocate (n, cvp) [optional]
  requires ! requires { a.allocate (n, cvp); }
  ||  requires { { a.allocate (n, cvp) } -> std::same_as<decltype (p)>; };

  // a.deallocate (p, n)
  {
    a.deallocate (p, n)
  }
  -> std::convertible_to<void>;

  // a.max_size () [optional]
  requires ! requires { a.max_size (); }
  ||  requires { { a.max_size () } -> std::same_as<decltype (n)>; };

  // a.construct (xp, args) [optional]
  requires ! requires { a.construct (xp); }
  ||  requires { { a.construct (xp) } -> std::convertible_to<void>; };

  // a.destroy (xp) [optional]
  requires ! requires { a.destroy (xp); }
  ||  requires { { a.destroy (xp) } -> std::convertible_to<void>; };

  /** Relationship between instances **/
  requires NoThrowConstructibleFrom<A, decltype (b)>;
  requires NoThrowConstructibleFrom<A, decltype (std::move (b))>;

  requires BoolConstant<typename std::allocator_traits<A>::is_always_equal>;

  /** Influence on container operations **/
  // a.select_on_container_copy_construction () [optional]
  requires ! requires { a.select_on_container_copy_construction (); }
  ||  requires {
    { a.select_on_container_copy_construction () } -> std::same_as<A>;
  };

  requires BoolConstant<
  typename std::allocator_traits<A>::propagate_on_container_copy_assignment>;

  requires BoolConstant<
  typename std::allocator_traits<A>::propagate_on_container_move_assignment>;

  requires BoolConstant<
  typename std::allocator_traits<A>::propagate_on_container_swap>;

  {
    a == b
  }
  -> std::same_as<bool>;
  {
    a != b
  }
  -> std::same_as<bool>;
}
&&  requires (A a1, A a2) {
  {
    a1 == a2
  }
  -> std::same_as<bool>;
  {
    a1 != a2
  }
  -> std::same_as<bool>;
};

static_assert (AllocatorFor<std::allocator<int>, int>,
               "std::allocator<int> failed to meet Allocator concept requirements.");

template <typename A>
concept Allocator = AllocatorFor<A, typename A::value_type>;

namespace small_vector {

// Basically, these shut off the concepts if we have an incomplete type.
// This namespace is only needed because of issues on Clang
// preventing us from short-circuiting for incomplete types.

template <typename T>
concept Destructible =
! concepts::Complete<T> || concepts::Destructible<T>;

template <typename T>
concept MoveAssignable =
! concepts::Complete<T> || concepts::MoveAssignable<T>;

template <typename T>
concept CopyAssignable =
! concepts::Complete<T> || concepts::CopyAssignable<T>;

template <typename T>
concept MoveConstructible =
! concepts::Complete<T> || concepts::MoveConstructible<T>;

template <typename T>
concept CopyConstructible =
! concepts::Complete<T> || concepts::CopyConstructible<T>;

template <typename T>
concept Swappable =
! concepts::Complete<T> || concepts::Swappable<T>;

template <typename T, typename SmallVector, typename Alloc>
concept DefaultInsertable =
! concepts::Complete<T> || concepts::DefaultInsertable<T, SmallVector, Alloc>;

template <typename T, typename SmallVector, typename Alloc>
concept MoveInsertable =
! concepts::Complete<T> || concepts::MoveInsertable<T, SmallVector, Alloc>;

template <typename T, typename SmallVector, typename Alloc>
concept CopyInsertable =
! concepts::Complete<T> || concepts::CopyInsertable<T, SmallVector, Alloc>;

template <typename T, typename SmallVector, typename Alloc>
concept Erasable =
! concepts::Complete<T> || concepts::Erasable<T, SmallVector, Alloc>;

template <typename T, typename SmallVector, typename Alloc, typename ...Args>
concept EmplaceConstructible =
! concepts::Complete<T> || concepts::EmplaceConstructible<T, SmallVector, Alloc, Args...>;

template <typename Alloc, typename T>
concept AllocatorFor =
! concepts::Complete<T> || concepts::AllocatorFor<Alloc, T>;

template <typename Alloc>
concept Allocator = AllocatorFor<Alloc, typename Alloc::value_type>;

} // namespace gch::concepts::small_vector

} // namespace gch::concepts

#endif

template <typename Allocator>
#ifdef PLUMED_GCH_LIB_CONCEPTS
requires concepts::small_vector::Allocator<Allocator>
#endif
struct default_buffer_size;

template <typename T,
          unsigned InlineCapacity = default_buffer_size<std::allocator<T>>::value,
          typename Allocator      = std::allocator<T>>
#ifdef PLUMED_GCH_LIB_CONCEPTS
requires concepts::small_vector::AllocatorFor<Allocator, T>
#endif
class small_vector;

template <typename Allocator>
#ifdef PLUMED_GCH_LIB_CONCEPTS
requires concepts::small_vector::Allocator<Allocator>
#endif
struct default_buffer_size {
private:
  template <typename, typename Enable = void>
  struct is_complete
: std::false_type {
  };

  template <typename U>
  struct is_complete<U, decltype (static_cast<void> (sizeof (U)))>
: std::true_type {
  };

public:
  using allocator_type     = Allocator;
  using value_type         = typename std::allocator_traits<allocator_type>::value_type;
  using empty_small_vector = small_vector<value_type, 0, allocator_type>;

  static_assert (is_complete<value_type>::value,
                 "Calculation of a default number of elements requires that `T` be complete.");

  static constexpr
  unsigned
  buffer_max = 256;

  static constexpr
  unsigned
  ideal_total = PLUMED_GCH_SMALL_VECTOR_DEFAULT_SIZE;

#ifndef PLUMED_GCH_UNRESTRICTED_DEFAULT_BUFFER_SIZE

  // FIXME: Some compilers will not emit the error from this static_assert
  //        while instantiating a small_vector, and attribute the mistake
  //        to some random other function.
  // static_assert (sizeof (value_type) <= buffer_max, "`sizeof (T)` too large");

#endif

  static constexpr
  unsigned
  ideal_buffer = ideal_total - sizeof (empty_small_vector);

  static_assert (sizeof (empty_small_vector) != 0,
                 "Empty `small_vector` should not have size 0.");

  static_assert (ideal_buffer < ideal_total,
                 "Empty `small_vector` is larger than ideal_total.");

  static constexpr
  unsigned
  value = (sizeof (value_type) <= ideal_buffer) ? (ideal_buffer / sizeof (value_type)) : 1;
};

#ifdef PLUMED_GCH_VARIABLE_TEMPLATES

template <typename Allocator>
PLUMED_GCH_INLINE_VARIABLE constexpr
unsigned
default_buffer_size_v = default_buffer_size<Allocator>::value;

#endif

template <typename Pointer, typename DifferenceType>
class small_vector_iterator {
public:
  using difference_type   = DifferenceType;
  using value_type        = typename std::iterator_traits<Pointer>::value_type;
  using pointer           = typename std::iterator_traits<Pointer>::pointer;
  using reference         = typename std::iterator_traits<Pointer>::reference;
  using iterator_category = typename std::iterator_traits<Pointer>::iterator_category;
#ifdef PLUMED_GCH_LIB_CONCEPTS
  using iterator_concept  = std::contiguous_iterator_tag;
#endif

//  small_vector_iterator            (void)                             = impl;
  small_vector_iterator            (const small_vector_iterator&)     = default;
  small_vector_iterator            (small_vector_iterator&&) noexcept = default;
  small_vector_iterator& operator= (const small_vector_iterator&)     = default;
  small_vector_iterator& operator= (small_vector_iterator&&) noexcept = default;
  ~small_vector_iterator           (void)                             = default;

#ifdef NDEBUG
  small_vector_iterator (void) = default;
#else
  constexpr
  small_vector_iterator (void) noexcept
    : m_ptr ()
  { }
#endif

  constexpr explicit
  small_vector_iterator (const Pointer& p) noexcept
    : m_ptr (p)
  { }

  template <typename U, typename D,
            typename std::enable_if<std::is_convertible<U, Pointer>::value>::type * = nullptr>
  constexpr PLUMED_GCH_IMPLICIT_CONVERSION
  small_vector_iterator (const small_vector_iterator<U, D>& other) noexcept
    : m_ptr (other.base ())
  { }

  PLUMED_GCH_CPP14_CONSTEXPR
  small_vector_iterator&
  operator++ (void) noexcept {
    ++m_ptr;
    return *this;
  }

  PLUMED_GCH_CPP14_CONSTEXPR
  small_vector_iterator
  operator++ (int) noexcept {
    return small_vector_iterator (m_ptr++);
  }

  PLUMED_GCH_CPP14_CONSTEXPR
  small_vector_iterator&
  operator-- (void) noexcept {
    --m_ptr;
    return *this;
  }

  PLUMED_GCH_CPP14_CONSTEXPR
  small_vector_iterator
  operator-- (int) noexcept {
    return small_vector_iterator (m_ptr--);
  }

  PLUMED_GCH_CPP14_CONSTEXPR
  small_vector_iterator&
  operator+= (difference_type n) noexcept {
    m_ptr += n;
    return *this;
  }

  constexpr
  small_vector_iterator
  operator+ (difference_type n) const noexcept {
    return small_vector_iterator (m_ptr + n);
  }

  PLUMED_GCH_CPP14_CONSTEXPR
  small_vector_iterator&
  operator-= (difference_type n) noexcept {
    m_ptr -= n;
    return *this;
  }

  constexpr
  small_vector_iterator
  operator- (difference_type n) const noexcept {
    return small_vector_iterator (m_ptr - n);
  }

  constexpr
  reference
  operator* (void) const noexcept {
#ifdef PLUMED_GCH_LIB_LAUNDER
    return launder_and_dereference (m_ptr);
#else
    return *m_ptr;
#endif
  }

  constexpr
  pointer
  operator-> (void) const noexcept {
    return get_pointer (m_ptr);
  }

  constexpr
  reference
  operator[] (difference_type n) const noexcept {
#ifdef PLUMED_GCH_LIB_LAUNDER
    return launder_and_dereference (m_ptr + n);
#else
    return m_ptr[n];
#endif
  }

  constexpr
  const Pointer&
  base (void) const noexcept {
    return m_ptr;
  }

private:
  template <typename Ptr = Pointer,
            typename std::enable_if<std::is_pointer<Ptr>::value, bool>::type = true>
  static constexpr
  pointer
  get_pointer (Pointer ptr) noexcept {
    return ptr;
  }

  template <typename Ptr = Pointer,
            typename std::enable_if<! std::is_pointer<Ptr>::value, bool>::type = false>
  static constexpr
  pointer
  get_pointer (Pointer ptr) noexcept {
    // Given the requirements for Allocator, Pointer must either be a raw pointer, or
    // have a defined operator-> which returns a raw pointer.
    return ptr.operator-> ();
  }

#ifdef PLUMED_GCH_LIB_LAUNDER

  template <typename Ptr = Pointer,
            typename std::enable_if<std::is_pointer<Ptr>::value, bool>::type = true>
  static constexpr
  reference
  launder_and_dereference (Pointer ptr) noexcept {
    return *std::launder (ptr);
  }

  template <typename Ptr = Pointer,
            typename std::enable_if<! std::is_pointer<Ptr>::value, bool>::type = false>
  static constexpr
  reference
  launder_and_dereference (Pointer ptr) noexcept {
    return *ptr;
  }

#endif

  Pointer m_ptr;
};

#ifdef PLUMED_GCH_LIB_THREE_WAY_COMPARISON

template <typename PointerLHS, typename DifferenceTypeLHS,
          typename PointerRHS, typename DifferenceTypeRHS>
constexpr
bool
operator== (const small_vector_iterator<PointerLHS, DifferenceTypeLHS>& lhs,
            const small_vector_iterator<PointerRHS, DifferenceTypeRHS>& rhs)
noexcept (noexcept (lhs.base () == rhs.base ()))
requires requires { { lhs.base () == rhs.base () } -> std::convertible_to<bool>; } {
  return lhs.base () == rhs.base ();
}

template <typename Pointer, typename DifferenceType>
constexpr
bool
operator== (const small_vector_iterator<Pointer, DifferenceType>& lhs,
            const small_vector_iterator<Pointer, DifferenceType>& rhs)
noexcept (noexcept (lhs.base () == rhs.base ()))
requires requires { { lhs.base () == rhs.base () } -> std::convertible_to<bool>; } {
  return lhs.base () == rhs.base ();
}

template <typename PointerLHS, typename DifferenceTypeLHS,
          typename PointerRHS, typename DifferenceTypeRHS>
requires std::three_way_comparable_with<PointerLHS, PointerRHS>
constexpr
auto
operator<=> (const small_vector_iterator<PointerLHS, DifferenceTypeLHS>& lhs,
             const small_vector_iterator<PointerRHS, DifferenceTypeRHS>& rhs)
noexcept (noexcept (lhs.base () <=> rhs.base ())) {
  return lhs.base () <=> rhs.base ();
}

template <typename Pointer, typename DifferenceType>
requires std::three_way_comparable<Pointer>
constexpr
auto
operator<=> (const small_vector_iterator<Pointer, DifferenceType>& lhs,
             const small_vector_iterator<Pointer, DifferenceType>& rhs)
noexcept (noexcept (lhs.base () <=> rhs.base ())) {
  return lhs.base () <=> rhs.base ();
}

template <typename PointerLHS, typename DifferenceTypeLHS,
          typename PointerRHS, typename DifferenceTypeRHS>
constexpr
auto
operator<=> (const small_vector_iterator<PointerLHS, DifferenceTypeLHS>& lhs,
             const small_vector_iterator<PointerRHS, DifferenceTypeRHS>& rhs)
noexcept (noexcept (lhs.base () < rhs.base ()) && noexcept (rhs.base () < lhs.base ())) {
  using ordering = std::weak_ordering;
  return (lhs.base () < rhs.base ()) ? ordering::less
         : (rhs.base () < lhs.base ()) ? ordering::greater
         : ordering::equivalent;
}

template <typename Pointer, typename DifferenceType>
constexpr
auto
operator<=> (const small_vector_iterator<Pointer, DifferenceType>& lhs,
             const small_vector_iterator<Pointer, DifferenceType>& rhs)
noexcept (noexcept (lhs.base () < rhs.base ()) && noexcept (rhs.base () < lhs.base ())) {
  using ordering = std::weak_ordering;
  return (lhs.base () < rhs.base ()) ? ordering::less
         : (rhs.base () < lhs.base ()) ? ordering::greater
         : ordering::equivalent;
}

#else

// Note: Passing this on from "Gaby" in stl_iterator.h -- templated
//       comparisons in generic code should have overloads for both
//       homogenous and heterogeneous types. This is because we get
//       ambiguous overload resolution when std::rel_ops is visible
//       (ie. `using namespace std::rel_ops`).

template <typename PointerLHS, typename DifferenceTypeLHS,
          typename PointerRHS, typename DifferenceTypeRHS>
constexpr
bool
operator== (const small_vector_iterator<PointerLHS, DifferenceTypeLHS>& lhs,
            const small_vector_iterator<PointerRHS, DifferenceTypeRHS>& rhs) noexcept {
  return lhs.base () == rhs.base ();
}

template <typename Pointer, typename DifferenceType>
constexpr
bool
operator== (const small_vector_iterator<Pointer, DifferenceType>& lhs,
            const small_vector_iterator<Pointer, DifferenceType>& rhs) noexcept {
  return lhs.base () == rhs.base ();
}

template <typename PointerLHS, typename DifferenceTypeLHS,
          typename PointerRHS, typename DifferenceTypeRHS>
constexpr
bool
operator!= (const small_vector_iterator<PointerLHS, DifferenceTypeLHS>& lhs,
            const small_vector_iterator<PointerRHS, DifferenceTypeRHS>& rhs) noexcept {
  return lhs.base () != rhs.base ();
}

template <typename Pointer, typename DifferenceType>
constexpr
bool
operator!= (const small_vector_iterator<Pointer, DifferenceType>& lhs,
            const small_vector_iterator<Pointer, DifferenceType>& rhs) noexcept {
  return lhs.base () != rhs.base ();
}

template <typename PointerLHS, typename DifferenceTypeLHS,
          typename PointerRHS, typename DifferenceTypeRHS>
constexpr
bool
operator< (const small_vector_iterator<PointerLHS, DifferenceTypeLHS>& lhs,
const small_vector_iterator<PointerRHS, DifferenceTypeRHS>& rhs) noexcept {
  return lhs.base () < rhs.base ();
}

template <typename Pointer, typename DifferenceType>
constexpr
bool
operator< (const small_vector_iterator<Pointer, DifferenceType>& lhs,
const small_vector_iterator<Pointer, DifferenceType>& rhs) noexcept {
  return lhs.base () < rhs.base ();
}

template <typename PointerLHS, typename DifferenceTypeLHS,
          typename PointerRHS, typename DifferenceTypeRHS>
constexpr
bool
operator>= (const small_vector_iterator<PointerLHS, DifferenceTypeLHS>& lhs,
            const small_vector_iterator<PointerRHS, DifferenceTypeRHS>& rhs) noexcept {
  return lhs.base () >= rhs.base ();
}

template <typename Pointer, typename DifferenceType>
constexpr
bool
operator>= (const small_vector_iterator<Pointer, DifferenceType>& lhs,
            const small_vector_iterator<Pointer, DifferenceType>& rhs) noexcept {
  return lhs.base () >= rhs.base ();
}

template <typename PointerLHS, typename DifferenceTypeLHS,
          typename PointerRHS, typename DifferenceTypeRHS>
constexpr
bool
operator> (const small_vector_iterator<PointerLHS, DifferenceTypeLHS>& lhs,
           const small_vector_iterator<PointerRHS, DifferenceTypeRHS>& rhs) noexcept {
  return lhs.base () > rhs.base ();
}

template <typename Pointer, typename DifferenceType>
constexpr
bool
operator> (const small_vector_iterator<Pointer, DifferenceType>& lhs,
           const small_vector_iterator<Pointer, DifferenceType>& rhs) noexcept {
  return lhs.base () > rhs.base ();
}

template <typename PointerLHS, typename DifferenceTypeLHS,
          typename PointerRHS, typename DifferenceTypeRHS>
constexpr
bool
operator<= (const small_vector_iterator<PointerLHS, DifferenceTypeLHS>& lhs,
const small_vector_iterator<PointerRHS, DifferenceTypeRHS>& rhs) noexcept {
  return lhs.base () <= rhs.base ();
}

template <typename Pointer, typename DifferenceType>
constexpr
bool
operator<= (const small_vector_iterator<Pointer, DifferenceType>& lhs,
const small_vector_iterator<Pointer, DifferenceType>& rhs) noexcept {
  return lhs.base () <= rhs.base ();
}

#endif

template <typename PointerLHS, typename PointerRHS, typename DifferenceType>
constexpr
DifferenceType
operator- (const small_vector_iterator<PointerLHS, DifferenceType>& lhs,
           const small_vector_iterator<PointerRHS, DifferenceType>& rhs) noexcept {
  return static_cast<DifferenceType> (lhs.base () - rhs.base ());
}

template <typename Pointer, typename DifferenceType>
constexpr
DifferenceType
operator- (const small_vector_iterator<Pointer, DifferenceType>& lhs,
           const small_vector_iterator<Pointer, DifferenceType>& rhs) noexcept {
  return static_cast<DifferenceType> (lhs.base () - rhs.base ());
}

template <typename Pointer, typename DifferenceType>
constexpr
small_vector_iterator<Pointer, DifferenceType>
operator+ (DifferenceType n, const small_vector_iterator<Pointer, DifferenceType>& it) noexcept {
  return it + n;
}

namespace detail {

#ifndef PLUMED_GCH_LIB_IS_SWAPPABLE

namespace small_vector_adl {

using std::swap;

template <typename T, typename Enable = void>
struct is_nothrow_swappable
: std::false_type {
};

template <typename T>
struct is_nothrow_swappable<T, decltype (swap (std::declval<T&> (), std::declval<T&> ()))>
: std::integral_constant<bool, noexcept (swap (std::declval<T&> (), std::declval<T&> ()))> {
};

}

#endif

template <typename T, unsigned InlineCapacity>
class inline_storage {
public:
  using value_ty = T;

  inline_storage            (void)                      = default;
  inline_storage            (const inline_storage&)     = delete;
  inline_storage            (inline_storage&&) noexcept = delete;
  inline_storage& operator= (const inline_storage&)     = delete;
  inline_storage& operator= (inline_storage&&) noexcept = delete;
  ~inline_storage           (void)                      = default;

  PLUMED_GCH_NODISCARD PLUMED_GCH_CPP14_CONSTEXPR
  value_ty *
  get_inline_ptr (void) noexcept {
    return static_cast<value_ty *> (static_cast<void *> (std::addressof (*m_data)));
  }

  PLUMED_GCH_NODISCARD constexpr
  const value_ty *
  get_inline_ptr (void) const noexcept {
    return static_cast<const value_ty *> (static_cast<const void *> (std::addressof (*m_data)));
  }

  static constexpr
  std::size_t
  element_size (void) noexcept {
    return sizeof (value_ty);
  }

  static constexpr
  std::size_t
  alignment (void) noexcept {
    return alignof (value_ty);
  }

  static constexpr
  unsigned
  num_elements (void) noexcept {
    return InlineCapacity;
  }

  static constexpr
  std::size_t
  num_bytes (void) noexcept {
    return num_elements () * element_size ();
  }

private:
  typename std::aligned_storage<element_size (), alignment ()>::type m_data[num_elements ()];
};

template <typename T>
class PLUMED_GCH_EMPTY_BASE inline_storage<T, 0> {
public:
  using value_ty = T;

  inline_storage            (void)                      = default;
  inline_storage            (const inline_storage&)     = delete;
  inline_storage            (inline_storage&&) noexcept = delete;
  inline_storage& operator= (const inline_storage&)     = delete;
  inline_storage& operator= (inline_storage&&) noexcept = delete;
  ~inline_storage           (void)                      = default;

  PLUMED_GCH_NODISCARD PLUMED_GCH_CPP14_CONSTEXPR
  value_ty *
  get_inline_ptr (void) noexcept {
    return nullptr;
  }

  PLUMED_GCH_NODISCARD constexpr
  const value_ty *
  get_inline_ptr (void) const noexcept {
    return nullptr;
  }

  static constexpr
  std::size_t
  element_size (void) noexcept {
    return sizeof (value_ty);
  }

  static constexpr
  std::size_t
  alignment (void) noexcept {
    return alignof (value_ty);
  }

  static constexpr
  unsigned
  num_elements (void) noexcept {
    return 0;
  }

  static constexpr
  std::size_t
  num_bytes (void) noexcept {
    return 0;
  }
};

template <typename Allocator, bool AvailableForEBO = std::is_empty<Allocator>::value
#ifdef PLUMED_GCH_LIB_IS_FINAL
          &&! std::is_final<Allocator>::value
#endif // If you are using this with C++11 just don't use an allocator marked as final :P
          >
class allocator_inliner;

template <typename Allocator>
class PLUMED_GCH_EMPTY_BASE allocator_inliner<Allocator, true>
: private Allocator {
  using alloc_traits = std::allocator_traits<Allocator>;

  static constexpr
  bool
  copy_assign_is_noop = ! alloc_traits::propagate_on_container_copy_assignment::value;

  static constexpr
  bool
  move_assign_is_noop = ! alloc_traits::propagate_on_container_move_assignment::value;

  static constexpr
  bool
  swap_is_noop = ! alloc_traits::propagate_on_container_swap::value;

  template <bool IsNoOp = copy_assign_is_noop,
            typename std::enable_if<IsNoOp, bool>::type = true>
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  maybe_assign (const allocator_inliner&) noexcept { }

  template <bool IsNoOp = copy_assign_is_noop,
            typename std::enable_if<! IsNoOp, bool>::type = false>
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  maybe_assign (const allocator_inliner& other)
  noexcept (noexcept (std::declval<Allocator&> ().operator= (other))) {
    Allocator::operator= (other);
  }

  template <bool IsNoOp = move_assign_is_noop,
            typename std::enable_if<IsNoOp, bool>::type = true>
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  maybe_assign (allocator_inliner&&) noexcept { }

  template <bool IsNoOp = move_assign_is_noop,
            typename std::enable_if<! IsNoOp, bool>::type = false>
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  maybe_assign (allocator_inliner&& other)
  noexcept (noexcept (std::declval<Allocator&> ().operator= (std::move (other)))) {
    Allocator::operator= (std::move (other));
  }

public:
  allocator_inliner            (void)                         = default;
  allocator_inliner            (const allocator_inliner&)     = default;
  allocator_inliner            (allocator_inliner&&) noexcept = default;
//    allocator_inliner& operator= (const allocator_inliner&)     = impl;
//    allocator_inliner& operator= (allocator_inliner&&) noexcept = impl;
  ~allocator_inliner           (void)                         = default;

  constexpr explicit
  allocator_inliner (const Allocator& alloc) noexcept
    : Allocator (alloc)
  { }

  PLUMED_GCH_CPP20_CONSTEXPR
  allocator_inliner&
  operator= (const allocator_inliner& other)
  noexcept (noexcept (std::declval<allocator_inliner&> ().maybe_assign (other))) {
    assert (&other != this
            &&  "`allocator_inliner` should not participate in self-copy-assignment.");
    maybe_assign (other);
    return *this;
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  allocator_inliner&
  operator= (allocator_inliner&& other)
  noexcept (noexcept (std::declval<allocator_inliner&> ().maybe_assign (std::move (other)))) {
    assert (&other != this
            &&  "`allocator_inliner` should not participate in self-move-assignment.");
    maybe_assign (std::move (other));
    return *this;
  }

  PLUMED_GCH_CPP14_CONSTEXPR
  Allocator&
  allocator_ref (void) noexcept {
    return *this;
  }

  constexpr
  const Allocator&
  allocator_ref (void) const noexcept {
    return *this;
  }

  template <bool IsNoOp = swap_is_noop,
            typename std::enable_if<IsNoOp, bool>::type = true>
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  swap (allocator_inliner&)
  { }

  template <bool IsNoOp = swap_is_noop,
            typename std::enable_if<! IsNoOp, bool>::type = false>
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  swap (allocator_inliner& other) {
    using std::swap;
    swap (static_cast<Allocator&> (*this), static_cast<Allocator&> (other));
  }
};

template <typename Allocator>
class allocator_inliner<Allocator, false> {
  using alloc_traits = std::allocator_traits<Allocator>;

  static constexpr
  bool
  copy_assign_is_noop = ! alloc_traits::propagate_on_container_copy_assignment::value;

  static constexpr
  bool
  move_assign_is_noop = ! alloc_traits::propagate_on_container_move_assignment::value;

  static constexpr
  bool
  swap_is_noop = ! alloc_traits::propagate_on_container_swap::value;

  template <bool IsNoOp = copy_assign_is_noop,
            typename std::enable_if<IsNoOp, bool>::type = true>
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  maybe_assign (const allocator_inliner&) noexcept { }

  template <bool IsNoOp = copy_assign_is_noop,
            typename std::enable_if<! IsNoOp, bool>::type = false>
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  maybe_assign (const allocator_inliner& other)
  noexcept (noexcept (std::declval<decltype (other.m_alloc)&> () = other.m_alloc)) {
    m_alloc = other.m_alloc;
  }

  template <bool IsNoOp = move_assign_is_noop,
            typename std::enable_if<IsNoOp, bool>::type = true>
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  maybe_assign (allocator_inliner&&) noexcept { }

  template <bool IsNoOp = move_assign_is_noop,
            typename std::enable_if<! IsNoOp, bool>::type = false>
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  maybe_assign (allocator_inliner&& other)
  noexcept (noexcept (std::declval<decltype (other.m_alloc)&> () = std::move (other.m_alloc))) {
    m_alloc = std::move (other.m_alloc);
  }

public:
  allocator_inliner            (void)                         = default;
  allocator_inliner            (const allocator_inliner&)     = default;
  allocator_inliner            (allocator_inliner&&) noexcept = default;
//    allocator_inliner& operator= (const allocator_inliner&)     = impl;
//    allocator_inliner& operator= (allocator_inliner&&) noexcept = impl;
  ~allocator_inliner           (void)                         = default;

  PLUMED_GCH_CPP20_CONSTEXPR explicit
  allocator_inliner (const Allocator& alloc) noexcept
    : m_alloc (alloc)
  { }

  PLUMED_GCH_CPP20_CONSTEXPR
  allocator_inliner&
  operator= (const allocator_inliner& other)
  noexcept (noexcept (std::declval<allocator_inliner&> ().maybe_assign (other))) {
    assert (&other != this
            &&  "`allocator_inliner` should not participate in self-copy-assignment.");
    maybe_assign (other);
    return *this;
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  allocator_inliner&
  operator= (allocator_inliner&& other)
  noexcept (noexcept (std::declval<allocator_inliner&> ().maybe_assign (std::move (other)))) {
    assert (&other != this
            &&  "`allocator_inliner` should not participate in self-move-assignment.");
    maybe_assign (std::move (other));
    return *this;
  }

  PLUMED_GCH_CPP14_CONSTEXPR
  Allocator&
  allocator_ref (void) noexcept {
    return m_alloc;
  }

  constexpr
  const Allocator&
  allocator_ref (void) const noexcept {
    return m_alloc;
  }

  template <bool IsNoOp = swap_is_noop,
            typename std::enable_if<IsNoOp, bool>::type = true>
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  swap (allocator_inliner&)
  { }

  template <bool IsNoOp = swap_is_noop,
            typename std::enable_if<! IsNoOp, bool>::type = false>
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  swap (allocator_inliner& other) {
    using std::swap;
    swap (m_alloc, other.m_alloc);
  }

private:
  Allocator m_alloc;
};

template <typename Allocator>
class PLUMED_GCH_EMPTY_BASE allocator_interface
: public allocator_inliner<Allocator> {
public:
  template <typename, typename = void>
  struct is_complete
: std::false_type {
  };

  template <typename U>
  struct is_complete<U, decltype (static_cast<void> (sizeof (U)))>
: std::true_type {
  };

  using size_type = typename std::allocator_traits<Allocator>::size_type;

  // If difference_type is larger than size_type then we need
  // to rectify that problem.
  using difference_type = typename std::conditional<
                          (
                            static_cast<std::size_t> ((std::numeric_limits<size_type>::max) ())
                            < // less-than
                            static_cast<std::size_t> ((std::numeric_limits<
                                typename std::allocator_traits<Allocator>::difference_type>::max) ())
                          ),
                          typename std::make_signed<size_type>::type,
                          typename std::allocator_traits<Allocator>::difference_type>::type;

private:
  using alloc_base = allocator_inliner<Allocator>;

protected:
  using alloc_ty     = Allocator;
  using alloc_traits = std::allocator_traits<alloc_ty>;
  using value_ty     = typename alloc_traits::value_type;
  using ptr          = typename alloc_traits::pointer;
  using cptr         = typename alloc_traits::const_pointer;
  using vptr         = typename alloc_traits::void_pointer;
  using cvptr        = typename alloc_traits::const_void_pointer;

  // Select the fastest types larger than the user-facing types. These are only intended for
  // internal computations, and should not have any memory footprint visible to consumers.
  using size_ty =
    typename std::conditional<
    (sizeof (size_type) <= sizeof (std::uint8_t)),
    std::uint_fast8_t,
    typename std::conditional<
    (sizeof (size_type) <= sizeof (std::uint16_t)),
    std::uint_fast16_t,
    typename std::conditional<
    (sizeof (size_type) <= sizeof (std::uint32_t)),
    std::uint_fast32_t,
    typename std::conditional<
    (sizeof (size_type) <= sizeof (std::uint64_t)),
    std::uint_fast64_t,
    size_type
    >::type
    >::type
    >::type
    >::type;

  using diff_ty =
    typename std::conditional<
    (sizeof (difference_type) <= sizeof (std::int8_t)),
    std::int_fast8_t,
    typename std::conditional<
    (sizeof (difference_type) <= sizeof (std::int16_t)),
    std::int_fast16_t,
    typename std::conditional<
    (sizeof (difference_type) <= sizeof (std::int32_t)),
    std::int_fast32_t,
    typename std::conditional<
    (sizeof (difference_type) <= sizeof (std::int64_t)),
    std::int_fast64_t,
    difference_type
    >::type
    >::type
    >::type
    >::type;

  using alloc_base::allocator_ref;

private:
  template <typename ...>
  using void_t = void;

  template <bool B>
  using bool_constant = std::integral_constant<bool, B>;

  template <typename V, typename Enable = void>
  struct is_trivially_destructible
: std::false_type {
  };

  template <typename V>
  struct is_trivially_destructible<V, typename std::enable_if<is_complete<V>::value>::type>
: std::is_trivially_destructible<V> {
  };

  template <typename Void, typename T, typename ...Args>
  struct is_trivially_constructible_impl
: std::false_type {
  };

  template <typename V, typename ...Args>
  struct is_trivially_constructible_impl<
    typename std::enable_if<is_complete<V>::value>::type,
    V, Args...>
: std::is_trivially_constructible<V, Args...> {
  };

  template <typename V, typename ...Args>
  struct is_trivially_constructible
: is_trivially_constructible_impl<void, V, Args...> {
  };

  template <typename T, typename Enable = void>
  struct underlying_if_enum {
    using type = T;
  };

  template <typename T>
  struct underlying_if_enum<T, typename std::enable_if<std::is_enum<T>::value>::type>
: std::underlying_type<T> {
  };

  template <typename T>
  using underlying_if_enum_t = typename underlying_if_enum<T>::type;

  template <typename, typename = void>
  struct has_ptr_traits_to_address
: std::false_type {
  };

  template <typename P>
  struct has_ptr_traits_to_address<P,
                                   void_t<decltype (std::pointer_traits<P>::to_address (std::declval<P> ()))>>
: std::true_type {
  };

  template <typename Void, typename A, typename V, typename ...Args>
  struct has_alloc_construct_check
: std::false_type {
  };

  template <typename A, typename V, typename ...Args>
  struct has_alloc_construct_check<
    void_t<decltype (std::declval<A&> ().construct (std::declval<V *> (),
                     std::declval<Args> ()...))>,
    A, V, Args...>
: std::true_type {
  };

  template <typename Void, typename A, typename V, typename ...Args>
  struct has_alloc_construct_impl
: std::false_type {
  };

  template <typename A, typename V, typename ...Args>
  struct has_alloc_construct_impl<typename std::enable_if<is_complete<V>::value>::type,
                                  A, V, Args...>
: has_alloc_construct_check<void, A, V, Args...> {
  };

  template <typename A, typename V, typename ...Args>
  struct has_alloc_construct
: has_alloc_construct_impl<void, A, V, Args...> {
  };

  template <typename A, typename V, typename ...Args>
  struct must_use_alloc_construct
: bool_constant<! std::is_same<A, std::allocator<V>>::value
                  &&  has_alloc_construct<A, V, Args...>::value> {
                  };

  template <typename Void, typename A, typename V>
  struct has_alloc_destroy_impl
: std::false_type {
  };

  template <typename A, typename V>
  struct has_alloc_destroy_impl<
    void_t<decltype (std::declval<A&> ().destroy (std::declval<V *> ()))>,
    A, V>
: std::true_type {
  };

  template <typename A, typename V, typename Enable = void>
  struct has_alloc_destroy
: std::false_type {
  };

  template <typename A, typename V>
  struct has_alloc_destroy<A, V, typename std::enable_if<is_complete<V>::value>::type>
: has_alloc_destroy_impl<void, A, V> {
  };

  template <typename A, typename V>
  struct must_use_alloc_destroy
: bool_constant<! std::is_same<A, std::allocator<V>>::value
                  &&  has_alloc_destroy<A, V>::value> {
                  };

public:
  allocator_interface (void)                           = default;
//    allocator_interface (const allocator_interface&)     = impl;
  allocator_interface (allocator_interface&&) noexcept = default;

  PLUMED_GCH_CPP20_CONSTEXPR
  allocator_interface&
  operator= (const allocator_interface&) = default;

  PLUMED_GCH_CPP20_CONSTEXPR
  allocator_interface&
  operator= (allocator_interface&&) noexcept = default;

  ~allocator_interface (void) = default;

  PLUMED_GCH_CPP20_CONSTEXPR
  allocator_interface (const allocator_interface& other) noexcept
    : alloc_base (alloc_traits::select_on_container_copy_construction (other.allocator_ref ()))
  { }

  constexpr explicit
  allocator_interface (const alloc_ty& alloc) noexcept
    : alloc_base (alloc)
  { }

  template <typename T>
  constexpr explicit
  allocator_interface (T&&, const alloc_ty& alloc) noexcept
    : allocator_interface (alloc)
  { }

  template <typename, typename, typename = void>
  struct is_memcpyable_integral
: std::false_type {
  };

  template <typename From, typename To>
  struct is_memcpyable_integral<From, To,
                                typename std::enable_if<is_complete<From>::value>::type> {
                                  using from = underlying_if_enum_t<From>;
                                  using to   = underlying_if_enum_t<To>;

                                  static constexpr
                                  bool
                                  value = (sizeof (from) == sizeof (to))
                                  &&  (std::is_same<bool, from>::value == std::is_same<bool, to>::value)
                                  &&  std::is_integral<from>::value
                                  &&  std::is_integral<to>::value;
                                                      };

  template <typename From, typename To>
  struct is_convertible_pointer
: bool_constant<std::is_pointer<From>::value
                  &&  std::is_pointer<To>::value
                  &&  std::is_convertible<From, To>::value> {
                  };

  // Memcpyable assignment.
  template <typename QualifiedFrom, typename QualifiedTo = value_ty, typename Enable = void>
  struct is_memcpyable
: std::false_type {
  };

  template <typename QualifiedFrom, typename QualifiedTo>
  struct is_memcpyable<QualifiedFrom, QualifiedTo,
                       typename std::enable_if<is_complete<QualifiedFrom>::value>::type> {
                         static_assert (! std::is_reference<QualifiedTo>::value,
                                        "QualifiedTo must not be a reference.");

                         using from = typename std::remove_reference<
                           typename std::remove_cv<QualifiedFrom>::type>::type;

                         using to = typename std::remove_cv<QualifiedTo>::type;

                         static constexpr
                         bool
                         value = std::is_trivially_assignable<QualifiedTo&, QualifiedFrom>::value
                         &&  std::is_trivially_copyable<to>::value
                         &&  (  std::is_same<typename std::remove_cv<from>::type, to>::value
                                ||  is_memcpyable_integral<from, to>::value
                                ||  is_convertible_pointer<from, to>::value);
                                                          };

  // Memcpyable construction.
  template <typename QualifiedFrom, typename QualifiedTo>
  struct is_uninitialized_memcpyable_impl {
    static_assert (! std::is_reference<QualifiedTo>::value,
                   "QualifiedTo must not be a reference.");

    using from = typename std::remove_reference<
      typename std::remove_cv<QualifiedFrom>::type>::type;

    using to = typename std::remove_cv<QualifiedTo>::type;

    static constexpr
    bool
    value = std::is_trivially_constructible<QualifiedTo, QualifiedFrom>::value
    &&  std::is_trivially_copyable<to>::value
    &&  (  std::is_same<typename std::remove_cv<from>::type, to>::value
           ||  is_memcpyable_integral<from, to>::value
           ||  is_convertible_pointer<from, to>::value)
    &&  (! must_use_alloc_construct<alloc_ty, value_ty, from>::value
         &&! must_use_alloc_destroy<alloc_ty, value_ty>::value);
                                   };

  template <typename To, typename ...Args>
  struct is_uninitialized_memcpyable
: std::false_type {
  };

  template <typename To, typename From>
  struct is_uninitialized_memcpyable<To, From>
: is_uninitialized_memcpyable_impl<From, To> {
  };

  template <typename Iterator>
  struct is_small_vector_iterator
: std::false_type {
  };

  template <typename ...Ts>
  struct is_small_vector_iterator<small_vector_iterator<Ts...>>
: std::true_type {
  };

  template <typename InputIt>
  struct is_contiguous_iterator
: bool_constant<
    std::is_same<InputIt, ptr>::value
    ||  std::is_same<InputIt, cptr>::value
    ||  is_small_vector_iterator<InputIt>::value
#ifdef PLUMED_GCH_LIB_CONCEPTS
    ||  std::contiguous_iterator<InputIt>
#endif
#ifdef PLUMED_GCH_STDLIB_INTEROP
    ||  std::is_same<InputIt, typename std::array<value_ty>::iterator>::value
    ||  std::is_same<InputIt, typename std::array<value_ty>::const_iterator>::value
    ||  (! std::is_same<value_ty, bool>
         &&  (  std::is_same<InputIt, typename std::vector<value_ty>::iterator>::value
                ||  std::is_same<InputIt, typename std::vector<value_ty>::const_iterator>::value)
                                )
    ||  std::is_same<InputIt,
                     decltype (std::begin (std::declval<std::valarray<value_ty>&> ()))>::value
    ||  std::is_same<InputIt,
                     decltype (std::begin (std::declval<const std::valarray<value_ty>&> ()))>::value
#endif
    > {
    };

  template <typename InputIt>
  struct is_memcpyable_iterator
: bool_constant<is_memcpyable<decltype (*std::declval<InputIt> ())>::value
                  &&  is_contiguous_iterator<InputIt>::value> {
                  };

  // Unwrap `move_iterator`s.
  template <typename InputIt>
  struct is_memcpyable_iterator<std::move_iterator<InputIt>>
: is_memcpyable_iterator<InputIt> {
  };

  template <typename InputIt, typename V = value_ty>
  struct is_uninitialized_memcpyable_iterator
: bool_constant<is_uninitialized_memcpyable<V, decltype (*std::declval<InputIt> ())>::value
                  &&  is_contiguous_iterator<InputIt>::value> {
                  };

  // unwrap move_iterators
  template <typename U, typename V>
  struct is_uninitialized_memcpyable_iterator<std::move_iterator<U>, V>
: is_uninitialized_memcpyable_iterator<U, V> {
  };

  PLUMED_GCH_NORETURN
  static PLUMED_GCH_CPP20_CONSTEXPR
  void
  throw_range_length_error (void) {
#ifdef PLUMED_GCH_EXCEPTIONS
    throw std::length_error ("The specified range is too long.");
#else
    std::fprintf (stderr, "[gch::small_vector] The specified range is too long.");
    std::abort ();
#endif
  }

  static constexpr
  value_ty *
  to_address (value_ty *p) noexcept {
    static_assert (! std::is_function<value_ty>::value, "value_ty is a function pointer.");
    return p;
  }

  static constexpr
  const value_ty *
  to_address (const value_ty *p) noexcept {
    static_assert (! std::is_function<value_ty>::value, "value_ty is a function pointer.");
    return p;
  }

  template <typename Pointer,
            typename std::enable_if<has_ptr_traits_to_address<Pointer>::value>::type * = nullptr>
  static constexpr
  auto
  to_address (const Pointer& p) noexcept
  -> decltype (std::pointer_traits<Pointer>::to_address (p)) {
    return std::pointer_traits<Pointer>::to_address (p);
  }

  template <typename Pointer,
            typename std::enable_if<! has_ptr_traits_to_address<Pointer>::value>::type * = nullptr>
  static constexpr
  auto
  to_address (const Pointer& p) noexcept
  -> decltype (to_address (p.operator-> ())) {
    return to_address (p.operator-> ());
  }

  template <typename Integer>
  PLUMED_GCH_NODISCARD
  static PLUMED_GCH_CONSTEVAL
  std::size_t
  numeric_max (void) noexcept {
    static_assert (0 <= (std::numeric_limits<Integer>::max) (), "Integer is nonpositive.");
    return static_cast<std::size_t> ((std::numeric_limits<Integer>::max) ());
  }

  PLUMED_GCH_NODISCARD
  static PLUMED_GCH_CPP17_CONSTEXPR
  size_ty
  internal_range_length (cptr first, cptr last) noexcept {
    // This is guaranteed to be less than or equal to max size_ty.
    return static_cast<size_ty> (last - first);
  }

  template <typename RandomIt>
  PLUMED_GCH_NODISCARD
  static PLUMED_GCH_CPP17_CONSTEXPR
  size_ty
  external_range_length_impl (RandomIt first, RandomIt last, std::random_access_iterator_tag) {
    assert (0 <= (last - first) && "Invalid range.");
    const auto len = static_cast<std::size_t> (last - first);
#ifndef NDEBUG
    if (numeric_max<size_ty> () < len) {
      throw_range_length_error ();
    }
#endif
    return static_cast<size_ty> (len);
  }

  template <typename ForwardIt>
  PLUMED_GCH_NODISCARD
  static PLUMED_GCH_CPP17_CONSTEXPR
  size_ty
  external_range_length_impl (ForwardIt first, ForwardIt last, std::forward_iterator_tag) {
#ifdef PLUMED_GCH_LIB_IS_CONSTANT_EVALUATED
    if (std::is_constant_evaluated ()) {
      // Make sure constexpr doesn't get broken by `using namespace std::rel_ops`.
      typename std::iterator_traits<ForwardIt>::difference_type len = 0;
      for (; ! (first == last); ++first) {
        ++len;
      }
      assert (static_cast<std::size_t> (len) <= numeric_max<size_ty> ());
      return static_cast<size_ty> (len);
    }
#endif

    const auto len = static_cast<std::size_t> (std::distance (first, last));
#ifndef NDEBUG
    if (numeric_max<size_ty> () < len) {
      throw_range_length_error ();
    }
#endif
    return static_cast<size_ty> (len);
  }

  template <typename ForwardIt,
            typename ItDiffT = typename std::iterator_traits<ForwardIt>::difference_type,
            typename std::enable_if<(numeric_max<size_ty> () < numeric_max<ItDiffT> ()),
                                     bool>::type = true>
            PLUMED_GCH_NODISCARD
            static PLUMED_GCH_CPP17_CONSTEXPR
            size_ty
  external_range_length (ForwardIt first, ForwardIt last) {
    using iterator_cat = typename std::iterator_traits<ForwardIt>::iterator_category;
    return external_range_length_impl (first, last, iterator_cat { });
  }

  template <typename ForwardIt,
            typename ItDiffT = typename std::iterator_traits<ForwardIt>::difference_type,
            typename std::enable_if<! (numeric_max<size_ty> () < numeric_max<ItDiffT> ()),
                                       bool>::type = false>
            PLUMED_GCH_NODISCARD
            static PLUMED_GCH_CPP17_CONSTEXPR
            size_ty
  external_range_length (ForwardIt first, ForwardIt last) noexcept {
#ifdef PLUMED_GCH_LIB_IS_CONSTANT_EVALUATED
    if (std::is_constant_evaluated ()) {
      // Make sure constexpr doesn't get broken by `using namespace std::rel_ops`.
      size_ty len = 0;
      for (; ! (first == last); ++first) {
        ++len;
      }
      return len;
    }
#endif

    return static_cast<size_ty> (std::distance (first, last));
  }

  template <typename Iterator,
            typename IteratorDiffT = typename std::iterator_traits<Iterator>::difference_type,
            typename Integer = IteratorDiffT>
  PLUMED_GCH_NODISCARD
  static PLUMED_GCH_CPP17_CONSTEXPR
  Iterator
  unchecked_next (Iterator pos, Integer n = 1) noexcept {
    unchecked_advance (pos, static_cast<IteratorDiffT> (n));
    return pos;
  }

  template <typename Iterator,
            typename IteratorDiffT = typename std::iterator_traits<Iterator>::difference_type,
            typename Integer = IteratorDiffT>
  PLUMED_GCH_NODISCARD
  static PLUMED_GCH_CPP17_CONSTEXPR
  Iterator
  unchecked_prev (Iterator pos, Integer n = 1) noexcept {
    unchecked_advance (pos, -static_cast<IteratorDiffT> (n));
    return pos;
  }

  template <typename Iterator,
            typename IteratorDiffT = typename std::iterator_traits<Iterator>::difference_type,
            typename Integer = IteratorDiffT>
  static PLUMED_GCH_CPP17_CONSTEXPR
  void
  unchecked_advance (Iterator& pos, Integer n) noexcept {
    std::advance (pos, static_cast<IteratorDiffT> (n));
  }

  PLUMED_GCH_NODISCARD PLUMED_GCH_CPP20_CONSTEXPR
  size_ty
  get_max_size (void) const noexcept {
    // This is protected from max/min macros.
    return (std::min) (static_cast<size_ty> (alloc_traits::max_size (allocator_ref ())),
                       static_cast<size_ty> (numeric_max<difference_type> ()));
  }

  PLUMED_GCH_NODISCARD PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  allocate (size_ty n) {
    return alloc_traits::allocate (allocator_ref (), static_cast<size_type> (n));
  }

  PLUMED_GCH_NODISCARD PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  allocate_with_hint (size_ty n, cptr hint) {
    return alloc_traits::allocate (allocator_ref (), static_cast<size_type> (n), hint);
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  deallocate (ptr p, size_ty n) {
    alloc_traits::deallocate (allocator_ref (), to_address (p),
                              static_cast<size_type> (n));
  }

  template <typename U,
            typename std::enable_if<
              is_uninitialized_memcpyable<value_ty, U>::value>::type * = nullptr>
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  construct (ptr p, U&& val) noexcept {
#ifdef PLUMED_GCH_LIB_IS_CONSTANT_EVALUATED
    if (std::is_constant_evaluated ()) {
      alloc_traits::construct (allocator_ref (), to_address (p), std::forward<U> (val));
      return;
    }
#endif
    std::memcpy (to_address (p), &val, sizeof (value_ty));
  }

  // basically alloc_traits::construct
  // all this is so we can replicate C++20 behavior in the other overload
  template <typename A = alloc_ty, typename V = value_ty, typename ...Args,
            typename std::enable_if<(  sizeof...(Args) != 1
                                       ||! is_uninitialized_memcpyable<V, Args...>::value)
                                    &&  has_alloc_construct<A, V, Args...>::value>::type * = nullptr>
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  construct (ptr p, Args&&... args)
  noexcept (noexcept (alloc_traits::construct (std::declval<alloc_ty&> (),
                      std::declval<value_ty *> (),
                      std::forward<Args> (args)...))) {
    alloc_traits::construct (allocator_ref (), to_address (p), std::forward<Args> (args)...);
  }

  template <typename A = alloc_ty, typename V = value_ty, typename ...Args,
            void_t<typename std::enable_if<(  sizeof...(Args) != 1
                   ||! is_uninitialized_memcpyable<V, Args...>::value)
                   &&! has_alloc_construct<A, V, Args...>::value>::type,
                   decltype (::new (std::declval<void *> ()) V (std::declval<Args> ()...))
                   > * = nullptr>
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  construct (ptr p, Args&&... args)
  noexcept (noexcept (::new (std::declval<void *> ()) value_ty (std::declval<Args> ()...))) {
    construct_at (to_address (p), std::forward<Args> (args)...);
  }

  template <typename A = alloc_ty, typename V = value_ty,
            typename std::enable_if<is_trivially_destructible<V>::value
                                    &&! must_use_alloc_destroy<A, V>::value>::type * = nullptr>
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  destroy (ptr) const noexcept
  { }

  template <typename A = alloc_ty, typename V = value_ty,
            typename std::enable_if<(! is_trivially_destructible<V>::value
                                     ||  must_use_alloc_destroy<A, V>::value)
                                    &&  has_alloc_destroy<A, V>::value>::type * = nullptr>
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  destroy (ptr p) noexcept {
    alloc_traits::destroy (allocator_ref (), to_address (p));
  }

  // defined so we match C++20 behavior in all cases.
  template <typename A = alloc_ty, typename V = value_ty,
            typename std::enable_if<(! is_trivially_destructible<V>::value
                                     ||  must_use_alloc_destroy<A, V>::value)
                                    &&! has_alloc_destroy<A, V>::value>::type * = nullptr>
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  destroy (ptr p) noexcept {
    destroy_at (to_address (p));
  }

  template <typename A = alloc_ty, typename V = value_ty,
            typename std::enable_if<is_trivially_destructible<V>::value
                                    &&! must_use_alloc_destroy<A, V>::value>::type * = nullptr>
  PLUMED_GCH_CPP14_CONSTEXPR
  void
  destroy_range (ptr, ptr) const noexcept
  { }

  template <typename A = alloc_ty, typename V = value_ty,
            typename std::enable_if<! is_trivially_destructible<V>::value
                                    ||  must_use_alloc_destroy<A, V>::value>::type * = nullptr>
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  destroy_range (ptr first, ptr last) noexcept {
    for (; ! (first == last); ++first) {
      destroy (first);
    }
  }

  // allowed if trivially copyable and we use the standard allocator
  // and InputIt is a contiguous iterator
  template <typename ForwardIt,
            typename std::enable_if<
              is_uninitialized_memcpyable_iterator<ForwardIt>::value, bool>::type = true>
  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  uninitialized_copy (ForwardIt first, ForwardIt last, ptr dest) noexcept {
    static_assert (std::is_constructible<value_ty, decltype (*first)>::value,
                   "`value_type` must be copy constructible.");

#ifdef PLUMED_GCH_LIB_IS_CONSTANT_EVALUATED
    if (std::is_constant_evaluated ()) {
      return default_uninitialized_copy (first, last, dest);
    }
#endif

    const size_ty num_copy = external_range_length (first, last);
    if (num_copy != 0) {
      std::memcpy (to_address (dest), to_address (first), num_copy * sizeof (value_ty));
    }
    return unchecked_next (dest, num_copy);
  }

  template <typename ForwardIt,
            typename std::enable_if<
              is_uninitialized_memcpyable_iterator<ForwardIt>::value, bool>::type = true>
  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  uninitialized_copy (std::move_iterator<ForwardIt> first,
                      std::move_iterator<ForwardIt> last,
                      ptr dest) noexcept {
    return uninitialized_copy (first.base (), last.base (), dest);
  }

  template <typename InputIt,
            typename std::enable_if<
              ! is_uninitialized_memcpyable_iterator<InputIt>::value, bool>::type = false>
  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  uninitialized_copy (InputIt first, InputIt last, ptr d_first) {
    return default_uninitialized_copy (first, last, d_first);
  }

  template <typename InputIt>
  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  default_uninitialized_copy (InputIt first, InputIt last, ptr d_first) {
    ptr d_last = d_first;
    PLUMED_GCH_TRY {
      // Note: Not != because `using namespace std::rel_ops` can break constexpr.
      for (; ! (first == last); ++first, static_cast<void> (++d_last)) {
        construct (d_last, *first);
      }
      return d_last;
    }
    PLUMED_GCH_CATCH (...) {
      destroy_range (d_first, d_last);
      PLUMED_GCH_THROW;
    }
  }

  template <typename A = alloc_ty, typename V = value_ty,
            typename std::enable_if<is_trivially_constructible<V>::value
                                    &&! must_use_alloc_construct<A, V>::value>::type * = nullptr>
  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  uninitialized_value_construct (ptr first, ptr last) {
#ifdef PLUMED_GCH_LIB_IS_CONSTANT_EVALUATED
    if (std::is_constant_evaluated ()) {
      return default_uninitialized_value_construct (first, last);
    }
#endif
    std::fill (first, last, value_ty ());
    return last;
  }

  template <typename A = alloc_ty, typename V = value_ty,
            typename std::enable_if<! is_trivially_constructible<V>::value
                                    ||  must_use_alloc_construct<A, V>::value>::type * = nullptr>
  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  uninitialized_value_construct (ptr first, ptr last) {
    return default_uninitialized_value_construct (first, last);
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  default_uninitialized_value_construct (ptr first, ptr last) {
    ptr curr = first;
    PLUMED_GCH_TRY {
      for (; ! (curr == last); ++curr) {
        construct (curr);
      }
      return curr;
    }
    PLUMED_GCH_CATCH (...) {
      destroy_range (first, curr);
      PLUMED_GCH_THROW;
    }
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  uninitialized_fill (ptr first, ptr last) {
    return uninitialized_value_construct (first, last);
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  uninitialized_fill (ptr first, ptr last, const value_ty& val) {
    ptr curr = first;
    PLUMED_GCH_TRY {
      for (; ! (curr == last); ++curr) {
        construct (curr, val);
      }
      return curr;
    }
    PLUMED_GCH_CATCH (...) {
      destroy_range (first, curr);
      PLUMED_GCH_THROW;
    }
  }

private:
  // If value_ty is an array, replicate C++20 behavior (I don't think that value_ty can
  // actually be an array because of the Erasable requirement, but there shouldn't
  // be any runtime cost for being defensive here).
  template <typename V = value_ty,
            typename std::enable_if<std::is_array<V>::value, bool>::type = true>
  static PLUMED_GCH_CPP20_CONSTEXPR
  void
  destroy_at (value_ty *p) noexcept {
    for (auto& e : *p) {
      destroy_at (std::addressof (e));
    }
  }

  template <typename V = value_ty,
            typename std::enable_if<! std::is_array<V>::value, bool>::type = false>
  static PLUMED_GCH_CPP20_CONSTEXPR
  void
  destroy_at (value_ty *p) noexcept {
    p->~value_ty ();
  }

  template <typename V = value_ty, typename ...Args>
  static PLUMED_GCH_CPP20_CONSTEXPR
  auto
  construct_at (value_ty *p, Args&&... args)
  noexcept (noexcept (::new (std::declval<void *> ()) V (std::declval<Args> ()...)))
  -> decltype (::new (std::declval<void *> ()) V (std::declval<Args> ()...)) {
#if defined (PLUMED_GCH_LIB_IS_CONSTANT_EVALUATED) && defined (PLUMED_GCH_LIB_CONSTEXPR_MEMORY)
    if (std::is_constant_evaluated ()) {
      return std::construct_at (p, std::forward<Args> (args)...);
    }
#endif
    void *vp = const_cast<void *> (static_cast<const volatile void *> (p));
    return ::new (vp) value_ty (std::forward<Args>(args)...);
  }
};

template <typename Pointer, typename SizeT>
class small_vector_data_base {
public:
  using ptr     = Pointer;
  using size_ty = SizeT;

  small_vector_data_base            (void)                              = default;
  small_vector_data_base            (const small_vector_data_base&)     = default;
  small_vector_data_base            (small_vector_data_base&&) noexcept = default;
  small_vector_data_base& operator= (const small_vector_data_base&)     = default;
  small_vector_data_base& operator= (small_vector_data_base&&) noexcept = default;
  ~small_vector_data_base           (void)                              = default;

  constexpr ptr     data_ptr (void) const noexcept {
    return m_data_ptr;
  }
  constexpr size_ty capacity (void) const noexcept {
    return m_capacity;
  }
  constexpr size_ty size     (void) const noexcept {
    return m_size;
  }

  PLUMED_GCH_CPP20_CONSTEXPR void set_data_ptr (ptr     data_ptr) noexcept {
    m_data_ptr = data_ptr;
  }
  PLUMED_GCH_CPP20_CONSTEXPR void set_capacity (size_ty capacity) noexcept {
    m_capacity = capacity;
  }
  PLUMED_GCH_CPP20_CONSTEXPR void set_size     (size_ty size)     noexcept {
    m_size     = size;
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  set (ptr data_ptr, size_ty capacity, size_ty size) {
    m_data_ptr = data_ptr;
    m_capacity = capacity;
    m_size     = size;
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  swap_data_ptr (small_vector_data_base& other) noexcept {
    using std::swap;
    swap (m_data_ptr, other.m_data_ptr);
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  swap_capacity (small_vector_data_base& other) noexcept {
    using std::swap;
    swap (m_capacity, other.m_capacity);
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  swap_size (small_vector_data_base& other) noexcept {
    using std::swap;
    swap (m_size, other.m_size);
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  swap (small_vector_data_base& other) noexcept {
    using std::swap;
    swap (m_data_ptr, other.m_data_ptr);
    swap (m_capacity, other.m_capacity);
    swap (m_size,     other.m_size);
  }

private:
  ptr     m_data_ptr;
  size_ty m_capacity;
  size_ty m_size;
};

template <typename Pointer, typename SizeT, typename T, unsigned InlineCapacity>
class small_vector_data
: public small_vector_data_base<Pointer, SizeT> {
public:
  using value_ty = T;

  small_vector_data            (void)                         = default;
  small_vector_data            (const small_vector_data&)     = delete;
  small_vector_data            (small_vector_data&&) noexcept = delete;
  small_vector_data& operator= (const small_vector_data&)     = delete;
  small_vector_data& operator= (small_vector_data&&) noexcept = delete;
  ~small_vector_data           (void)                         = default;

  PLUMED_GCH_CPP14_CONSTEXPR
  value_ty *
  storage (void) noexcept {
    return m_storage.get_inline_ptr ();
  }

  constexpr
  const value_ty *
  storage (void) const noexcept {
    return m_storage.get_inline_ptr ();
  }

private:
  inline_storage<value_ty, InlineCapacity> m_storage;
};

template <typename Pointer, typename SizeT, typename T>
class PLUMED_GCH_EMPTY_BASE small_vector_data<Pointer, SizeT, T, 0>
: public  small_vector_data_base<Pointer, SizeT>,
private inline_storage<T, 0> {
  using base = inline_storage<T, 0>;

public:
  using value_ty = T;

  small_vector_data            (void)                         = default;
  small_vector_data            (const small_vector_data&)     = delete;
  small_vector_data            (small_vector_data&&) noexcept = delete;
  small_vector_data& operator= (const small_vector_data&)     = delete;
  small_vector_data& operator= (small_vector_data&&) noexcept = delete;
  ~small_vector_data           (void)                         = default;

  PLUMED_GCH_CPP14_CONSTEXPR
  value_ty *
  storage (void) noexcept {
    return base::get_inline_ptr ();
  }

  constexpr
  const value_ty *
  storage (void) const noexcept {
    return base::get_inline_ptr ();
  }
};

template <typename Allocator, unsigned InlineCapacity>
class small_vector_base
: public allocator_interface<Allocator> {
public:
  using size_type       = typename allocator_interface<Allocator>::size_type;
  using difference_type = typename allocator_interface<Allocator>::difference_type;

  template <typename SameAllocator, unsigned DifferentInlineCapacity>
  friend class small_vector_base;

protected:
  using alloc_interface = allocator_interface<Allocator>;
  using alloc_traits    = typename alloc_interface::alloc_traits;
  using alloc_ty        = Allocator;

  using value_ty        = typename alloc_interface::value_ty;
  using ptr             = typename alloc_interface::ptr;
  using cptr            = typename alloc_interface::cptr;
  using size_ty         = typename alloc_interface::size_ty;
  using diff_ty         = typename alloc_interface::diff_ty;

  static_assert (alloc_interface::template is_complete<value_ty>::value || InlineCapacity == 0,
                 "`value_type` must be complete for instantiation of a non-zero number "
                 "of inline elements.");

  template <typename T>
  using is_complete = typename alloc_interface::template is_complete<T>;

  using alloc_interface::allocator_ref;
  using alloc_interface::construct;
  using alloc_interface::deallocate;
  using alloc_interface::destroy;
  using alloc_interface::destroy_range;
  using alloc_interface::external_range_length;
  using alloc_interface::get_max_size;
  using alloc_interface::internal_range_length;
  using alloc_interface::to_address;
  using alloc_interface::unchecked_advance;
  using alloc_interface::unchecked_next;
  using alloc_interface::unchecked_prev;
  using alloc_interface::uninitialized_copy;
  using alloc_interface::uninitialized_fill;
  using alloc_interface::uninitialized_value_construct;

  template <typename Integer>
  PLUMED_GCH_NODISCARD
  static PLUMED_GCH_CONSTEVAL
  std::size_t
  numeric_max (void) noexcept {
    return alloc_interface::template numeric_max<Integer> ();
  }

  PLUMED_GCH_NODISCARD
  static PLUMED_GCH_CONSTEVAL
  size_ty
  get_inline_capacity (void) noexcept {
    return static_cast<size_ty> (InlineCapacity);
  }

  template <typename ...>
  using void_t = void;

  template <bool B>
  using bool_constant = std::integral_constant<bool, B>;

  template <typename Void, typename AI, typename V, typename ...Args>
  struct is_emplace_constructible_impl
: std::false_type {
    using nothrow = std::false_type;
  };

  template <typename AI, typename V, typename ...Args>
  struct is_emplace_constructible_impl<
    void_t<typename std::enable_if<is_complete<V>::value>::type,
           decltype (std::declval<AI&> ().construct (std::declval<V *> (),
                     std::declval<Args> ()...))>,
    AI, V, Args...>
: std::true_type {
    using nothrow =
    bool_constant<noexcept (std::declval<AI&> ().construct (std::declval<V *> (),
                            std::declval<Args> ()...))>;
                                                           };

  template <typename ...Args>
  struct is_emplace_constructible
: is_emplace_constructible_impl<void, alloc_interface, value_ty, Args...> {
  };

  template <typename ...Args>
  struct is_nothrow_emplace_constructible
: is_emplace_constructible_impl<void, alloc_interface, value_ty, Args...>::nothrow {
  };

  template <typename V = value_ty>
  struct is_explicitly_move_insertable
: is_emplace_constructible<V&&> {
  };

  template <typename V = value_ty>
  struct is_explicitly_nothrow_move_insertable
: is_nothrow_emplace_constructible<V&&> {
  };

  template <typename V = value_ty>
  struct is_explicitly_copy_insertable
: std::integral_constant<bool, is_emplace_constructible<V&>::value
                           &&  is_emplace_constructible<const V&>::value> {
                           };

  template <typename V = value_ty>
  struct is_explicitly_nothrow_copy_insertable
: std::integral_constant<bool, is_nothrow_emplace_constructible<V&>::value
                           &&  is_nothrow_emplace_constructible<const V&>::value> {
                           };

  template <typename AI, typename Enable = void>
  struct is_eraseable
: std::false_type {
  };

  template <typename AI>
  struct is_eraseable<AI,
                      void_t<decltype (std::declval<AI&> ().destroy (std::declval<value_ty *> ()))>>
: std::true_type {
  };

  template <typename V>
  struct relocate_with_move
#ifdef PLUMED_GCH_NO_STRONG_EXCEPTION_GUARANTEES
: std::true_type
#else
: bool_constant<std::is_nothrow_move_constructible<V>::value
                  ||! is_explicitly_copy_insertable<V>::value>
#endif
  { };

  template <typename A>
  struct allocations_are_movable
: bool_constant<std::is_same<std::allocator<value_ty>, A>::value
                  ||  std::allocator_traits<A>::propagate_on_container_move_assignment::value
#ifdef PLUMED_GCH_LIB_IS_ALWAYS_EQUAL
                  ||  std::allocator_traits<A>::is_always_equal::value
#endif
                  > {
                  };

  template <typename A>
  struct allocations_are_swappable
: bool_constant<std::is_same<std::allocator<value_ty>, A>::value
                  ||  std::allocator_traits<A>::propagate_on_container_swap::value
#ifdef PLUMED_GCH_LIB_IS_ALWAYS_EQUAL
                  ||  std::allocator_traits<A>::is_always_equal::value
#endif
                  > {
                  };

  template <typename ...Args>
  using is_memcpyable = typename alloc_interface::template is_memcpyable<Args...>;

  template <typename ...Args>
  using is_memcpyable_iterator =
  typename alloc_interface::template is_memcpyable_iterator<Args...>;

  PLUMED_GCH_NORETURN
  static PLUMED_GCH_CPP20_CONSTEXPR
  void
  throw_overflow_error (void) {
#ifdef PLUMED_GCH_EXCEPTIONS
    throw std::overflow_error ("The requested conversion would overflow.");
#else
    std::fprintf (stderr, "[gch::small_vector] The requested conversion would overflow.\n");
    std::abort ();
#endif
  }

  PLUMED_GCH_NORETURN
  static PLUMED_GCH_CPP20_CONSTEXPR
  void
  throw_index_error (void) {
#ifdef PLUMED_GCH_EXCEPTIONS
    throw std::out_of_range ("The requested index was out of range.");
#else
    std::fprintf (stderr, "[gch::small_vector] The requested index was out of range.\n");
    std::abort ();
#endif
  }

  PLUMED_GCH_NORETURN
  static PLUMED_GCH_CPP20_CONSTEXPR
  void
  throw_increment_error (void) {
#ifdef PLUMED_GCH_EXCEPTIONS
    throw std::domain_error ("The requested increment was outside of the allowed range.");
#else
    std::fprintf (
      stderr,
      "[gch::small_vector] The requested increment was outside of the allowed range.\n");
    std::abort ();
#endif
  }

  PLUMED_GCH_NORETURN
  static PLUMED_GCH_CPP20_CONSTEXPR
  void
  throw_allocation_size_error (void) {
#ifdef PLUMED_GCH_EXCEPTIONS
    throw std::length_error ("The required allocation exceeds the maximum size.");
#else
    std::fprintf (
      stderr,
      "[gch::small_vector] The required allocation exceeds the maximum size.\n");
    std::abort ();
#endif
  }

  PLUMED_GCH_NODISCARD PLUMED_GCH_CPP14_CONSTEXPR
  ptr
  ptr_cast (const small_vector_iterator<cptr, diff_ty>& it) noexcept {
    return unchecked_next (begin_ptr (), it.base () - begin_ptr ());
  }

private:
  class stack_temporary {
  public:
    stack_temporary            (void)                       = delete;
    stack_temporary            (const stack_temporary&)     = delete;
    stack_temporary            (stack_temporary&&) noexcept = delete;
    stack_temporary& operator= (const stack_temporary&)     = delete;
    stack_temporary& operator= (stack_temporary&&) noexcept = delete;
//      ~stack_temporary           (void)                       = impl;

    template <typename ...Args>
    PLUMED_GCH_CPP20_CONSTEXPR explicit
    stack_temporary (alloc_interface& alloc_iface, Args&&... args)
      : m_interface (alloc_iface) {
      m_interface.construct (get_pointer (), std::forward<Args> (args)...);
    }

    PLUMED_GCH_CPP20_CONSTEXPR
    ~stack_temporary (void) {
      m_interface.destroy (get_pointer ());
    }

    PLUMED_GCH_NODISCARD PLUMED_GCH_CPP20_CONSTEXPR
    const value_ty&
    get (void) const noexcept {
      return *get_pointer ();
    }

    PLUMED_GCH_NODISCARD PLUMED_GCH_CPP20_CONSTEXPR
    value_ty&&
    release (void) noexcept {
      return std::move (*get_pointer ());
    }

  private:
    PLUMED_GCH_NODISCARD PLUMED_GCH_CPP20_CONSTEXPR
    cptr
    get_pointer (void) const noexcept {
      return static_cast<cptr> (static_cast<const void *> (std::addressof (m_data)));
    }

    PLUMED_GCH_NODISCARD PLUMED_GCH_CPP20_CONSTEXPR
    ptr
    get_pointer (void) noexcept {
      return static_cast<ptr> (static_cast<void *> (std::addressof (m_data)));
    }

    alloc_interface& m_interface;
    typename std::aligned_storage<sizeof (value_ty), alignof (value_ty)>::type m_data;
  };

#ifdef PLUMED_GCH_LIB_IS_CONSTANT_EVALUATED

  class heap_temporary {
  public:
    heap_temporary            (void)                      = delete;
    heap_temporary            (const heap_temporary&)     = delete;
    heap_temporary            (heap_temporary&&) noexcept = delete;
    heap_temporary& operator= (const heap_temporary&)     = delete;
    heap_temporary& operator= (heap_temporary&&) noexcept = delete;
//      ~heap_temporary           (void)                      = impl;

    template <typename ...Args>
    PLUMED_GCH_CPP20_CONSTEXPR explicit
    heap_temporary (alloc_interface& alloc_iface, Args&&... args)
      : m_interface (alloc_iface),
        m_data_ptr  (alloc_iface.allocate (sizeof (value_ty))) {
      PLUMED_GCH_TRY {
        m_interface.construct (m_data_ptr, std::forward<Args> (args)...);
      }
      PLUMED_GCH_CATCH (...) {
        m_interface.deallocate (m_data_ptr, sizeof (value_ty));
        PLUMED_GCH_THROW;
      }
    }

    PLUMED_GCH_CPP20_CONSTEXPR
    ~heap_temporary (void) {
      m_interface.destroy (m_data_ptr);
      m_interface.deallocate (m_data_ptr, sizeof (value_ty));
    }

    PLUMED_GCH_NODISCARD PLUMED_GCH_CPP20_CONSTEXPR
    const value_ty&
    get (void) const noexcept {
      return *m_data_ptr;
    }

    PLUMED_GCH_NODISCARD PLUMED_GCH_CPP20_CONSTEXPR
    value_ty&&
    release (void) noexcept {
      return std::move (*m_data_ptr);
    }

  private:
    alloc_interface& m_interface;
    ptr              m_data_ptr;
  };

#endif

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  wipe (void) {
    destroy_range (begin_ptr (), end_ptr ());
    if (has_allocation ()) {
      deallocate (data_ptr (), get_capacity ());
    }
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  set_data_ptr (ptr data_ptr) noexcept {
    m_data.set_data_ptr (data_ptr);
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  set_capacity (size_ty capacity) noexcept {
    m_data.set_capacity (static_cast<size_type> (capacity));
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  set_size (size_ty size) noexcept {
    m_data.set_size (static_cast<size_type> (size));
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  set_data (ptr data_ptr, size_ty capacity, size_ty size) noexcept {
    m_data.set (data_ptr, static_cast<size_type> (capacity), static_cast<size_type> (size));
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  swap_data_ptr (small_vector_base& other) noexcept {
    m_data.swap_data_ptr (other.m_data);
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  swap_capacity (small_vector_base& other) noexcept {
    m_data.swap_capacity (other.m_data);
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  swap_size (small_vector_base& other) noexcept {
    m_data.swap_size (other.m_data);
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  swap_allocation (small_vector_base& other) noexcept {
    m_data.swap (other.m_data);
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  reset_data (ptr data_ptr, size_ty capacity, size_ty size) {
    wipe ();
    m_data.set (data_ptr, static_cast<size_type> (capacity), static_cast<size_type> (size));
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  increase_size (size_ty n) noexcept {
    m_data.set_size (get_size () + n);
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  decrease_size (size_ty n) noexcept {
    m_data.set_size (get_size () - n);
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  unchecked_allocate (size_ty n) {
    assert (InlineCapacity < n && "Allocated capacity should be greater than InlineCapacity.");
    return alloc_interface::allocate (n);
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  unchecked_allocate (size_ty n, cptr hint) {
    assert (InlineCapacity < n && "Allocated capacity should be greater than InlineCapacity.");
    return alloc_interface::allocate_with_hint (n, hint);
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  checked_allocate (size_ty n) {
    if (get_max_size () < n) {
      throw_allocation_size_error ();
    }
    return unchecked_allocate (n);
  }

protected:
  PLUMED_GCH_NODISCARD PLUMED_GCH_CPP14_CONSTEXPR
  size_ty
  unchecked_calculate_new_capacity (const size_ty minimum_required_capacity) const noexcept {
    const size_ty current_capacity = get_capacity ();

    assert (current_capacity < minimum_required_capacity);

    if (get_max_size () - current_capacity <= current_capacity) {
      return get_max_size ();
    }

    // Note: This growth factor might be theoretically superior, but in testing it falls flat:
    // size_ty new_capacity = current_capacity + (current_capacity / 2);

    const size_ty new_capacity = 2 * current_capacity;
    if (new_capacity < minimum_required_capacity) {
      return minimum_required_capacity;
    }
    return new_capacity;
  }

  PLUMED_GCH_NODISCARD PLUMED_GCH_CPP14_CONSTEXPR
  size_ty
  checked_calculate_new_capacity (const size_ty minimum_required_capacity) const {
    if (get_max_size () < minimum_required_capacity) {
      throw_allocation_size_error ();
    }
    return unchecked_calculate_new_capacity (minimum_required_capacity);
  }

  template <unsigned I>
  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector_base&
  copy_assign_default (const small_vector_base<Allocator, I>& other) {
    if (get_capacity () < other.get_size ()) {
      // Reallocate.
      size_ty new_capacity = unchecked_calculate_new_capacity (other.get_size ());
      ptr     new_data_ptr = unchecked_allocate (new_capacity, other.allocation_end_ptr ());

      PLUMED_GCH_TRY {
        uninitialized_copy (other.begin_ptr (), other.end_ptr (), new_data_ptr);
      }
      PLUMED_GCH_CATCH (...) {
        deallocate (new_data_ptr, new_capacity);
        PLUMED_GCH_THROW;
      }

      reset_data (new_data_ptr, new_capacity, other.get_size ());
    } else {
      if (get_size () < other.get_size ()) {
        // No reallocation, partially in uninitialized space.
        std::copy_n (other.begin_ptr (), get_size (), begin_ptr ());
        uninitialized_copy (
          unchecked_next (other.begin_ptr (), get_size ()),
          other.end_ptr (),
          end_ptr ());
      } else {
        destroy_range (copy_range (other.begin_ptr (), other.end_ptr (), begin_ptr ()),
                       end_ptr ());
      }

      // data_ptr and capacity do not change in this case.
      set_size (other.get_size ());
    }

    alloc_interface::operator= (other);
    return *this;
  }

  template <unsigned I, typename AT = alloc_traits,
            typename std::enable_if<AT::propagate_on_container_copy_assignment::value
#ifdef PLUMED_GCH_LIB_IS_ALWAYS_EQUAL
                                    &&! AT::is_always_equal::value
#endif
                                    >::type * = nullptr>
  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector_base&
  copy_assign (const small_vector_base<Allocator, I>& other) {
    if (other.allocator_ref () == allocator_ref ()) {
      return copy_assign_default (other);
    }

    if (InlineCapacity < other.get_size ()) {
      alloc_interface new_alloc (other);

      const size_ty new_capacity = other.get_size ();
      const ptr new_data_ptr = new_alloc.allocate_with_hint (
                                 new_capacity,
                                 other.allocation_end_ptr ());

      PLUMED_GCH_TRY {
        uninitialized_copy (other.begin_ptr (), other.end_ptr (), new_data_ptr);
      }
      PLUMED_GCH_CATCH (...) {
        new_alloc.deallocate (new_data_ptr, new_capacity);
        PLUMED_GCH_THROW;
      }

      reset_data (new_data_ptr, new_capacity, other.get_size ());
      alloc_interface::operator= (new_alloc);
    } else {
      if (has_allocation ()) {
#ifdef PLUMED_GCH_LIB_IS_CONSTANT_EVALUATED
        ptr new_data_ptr;
        if (std::is_constant_evaluated ()) {
          alloc_interface new_alloc (other);
          new_data_ptr = new_alloc.allocate (InlineCapacity);
        } else {
          new_data_ptr = storage_ptr ();
        }
#else
        const ptr new_data_ptr = storage_ptr ();
#endif

        uninitialized_copy (other.begin_ptr (), other.end_ptr (), new_data_ptr);
        destroy_range (begin_ptr (), end_ptr ());
        deallocate (data_ptr (), get_capacity ());
        set_data_ptr (new_data_ptr);
        set_capacity (InlineCapacity);
      } else if (get_size () < other.get_size ()) {
        std::copy_n (other.begin_ptr (), get_size (), begin_ptr ());
        uninitialized_copy (
          unchecked_next (other.begin_ptr (), get_size ()),
          other.end_ptr (),
          end_ptr ());
      } else {
        destroy_range (copy_range (other.begin_ptr (), other.end_ptr (), begin_ptr ()),
                       end_ptr ());
      }
      set_size (other.get_size ());
      alloc_interface::operator= (other);
    }

    return *this;
  }

  template <unsigned I, typename AT = alloc_traits,
            typename std::enable_if<! AT::propagate_on_container_copy_assignment::value
#ifdef PLUMED_GCH_LIB_IS_ALWAYS_EQUAL
                                    ||  AT::is_always_equal::value
#endif
                                    >::type * = nullptr>
  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector_base&
  copy_assign (const small_vector_base<Allocator, I>& other) {
    return copy_assign_default (other);
  }

  template <unsigned I>
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  move_allocation_pointer (small_vector_base<alloc_ty, I>&& other) noexcept {
    reset_data (other.data_ptr (), other.get_capacity (), other.get_size ());
    other.set_default ();
  }

  template <unsigned N = InlineCapacity, typename std::enable_if<N == 0>::type * = nullptr>
  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector_base&
  move_assign_default (small_vector_base&& other) noexcept {
    move_allocation_pointer (std::move (other));
    alloc_interface::operator= (std::move (other));
    return *this;
  }

  template <unsigned LessEqualI,
            typename std::enable_if<(LessEqualI <= InlineCapacity)>::type * = nullptr>
            PLUMED_GCH_CPP20_CONSTEXPR
            small_vector_base&
            move_assign_default (small_vector_base<Allocator, LessEqualI>&& other)
            noexcept (std::is_nothrow_move_assignable<value_ty>::value
  &&  std::is_nothrow_move_constructible<value_ty>::value) {
    // We only move the allocation pointer over if it has strictly greater capacity than
    // the inline capacity of `*this` because allocations can never have a smaller capacity
    // than the inline capacity.
    if (InlineCapacity < other.get_capacity ()) {
      move_allocation_pointer (std::move (other));
    } else {
      // We are guaranteed to have sufficient capacity to store the elements.
      if (InlineCapacity < get_capacity ()) {
#ifdef PLUMED_GCH_LIB_IS_CONSTANT_EVALUATED
        ptr new_data_ptr;
        if (std::is_constant_evaluated ()) {
          new_data_ptr = other.allocate (InlineCapacity);
        } else {
          new_data_ptr = storage_ptr ();
        }
#else
        const ptr new_data_ptr = storage_ptr ();
#endif

        uninitialized_move (other.begin_ptr (), other.end_ptr (), new_data_ptr);
        destroy_range (begin_ptr (), end_ptr ());
        deallocate (data_ptr (), get_capacity ());
        set_data_ptr (new_data_ptr);
        set_capacity (InlineCapacity);
      } else if (get_size () < other.get_size ()) {
        // There are more elements in `other`.
        // Overwrite the existing range and uninitialized move the rest.
        ptr other_pivot = unchecked_next (other.begin_ptr (), get_size ());
        std::move (other.begin_ptr (), other_pivot, begin_ptr ());
        uninitialized_move (other_pivot, other.end_ptr (), end_ptr ());
      } else {
        // There are the same number or fewer elements in `other`.
        // Overwrite part of the existing range and destroy the rest.
        ptr new_end = std::move (other.begin_ptr (), other.end_ptr (), begin_ptr ());
        destroy_range (new_end, end_ptr ());
      }

      set_size (other.get_size ());

      // Note: We do not need to deallocate any allocations in `other` because the value of
      //       an object meeting the Allocator named requirements does not change value after
      //       a move.
    }

    alloc_interface::operator= (std::move (other));
    return *this;
  }

  template <unsigned GreaterI,
            typename std::enable_if<(InlineCapacity < GreaterI)>::type * = nullptr>
            PLUMED_GCH_CPP20_CONSTEXPR
            small_vector_base&
  move_assign_default (small_vector_base<Allocator, GreaterI>&& other) {
    if (other.has_allocation ()) {
      move_allocation_pointer (std::move (other));
    } else if (get_capacity () < other.get_size ()
               ||  (has_allocation () && ! (other.allocator_ref () == allocator_ref ()))) {
      // Reallocate.

      // The compiler should be able to optimize this.
      size_ty new_capacity =
        get_capacity () < other.get_size ()
        ? unchecked_calculate_new_capacity (other.get_size ())
        : get_capacity ();

      ptr new_data_ptr = other.allocate_with_hint (new_capacity, other.allocation_end_ptr ());

      PLUMED_GCH_TRY {
        uninitialized_move (other.begin_ptr (), other.end_ptr (), new_data_ptr);
      }
      PLUMED_GCH_CATCH (...) {
        other.deallocate (new_data_ptr, new_capacity);
        PLUMED_GCH_THROW;
      }

      reset_data (new_data_ptr, new_capacity, other.get_size ());
    } else {
      if (get_size () < other.get_size ()) {
        // There are more elements in `other`.
        // Overwrite the existing range and uninitialized move the rest.
        ptr other_pivot = unchecked_next (other.begin_ptr (), get_size ());
        std::move (other.begin_ptr (), other_pivot, begin_ptr ());
        uninitialized_move (other_pivot, other.end_ptr (), end_ptr ());
      } else {
        // fewer elements in other
        // overwrite part of the existing range and destroy the rest
        ptr new_end = std::move (other.begin_ptr (), other.end_ptr (), begin_ptr ());
        destroy_range (new_end, end_ptr ());
      }

      // `data_ptr` and `capacity` do not change in this case.
      set_size (other.get_size ());
    }

    alloc_interface::operator= (std::move (other));
    return *this;
  }

  template <unsigned I>
  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector_base&
  move_assign_unequal_no_propagate (small_vector_base<Allocator, I>&& other) {
    if (get_capacity () < other.get_size ()) {
      // Reallocate.
      size_ty new_capacity = unchecked_calculate_new_capacity (other.get_size ());
      ptr     new_data_ptr = unchecked_allocate (new_capacity, other.allocation_end_ptr ());

      PLUMED_GCH_TRY {
        uninitialized_move (other.begin_ptr (), other.end_ptr (), new_data_ptr);
      }
      PLUMED_GCH_CATCH (...) {
        deallocate (new_data_ptr, new_capacity);
        PLUMED_GCH_THROW;
      }

      reset_data (new_data_ptr, new_capacity, other.get_size ());
    } else {
      if (get_size () < other.get_size ()) {
        // There are more elements in `other`.
        // Overwrite the existing range and uninitialized move the rest.
        ptr other_pivot = unchecked_next (other.begin_ptr (), get_size ());
        std::move (other.begin_ptr (), other_pivot, begin_ptr ());
        uninitialized_move (other_pivot, other.end_ptr (), end_ptr ());
      } else {
        // There are fewer elements in `other`.
        // Overwrite part of the existing range and destroy the rest.
        destroy_range (
          std::move (other.begin_ptr (), other.end_ptr (), begin_ptr ()),
          end_ptr ());
      }

      // data_ptr and capacity do not change in this case
      set_size (other.get_size ());
    }

    alloc_interface::operator= (std::move (other));
    return *this;
  }

  template <unsigned I, typename A = alloc_ty,
            typename std::enable_if<allocations_are_movable<A>::value>::type * = nullptr>
  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector_base&
  move_assign (small_vector_base<Allocator, I>&& other)
  noexcept (noexcept (
              std::declval<small_vector_base&> ().move_assign_default (std::move (other)))) {
    return move_assign_default (std::move (other));
  }

  template <unsigned I, typename A = alloc_ty,
            typename std::enable_if<! allocations_are_movable<A>::value>::type * = nullptr>
  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector_base&
  move_assign (small_vector_base<Allocator, I>&& other) {
    if (other.allocator_ref () == allocator_ref ()) {
      return move_assign_default (std::move (other));
    }
    return move_assign_unequal_no_propagate (std::move (other));
  }

  template <unsigned I = InlineCapacity,
            typename std::enable_if<I == 0>::type * = nullptr>
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  move_initialize (small_vector_base&& other) noexcept {
    set_data (other.data_ptr (), other.get_capacity (), other.get_size ());
    other.set_default ();
  }

  template <unsigned LessEqualI,
            typename std::enable_if<(LessEqualI <= InlineCapacity)>::type * = nullptr>
            PLUMED_GCH_CPP20_CONSTEXPR
            void
            move_initialize (small_vector_base<Allocator, LessEqualI>&& other)
  noexcept (std::is_nothrow_move_constructible<value_ty>::value) {
    if (InlineCapacity < other.get_capacity ()) {
      set_data (other.data_ptr (), other.get_capacity (), other.get_size ());
      other.set_default ();
    } else {
      set_to_inline_storage ();
      uninitialized_move (other.begin_ptr (), other.end_ptr (), data_ptr ());
      set_size (other.get_size ());
    }
  }

  template <unsigned GreaterI,
            typename std::enable_if<(InlineCapacity < GreaterI)>::type * = nullptr>
            PLUMED_GCH_CPP20_CONSTEXPR
            void
  move_initialize (small_vector_base<Allocator, GreaterI>&& other) {
    if (other.has_allocation ()) {
      set_data (other.data_ptr (), other.get_capacity (), other.get_size ());
      other.set_default ();
    } else {
      if (InlineCapacity < other.get_size ()) {
        // We may throw in this case.
        set_data_ptr (unchecked_allocate (other.get_size (), other.allocation_end_ptr ()));
        set_capacity (other.get_size ());

        PLUMED_GCH_TRY {
          uninitialized_move (other.begin_ptr (), other.end_ptr (), data_ptr ());
        }
        PLUMED_GCH_CATCH (...) {
          deallocate (data_ptr (), get_capacity ());
          PLUMED_GCH_THROW;
        }
      } else {
        set_to_inline_storage ();
        uninitialized_move (other.begin_ptr (), other.end_ptr (), data_ptr ());
      }

      set_size (other.get_size ());
    }
  }

public:
//    small_vector_base            (void)                         = impl;
  small_vector_base            (const small_vector_base&)     = delete;
  small_vector_base            (small_vector_base&&) noexcept = delete;
  small_vector_base& operator= (const small_vector_base&)     = delete;
  small_vector_base& operator= (small_vector_base&&) noexcept = delete;
//    ~small_vector_base           (void)                         = impl;

  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector_base (void) noexcept {
    set_default ();
  }

  static constexpr struct bypass_tag { } bypass { };

  template <unsigned I, typename ...MaybeAlloc>
  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector_base (bypass_tag,
                     const small_vector_base<Allocator, I>& other,
                     const MaybeAlloc&... alloc)
    : alloc_interface (other, alloc...) {
    if (InlineCapacity < other.get_size ()) {
      set_data_ptr (unchecked_allocate (other.get_size (), other.allocation_end_ptr ()));
      set_capacity (other.get_size ());

      PLUMED_GCH_TRY {
        uninitialized_copy (other.begin_ptr (), other.end_ptr (), data_ptr ());
      }
      PLUMED_GCH_CATCH (...) {
        deallocate (data_ptr (), get_capacity ());
        PLUMED_GCH_THROW;
      }
    } else {
      set_to_inline_storage ();
      uninitialized_copy (other.begin_ptr (), other.end_ptr (), data_ptr ());
    }

    set_size (other.get_size ());
  }

  template <unsigned I>
  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector_base (bypass_tag, small_vector_base<Allocator, I>&& other)
  noexcept (std::is_nothrow_move_constructible<value_ty>::value
            ||  (I == 0 && I == InlineCapacity))
    : alloc_interface (std::move (other)) {
    move_initialize (std::move (other));
  }

  template <unsigned I, typename A = alloc_ty,
            typename std::enable_if<std::is_same<std::allocator<value_ty>, A>::value
#ifdef PLUMED_GCH_LIB_IS_ALWAYS_EQUAL
                                    ||  std::allocator_traits<A>::is_always_equal::value
#endif
                                    >::type * = nullptr>
  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector_base (bypass_tag, small_vector_base<Allocator, I>&& other, const alloc_ty&)
  noexcept (noexcept (small_vector_base (bypass, std::move (other))))
    : small_vector_base (bypass, std::move (other))
  { }

  template <unsigned I, typename A = alloc_ty,
            typename std::enable_if<! (std::is_same<std::allocator<value_ty>, A>::value
#ifdef PLUMED_GCH_LIB_IS_ALWAYS_EQUAL
                                       ||  std::allocator_traits<A>::is_always_equal::value
#endif
                                                                )>::type * = nullptr>
  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector_base (bypass_tag, small_vector_base<Allocator, I>&& other, const alloc_ty& alloc)
    : alloc_interface (alloc) {
    if (other.allocator_ref () == alloc) {
      move_initialize (std::move (other));
      return;
    }

    if (InlineCapacity < other.get_size ()) {
      // We may throw in this case.
      set_data_ptr (unchecked_allocate (other.get_size (), other.allocation_end_ptr ()));
      set_capacity (other.get_size ());

      PLUMED_GCH_TRY {
        uninitialized_move (other.begin_ptr (), other.end_ptr (), data_ptr ());
      }
      PLUMED_GCH_CATCH (...) {
        deallocate (data_ptr (), get_capacity ());
        PLUMED_GCH_THROW;
      }
    } else {
      set_to_inline_storage ();
      uninitialized_move (other.begin_ptr (), other.end_ptr (), data_ptr ());
    }

    set_size (other.get_size ());
  }

  PLUMED_GCH_CPP20_CONSTEXPR explicit
  small_vector_base (const alloc_ty& alloc) noexcept
    : alloc_interface (alloc) {
    set_default ();
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector_base (size_ty count, const alloc_ty& alloc)
    : alloc_interface (alloc) {
    if (InlineCapacity < count) {
      set_data_ptr (checked_allocate (count));
      set_capacity (count);
    } else {
      set_to_inline_storage ();
    }

    PLUMED_GCH_TRY {
      uninitialized_value_construct (begin_ptr (), unchecked_next (begin_ptr (), count));
    }
    PLUMED_GCH_CATCH (...) {
      if (has_allocation ()) {
        deallocate (data_ptr (), get_capacity ());
      }
      PLUMED_GCH_THROW;
    }
    set_size (count);
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector_base (size_ty count, const value_ty& val, const alloc_ty& alloc)
    : alloc_interface (alloc) {
    if (InlineCapacity < count) {
      set_data_ptr (checked_allocate (count));
      set_capacity (count);
    } else {
      set_to_inline_storage ();
    }

    PLUMED_GCH_TRY {
      uninitialized_fill (begin_ptr (), unchecked_next (begin_ptr (), count), val);
    }
    PLUMED_GCH_CATCH (...) {
      if (has_allocation ()) {
        deallocate (data_ptr (), get_capacity ());
      }
      PLUMED_GCH_THROW;
    }
    set_size (count);
  }

  template <typename Generator>
  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector_base (size_ty count, Generator& g, const alloc_ty& alloc)
    : alloc_interface (alloc) {
    if (InlineCapacity < count) {
      set_data_ptr (checked_allocate (count));
      set_capacity (count);
    } else {
      set_to_inline_storage ();
    }

    ptr curr = begin_ptr ();
    const ptr new_end = unchecked_next (begin_ptr (), count);
    PLUMED_GCH_TRY {
      for (; ! (curr == new_end); ++curr) {
        construct (curr, g ());
      }
    }
    PLUMED_GCH_CATCH (...) {
      destroy_range (begin_ptr (), curr);
      if (has_allocation ()) {
        deallocate (data_ptr (), get_capacity ());
      }
      PLUMED_GCH_THROW;
    }
    set_size (count);
  }

#ifdef PLUMED_GCH_LIB_CONCEPTS
  template <std::input_iterator InputIt>
#else
  template <typename InputIt>
#endif
  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector_base (InputIt first, InputIt last, std::input_iterator_tag,
                     const alloc_ty& alloc)
    : small_vector_base (alloc) {
    using iterator_cat = typename std::iterator_traits<InputIt>::iterator_category;
    append_range (first, last, iterator_cat { });
  }

#ifdef PLUMED_GCH_LIB_CONCEPTS
  template <std::forward_iterator ForwardIt>
#else
  template <typename ForwardIt>
#endif
  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector_base (ForwardIt first, ForwardIt last, std::forward_iterator_tag,
                     const alloc_ty& alloc)
    : alloc_interface (alloc) {
    size_ty count = external_range_length (first, last);
    if (InlineCapacity < count) {
      set_data_ptr (unchecked_allocate (count));
      set_capacity (count);
      PLUMED_GCH_TRY {
        uninitialized_copy (first, last, begin_ptr ());
      }
      PLUMED_GCH_CATCH (...) {
        deallocate (data_ptr (), get_capacity ());
        PLUMED_GCH_THROW;
      }
    } else {
      set_to_inline_storage ();
      uninitialized_copy (first, last, begin_ptr ());
    }

    set_size (count);
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  ~small_vector_base (void) noexcept {
    assert (InlineCapacity <= get_capacity () && "Invalid capacity.");
    wipe ();
  }

protected:

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  set_to_inline_storage (void) {
    set_capacity (InlineCapacity);
#ifdef PLUMED_GCH_LIB_IS_CONSTANT_EVALUATED
    if (std::is_constant_evaluated ()) {
      return set_data_ptr (alloc_interface::allocate (InlineCapacity));
    }
#endif
    set_data_ptr (storage_ptr ());
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  assign_with_copies (size_ty count, const value_ty& val) {
    if (get_capacity () < count) {
      size_ty new_capacity = checked_calculate_new_capacity (count);
      ptr     new_begin    = unchecked_allocate (new_capacity);

      PLUMED_GCH_TRY {
        uninitialized_fill (new_begin, unchecked_next (new_begin, count), val);
      }
      PLUMED_GCH_CATCH (...) {
        deallocate (new_begin, new_capacity);
        PLUMED_GCH_THROW;
      }

      reset_data (new_begin, new_capacity, count);
    } else if (get_size () < count) {
      std::fill (begin_ptr (), end_ptr (), val);
      uninitialized_fill (end_ptr (), unchecked_next (begin_ptr (), count), val);
      set_size (count);
    } else {
      erase_range (std::fill_n (begin_ptr (), count, val), end_ptr ());
    }
  }

  template <typename InputIt,
            typename std::enable_if<std::is_assignable<
                                      value_ty&,
                                      decltype (*std::declval<InputIt> ())>::value>::type * = nullptr>
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  assign_with_range (InputIt first, InputIt last, std::input_iterator_tag) {
    using iterator_cat = typename std::iterator_traits<InputIt>::iterator_category;

    ptr curr = begin_ptr ();
    for (; ! (end_ptr () == curr || first == last); ++curr, static_cast<void> (++first)) {
      *curr = *first;
    }

    if (first == last) {
      erase_to_end (curr);
    } else
      append_range (first, last, iterator_cat { });
  }

  template <typename ForwardIt,
            typename std::enable_if<std::is_assignable<
                                      value_ty&,
                                      decltype (*std::declval<ForwardIt> ())>::value>::type * = nullptr>
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  assign_with_range (ForwardIt first, ForwardIt last, std::forward_iterator_tag) {
    const size_ty count = external_range_length (first, last);
    if (get_capacity () < count) {
      size_ty new_capacity = checked_calculate_new_capacity (count);
      ptr     new_begin    = unchecked_allocate (new_capacity);

      PLUMED_GCH_TRY {
        uninitialized_copy (first, last, new_begin);
      }
      PLUMED_GCH_CATCH (...) {
        deallocate (new_begin, new_capacity);
        PLUMED_GCH_THROW;
      }

      reset_data (new_begin, new_capacity, count);
    } else if (get_size () < count) {
      ForwardIt pivot = copy_n_return_in (first, get_size (), begin_ptr ());
      uninitialized_copy (pivot, last, end_ptr ());
      set_size (count);
    } else {
      erase_range (copy_range (first, last, begin_ptr ()), end_ptr ());
    }
  }

  template <typename InputIt,
            typename std::enable_if<! std::is_assignable<
                                      value_ty&,
                                      decltype (*std::declval<InputIt> ())>::value>::type * = nullptr>
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  assign_with_range (InputIt first, InputIt last, std::input_iterator_tag) {
    using iterator_cat = typename std::iterator_traits<InputIt>::iterator_category;

    // If not assignable then destroy all elements and append.
    erase_all ();
    append_range (first, last, iterator_cat { });
  }

  // Ie. move-if-noexcept.
  struct strong_exception_policy {
  };

  template <typename Policy = void, typename V = value_ty,
            typename std::enable_if<is_explicitly_move_insertable<V>::value
                                    &&  (! std::is_same<Policy, strong_exception_policy>::value
                                         ||  relocate_with_move<V>::value),
                                    bool>::type = true>
  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  uninitialized_move (ptr first, ptr last, ptr d_first)
  noexcept (std::is_nothrow_move_constructible<value_ty>::value) {
    return uninitialized_copy (std::make_move_iterator (first),
                               std::make_move_iterator (last),
                               d_first);
  }

  template <typename Policy = void, typename V = value_ty,
            typename std::enable_if<! is_explicitly_move_insertable<V>::value
                                    ||  (  std::is_same<Policy, strong_exception_policy>::value
                                        &&! relocate_with_move<V>::value),
                                    bool>::type = false>
  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  uninitialized_move (ptr first, ptr last, ptr d_first)
  noexcept (alloc_interface::template is_uninitialized_memcpyable_iterator<ptr>::value) {
    return uninitialized_copy (first, last, d_first);
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  shift_into_uninitialized (ptr pos, size_ty n_shift) {
    // Shift elements over to the right into uninitialized space.
    // Returns the start of the shifted range.
    // Precondition: shift < end_ptr () - pos
    assert (n_shift != 0 && "The value of `n_shift` should not be 0.");

    const ptr original_end = end_ptr ();
    const ptr pivot        = unchecked_prev (original_end, n_shift);

    uninitialized_move (pivot, original_end, original_end);
    increase_size (n_shift);
    return move_right (pos, pivot, original_end);
  }

  template <typename ...Args>
  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  append_element (Args&&... args) {
    if (get_size () < get_capacity ()) {
      return emplace_into_current_end (std::forward<Args> (args)...);
    }
    return emplace_into_reallocation_end (std::forward<Args> (args)...);
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  append_copies (size_ty count, const value_ty& val) {
    if (num_uninitialized () < count) {
      // Reallocate.
      if (get_max_size () - get_size () < count) {
        throw_allocation_size_error ();
      }

      size_ty original_size = get_size ();
      size_ty new_size      = get_size () + count;

      // The check is handled by the if-guard.
      size_ty new_capacity = unchecked_calculate_new_capacity (new_size);
      ptr     new_data_ptr = unchecked_allocate (new_capacity, allocation_end_ptr ());
      ptr     new_last     = unchecked_next (new_data_ptr, original_size);

      PLUMED_GCH_TRY {
        new_last = uninitialized_fill (new_last, unchecked_next (new_last, count), val);
        uninitialized_move (begin_ptr (), end_ptr (), new_data_ptr);
      }
      PLUMED_GCH_CATCH (...) {
        destroy_range (unchecked_next (new_data_ptr, original_size), new_last);
        deallocate (new_data_ptr, new_capacity);
        PLUMED_GCH_THROW;
      }

      reset_data (new_data_ptr, new_capacity, new_size);
      return unchecked_next (new_data_ptr, original_size);
    } else {
      const ptr ret = end_ptr ();
      uninitialized_fill (ret, unchecked_next (ret, count), val);
      increase_size (count);
      return ret;
    }
  }

  template <typename MovePolicy, typename InputIt,
            typename std::enable_if<
              std::is_same<MovePolicy, strong_exception_policy>::value, bool>::type = true>
  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  append_range (InputIt first, InputIt last, std::input_iterator_tag) {
    // Append with a strong exception guarantee.
    size_ty original_size = get_size ();
    for (; ! (first == last); ++first) {
      PLUMED_GCH_TRY {
        append_element (*first);
      }
      PLUMED_GCH_CATCH (...) {
        erase_range (unchecked_next (begin_ptr (), original_size), end_ptr ());
        PLUMED_GCH_THROW;
      }
    }
    return unchecked_next (begin_ptr (), original_size);
  }

  template <typename MovePolicy = void, typename InputIt,
            typename std::enable_if<
              ! std::is_same<MovePolicy, strong_exception_policy>::value, bool>::type = false>
  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  append_range (InputIt first, InputIt last, std::input_iterator_tag) {
    size_ty original_size = get_size ();
    for (; ! (first == last); ++first) {
      append_element (*first);
    }
    return unchecked_next (begin_ptr (), original_size);
  }

  template <typename MovePolicy = void, typename ForwardIt>
  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  append_range (ForwardIt first, ForwardIt last, std::forward_iterator_tag) {
    const size_ty num_insert = external_range_length (first, last);

    if (num_uninitialized () < num_insert) {
      // Reallocate.
      if (get_max_size () - get_size () < num_insert) {
        throw_allocation_size_error ();
      }

      size_ty original_size = get_size ();
      size_ty new_size      = get_size () + num_insert;

      // The check is handled by the if-guard.
      size_ty new_capacity = unchecked_calculate_new_capacity (new_size);
      ptr     new_data_ptr = unchecked_allocate (new_capacity, allocation_end_ptr ());
      ptr     new_last     = unchecked_next (new_data_ptr, original_size);

      PLUMED_GCH_TRY {
        new_last = uninitialized_copy (first, last, new_last);
        uninitialized_move<MovePolicy> (begin_ptr (), end_ptr (), new_data_ptr);
      }
      PLUMED_GCH_CATCH (...) {
        destroy_range (unchecked_next (new_data_ptr, original_size), new_last);
        deallocate (new_data_ptr, new_capacity);
        PLUMED_GCH_THROW;
      }

      reset_data (new_data_ptr, new_capacity, new_size);
      return unchecked_next (new_data_ptr, original_size);
    } else {
      ptr ret = end_ptr ();
      uninitialized_copy (first, last, ret);
      increase_size (num_insert);
      return ret;
    }
  }

  template <typename ...Args>
  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  emplace_at (ptr pos, Args&&... args) {
    assert (get_size () <= get_capacity () && "size was greater than capacity");

    if (get_size () < get_capacity ()) {
      return emplace_into_current (pos, std::forward<Args> (args)...);
    }
    return emplace_into_reallocation (pos, std::forward<Args> (args)...);
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  insert_copies (ptr pos, size_ty count, const value_ty& val) {
    if (0 == count) {
      return pos;
    }

    if (end_ptr () == pos) {
      if (1 == count) {
        return append_element (val);
      }
      return append_copies (count, val);
    }

    if (num_uninitialized () < count) {
      // Reallocate.
      if (get_max_size () - get_size () < count) {
        throw_allocation_size_error ();
      }

      const size_ty offset = internal_range_length (begin_ptr (), pos);

      const size_ty new_size = get_size () + count;

      // The check is handled by the if-guard.
      const size_ty new_capacity = unchecked_calculate_new_capacity (new_size);
      ptr new_data_ptr           = unchecked_allocate (new_capacity, allocation_end_ptr ());
      ptr new_first              = unchecked_next (new_data_ptr, offset);
      ptr new_last               = new_first;

      PLUMED_GCH_TRY {
        uninitialized_fill (new_first, unchecked_next (new_first, count), val);
        unchecked_advance  (new_last, count);

        uninitialized_move (begin_ptr (), pos, new_data_ptr);
        new_first = new_data_ptr;
        uninitialized_move (pos, end_ptr (), new_last);
      }
      PLUMED_GCH_CATCH (...) {
        destroy_range (new_first, new_last);
        deallocate (new_data_ptr, new_capacity);
        PLUMED_GCH_THROW;
      }

      reset_data (new_data_ptr, new_capacity, new_size);
      return unchecked_next (begin_ptr (), offset);
    } else {
      // If we have fewer to insert than tailing elements after `pos`, we shift into
      // uninitialized and then copy over.

      const size_ty tail_size = internal_range_length (pos, end_ptr ());
      if (tail_size < count) {
        // The number inserted is larger than the number after `pos`,
        // so part of the input will be used to construct new elements,
        // and another part of it will assign existing ones.
        // In order:
        //   Construct new elements immediately after end_ptr () using the input.
        //   Move-construct existing elements over to the tail.
        //   Assign existing elements using the input.

        ptr original_end = end_ptr ();

        // Place a portion of the input into the uninitialized section.
        size_ty num_val_tail = count - tail_size;

#ifdef PLUMED_GCH_LIB_IS_CONSTANT_EVALUATED
        if (std::is_constant_evaluated ()) {
          uninitialized_fill (end_ptr (), unchecked_next (end_ptr (), num_val_tail), val);
          increase_size (num_val_tail);

          const heap_temporary tmp (*this, val);

          uninitialized_move (pos, original_end, end_ptr ());
          increase_size (tail_size);

          std::fill_n (pos, tail_size, tmp.get ());

          return pos;
        }
#endif

        uninitialized_fill (end_ptr (), unchecked_next (end_ptr (), num_val_tail), val);
        increase_size (num_val_tail);

        PLUMED_GCH_TRY {
          // We need to handle possible aliasing here.
          const stack_temporary tmp (*this, val);

          // Now, move the tail to the end.
          uninitialized_move (pos, original_end, end_ptr ());
          increase_size (tail_size);

          PLUMED_GCH_TRY
          {
            // Finally, try to copy the rest of the elements over.
            std::fill_n (pos, tail_size, tmp.get ());
          }
          PLUMED_GCH_CATCH (...) {
            // Attempt to roll back and destroy the tail if we fail.
            ptr inserted_end = unchecked_prev (end_ptr (), tail_size);
            move_left (inserted_end, end_ptr (), pos);
            destroy_range (inserted_end, end_ptr ());
            decrease_size (tail_size);
            PLUMED_GCH_THROW;
          }
        }
        PLUMED_GCH_CATCH (...) {
          // Destroy the elements constructed from the input.
          destroy_range (original_end, end_ptr ());
          decrease_size (internal_range_length (original_end, end_ptr ()));
          PLUMED_GCH_THROW;
        }
      } else {
#ifdef PLUMED_GCH_LIB_IS_CONSTANT_EVALUATED
        if (std::is_constant_evaluated ()) {
          const heap_temporary tmp (*this, val);

          ptr inserted_end = shift_into_uninitialized (pos, count);
          std::fill (pos, inserted_end, tmp.get ());

          return pos;
        }
#endif
        const stack_temporary tmp (*this, val);

        ptr inserted_end = shift_into_uninitialized (pos, count);

        // Attempt to copy over the elements.
        // If we fail we'll attempt a full roll-back.
        PLUMED_GCH_TRY {
          std::fill (pos, inserted_end, tmp.get ());
        }
        PLUMED_GCH_CATCH (...) {
          ptr original_end = move_left (inserted_end, end_ptr (), pos);
          destroy_range (original_end, end_ptr ());
          decrease_size (count);
          PLUMED_GCH_THROW;
        }
      }
      return pos;
    }
  }

  template <typename ForwardIt>
  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  insert_range_helper (ptr pos, ForwardIt first, ForwardIt last) {
    assert (! (first == last) && "The range should not be empty.");
    assert (! (end_ptr () == pos) && "`pos` should not be at the end.");

    const size_ty num_insert = external_range_length (first, last);
    if (num_uninitialized () < num_insert) {
      // Reallocate.
      if (get_max_size () - get_size () < num_insert) {
        throw_allocation_size_error ();
      }

      const size_ty offset   = internal_range_length (begin_ptr (), pos);
      const size_ty new_size = get_size () + num_insert;

      // The check is handled by the if-guard.
      const size_ty new_capacity = unchecked_calculate_new_capacity (new_size);
      const ptr     new_data_ptr = unchecked_allocate (new_capacity, allocation_end_ptr ());
      ptr           new_first    = unchecked_next (new_data_ptr, offset);
      ptr           new_last     = new_first;

      PLUMED_GCH_TRY {
        uninitialized_copy (first, last, new_first);
        unchecked_advance  (new_last, num_insert);

        uninitialized_move (begin_ptr (), pos, new_data_ptr);
        new_first = new_data_ptr;
        uninitialized_move (pos, end_ptr (), new_last);
      }
      PLUMED_GCH_CATCH (...) {
        destroy_range (new_first, new_last);
        deallocate (new_data_ptr, new_capacity);
        PLUMED_GCH_THROW;
      }

      reset_data (new_data_ptr, new_capacity, new_size);
      return unchecked_next (begin_ptr (), offset);
    } else {
      // if we have fewer to insert than tailing elements after
      // `pos` we shift into uninitialized and then copy over
      const size_ty tail_size = internal_range_length (pos, end_ptr ());
      if (tail_size < num_insert) {
        // Use the same method as insert_copies.
        ptr original_end = end_ptr ();
        ForwardIt pivot  = unchecked_next (first, tail_size);

        // Place a portion of the input into the uninitialized section.
        uninitialized_copy (pivot, last, end_ptr ());
        increase_size (num_insert - tail_size);

        PLUMED_GCH_TRY {
          // Now move the tail to the end.
          uninitialized_move (pos, original_end, end_ptr ());
          increase_size (tail_size);

          PLUMED_GCH_TRY
          {
            // Finally, try to copy the rest of the elements over.
            copy_range (first, pivot, pos);
          }
          PLUMED_GCH_CATCH (...) {
            // Attempt to roll back and destroy the tail if we fail.
            ptr inserted_end = unchecked_prev (end_ptr (), tail_size);
            move_left (inserted_end, end_ptr (), pos);
            destroy_range (inserted_end, end_ptr ());
            decrease_size (tail_size);
            PLUMED_GCH_THROW;
          }
        }
        PLUMED_GCH_CATCH (...) {
          // If we throw, destroy the first copy we made.
          destroy_range (original_end, end_ptr ());
          decrease_size (internal_range_length (original_end, end_ptr ()));
          PLUMED_GCH_THROW;
        }
      } else {
        shift_into_uninitialized (pos, num_insert);

        // Attempt to copy over the elements.
        // If we fail we'll attempt a full roll-back.
        PLUMED_GCH_TRY {
          copy_range (first, last, pos);
        }
        PLUMED_GCH_CATCH (...) {
          ptr inserted_end = unchecked_next (pos, num_insert);
          ptr original_end = move_left (inserted_end, end_ptr (), pos);
          destroy_range (original_end, end_ptr ());
          decrease_size (num_insert);
          PLUMED_GCH_THROW;
        }
      }
      return pos;
    }
  }

  template <typename InputIt>
  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  insert_range (ptr pos, InputIt first, InputIt last, std::input_iterator_tag) {
    assert (! (first == last) && "The range should not be empty.");

    // Ensure we use this specific overload to give a strong exception guarantee for 1 element.
    if (end_ptr () == pos)
      return append_range (first, last, std::input_iterator_tag { });

    using iterator_cat = typename std::iterator_traits<InputIt>::iterator_category;
    small_vector_base tmp (first, last, iterator_cat { }, allocator_ref ());

    return insert_range_helper (
             pos,
             std::make_move_iterator (tmp.begin_ptr ()),
             std::make_move_iterator (tmp.end_ptr ()));
  }

  template <typename ForwardIt>
  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  insert_range (ptr pos, ForwardIt first, ForwardIt last, std::forward_iterator_tag) {
    if (! (end_ptr () == pos)) {
      return insert_range_helper (pos, first, last);
    }

    if (unchecked_next (first) == last) {
      return append_element (*first);
    }

    using iterator_cat = typename std::iterator_traits<ForwardIt>::iterator_category;
    return append_range (first, last, iterator_cat { });
  }

  template <typename ...Args>
  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  emplace_into_current_end (Args&&... args) {
    construct (end_ptr (), std::forward<Args> (args)...);
    increase_size (1);
    return unchecked_prev (end_ptr ());
  }

  template <typename V = value_ty,
            typename std::enable_if<
              std::is_nothrow_move_constructible<V>::value>::type * = nullptr>
  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  emplace_into_current (ptr pos, value_ty&& val) {
    if (pos == end_ptr ()) {
      return emplace_into_current_end (std::move (val));
    }

    // In the special case of value_ty&& we don't make a copy because behavior is unspecified
    // when it is an internal element. Hence, we'll take the opportunity to optimize and assume
    // that it isn't an internal element.
    shift_into_uninitialized (pos, 1);
    destroy (pos);
    construct (pos, std::move (val));
    return pos;
  }

  template <typename ...Args>
  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  emplace_into_current (ptr pos, Args&&... args) {
    if (pos == end_ptr ()) {
      return emplace_into_current_end (std::forward<Args> (args)...);
    }

#ifdef PLUMED_GCH_LIB_IS_CONSTANT_EVALUATED
    if (std::is_constant_evaluated ()) {
      heap_temporary tmp (*this, std::forward<Args> (args)...);
      shift_into_uninitialized (pos, 1);
      *pos = tmp.release ();
      return pos;
    }
#endif

    // This is necessary because of possible aliasing.
    stack_temporary tmp (*this, std::forward<Args> (args)...);
    shift_into_uninitialized (pos, 1);
    *pos = tmp.release ();
    return pos;
  }

  template <typename ...Args>
  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  emplace_into_reallocation_end (Args&&... args) {
    // Appending; strong exception guarantee.
    if (get_max_size () == get_size ()) {
      throw_allocation_size_error ();
    }

    const size_ty new_size = get_size () + 1;

    // The check is handled by the if-guard.
    const size_ty new_capacity = unchecked_calculate_new_capacity (new_size);
    const ptr     new_data_ptr = unchecked_allocate (new_capacity, allocation_end_ptr ());
    const ptr     emplace_pos  = unchecked_next (new_data_ptr, get_size ());

    PLUMED_GCH_TRY {
      construct (emplace_pos, std::forward<Args> (args)...);
      PLUMED_GCH_TRY
      {
        uninitialized_move<strong_exception_policy> (begin_ptr (), end_ptr (), new_data_ptr);
      }
      PLUMED_GCH_CATCH (...) {
        destroy (emplace_pos);
        PLUMED_GCH_THROW;
      }
    }
    PLUMED_GCH_CATCH (...) {
      deallocate (new_data_ptr, new_capacity);
      PLUMED_GCH_THROW;
    }

    reset_data (new_data_ptr, new_capacity, new_size);
    return emplace_pos;
  }

  template <typename ...Args>
  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  emplace_into_reallocation (ptr pos, Args&&... args) {
    const size_ty offset = internal_range_length (begin_ptr (), pos);
    if (offset == get_size ()) {
      return emplace_into_reallocation_end (std::forward<Args> (args)...);
    }

    if (get_max_size () == get_size ()) {
      throw_allocation_size_error ();
    }

    const size_ty new_size = get_size () + 1;

    // The check is handled by the if-guard.
    const size_ty new_capacity = unchecked_calculate_new_capacity (new_size);
    const ptr     new_data_ptr = unchecked_allocate (new_capacity, allocation_end_ptr ());
    ptr           new_first    = unchecked_next (new_data_ptr, offset);
    ptr           new_last     = new_first;

    PLUMED_GCH_TRY {
      construct (new_first, std::forward<Args> (args)...);
      unchecked_advance (new_last, 1);

      uninitialized_move (begin_ptr (), pos, new_data_ptr);
      new_first = new_data_ptr;
      uninitialized_move (pos, end_ptr (), new_last);
    }
    PLUMED_GCH_CATCH (...) {
      destroy_range (new_first, new_last);
      deallocate (new_data_ptr, new_capacity);
      PLUMED_GCH_THROW;
    }

    reset_data (new_data_ptr, new_capacity, new_size);
    return unchecked_next (begin_ptr (), offset);
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  shrink_to_size (void) {
    if (! has_allocation () || get_size () == get_capacity ()) {
      return begin_ptr ();
    }

    // The rest runs only if allocated.

    size_ty new_capacity;
    ptr     new_data_ptr;

    if (InlineCapacity < get_size ()) {
      new_capacity = get_size ();
      new_data_ptr = unchecked_allocate (new_capacity, allocation_end_ptr ());
    } else {
      // We move to inline storage.
      new_capacity = InlineCapacity;
#ifdef PLUMED_GCH_LIB_IS_CONSTANT_EVALUATED
      if (std::is_constant_evaluated ()) {
        new_data_ptr = alloc_interface::allocate (InlineCapacity);
      } else {
        new_data_ptr = storage_ptr ();
      }
#else
      new_data_ptr = storage_ptr ();
#endif
    }

    uninitialized_move (begin_ptr (), end_ptr (), new_data_ptr);

    destroy_range (begin_ptr (), end_ptr ());
    deallocate (data_ptr (), get_capacity ());

    set_data_ptr (new_data_ptr);
    set_capacity (new_capacity);

    return begin_ptr ();
  }

  template <typename ...ValueT>
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  resize_with (size_ty new_size, const ValueT&... val) {
    // ValueT... should either be value_ty or empty.

    if (new_size == 0) {
      erase_all ();
    }

    if (get_capacity () < new_size) {
      // Reallocate.

      if (get_max_size () < new_size) {
        throw_allocation_size_error ();
      }

      const size_ty original_size = get_size ();

      // The check is handled by the if-guard.
      const size_ty new_capacity = unchecked_calculate_new_capacity (new_size);
      ptr           new_data_ptr = unchecked_allocate (new_capacity, allocation_end_ptr ());
      ptr           new_last     = unchecked_next (new_data_ptr, original_size);

      PLUMED_GCH_TRY {
        new_last = uninitialized_fill (
          new_last,
          unchecked_next (new_data_ptr, new_size),
          val...);

        // Strong exception guarantee.
        uninitialized_move<strong_exception_policy> (begin_ptr (), end_ptr (), new_data_ptr);
      }
      PLUMED_GCH_CATCH (...) {
        destroy_range (unchecked_next (new_data_ptr, original_size), new_last);
        deallocate (new_data_ptr, new_capacity);
        PLUMED_GCH_THROW;
      }

      reset_data (new_data_ptr, new_capacity, new_size);
    } else if (get_size () < new_size) {
      // Construct in the uninitialized section.
      uninitialized_fill (end_ptr (), unchecked_next (begin_ptr (), new_size), val...);
      set_size (new_size);
    } else {
      erase_range (unchecked_next (begin_ptr (), new_size), end_ptr ());
    }

    // Do nothing if the count is the same as the current size.
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  request_capacity (size_ty request) {
    if (request <= get_capacity ()) {
      return;
    }

    size_ty new_capacity = checked_calculate_new_capacity (request);
    ptr     new_begin    = unchecked_allocate (new_capacity);

    PLUMED_GCH_TRY {
      uninitialized_move<strong_exception_policy> (begin_ptr (), end_ptr (), new_begin);
    }
    PLUMED_GCH_CATCH (...) {
      deallocate (new_begin, new_capacity);
      PLUMED_GCH_THROW;
    }

    wipe ();

    set_data_ptr (new_begin);
    set_capacity (new_capacity);
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  erase_at (ptr pos) {
    move_left (unchecked_next (pos), end_ptr (), pos);
    erase_last ();
    return pos;
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  erase_last (void) {
    decrease_size (1);

    // The element located at end_ptr is still alive since the size decreased.
    destroy (end_ptr ());
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  erase_range (ptr first, ptr last) {
    if (! (first == last)) {
      erase_to_end (move_left (last, end_ptr (), first));
    }
    return first;
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  erase_to_end (ptr pos) {
    assert (0 <= (end_ptr () - pos) && "`pos` was in the uninitialized range");
    if (size_ty change = internal_range_length (pos, end_ptr ())) {
      decrease_size (change);
      destroy_range (pos, unchecked_next (pos, change));
    }
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  erase_all (void) {
    ptr curr_end = end_ptr ();
    set_size (0);
    destroy_range (begin_ptr (), curr_end);
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  swap_elements (small_vector_base& other)
  noexcept (std::is_nothrow_move_constructible<value_ty>::value
#ifdef PLUMED_GCH_LIB_IS_SWAPPABLE
            &&  std::is_nothrow_swappable<value_ty>::value
#else
            &&  detail::small_vector_adl::is_nothrow_swappable<value_ty>::value
#endif
           ) {
    assert (get_size () <= other.get_size ());

    const ptr other_tail = std::swap_ranges (begin_ptr (), end_ptr (), other.begin_ptr ());
    uninitialized_move (other_tail, other.end_ptr (), end_ptr ());
    destroy_range (other_tail, other.end_ptr ());

    swap_size (other);
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  swap_default (small_vector_base& other)
  noexcept (std::is_nothrow_move_constructible<value_ty>::value
#ifdef PLUMED_GCH_LIB_IS_SWAPPABLE
            &&  std::is_nothrow_swappable<value_ty>::value
#else
            &&  detail::small_vector_adl::is_nothrow_swappable<value_ty>::value
#endif
           ) {
    // This function is used when:
    //   We are using the standard allocator.
    //   The allocators propagate and are equal.
    //   The allocators are always equal.
    //   The allocators do not propagate and are equal.
    //   The allocators propagate and are not equal.

    // Not handled:
    //   The allocators do not propagate and are not equal.

    assert (get_capacity () <= other.get_capacity ());

    if (has_allocation ()) { // Implies that `other` also has an allocation.
      swap_allocation (other);
    } else if (other.has_allocation ()) {
      // Note: This will never be constant evaluated because both are always allocated.
      uninitialized_move (begin_ptr (), end_ptr (), other.storage_ptr ());
      destroy_range (begin_ptr (), end_ptr ());

      set_data_ptr (other.data_ptr ());
      set_capacity (other.get_capacity ());

      other.set_data_ptr (other.storage_ptr ());
      other.set_capacity (InlineCapacity);

      swap_size (other);
    } else if (get_size () < other.get_size ()) {
      swap_elements (other);
    } else {
      other.swap_elements (*this);
    }

    alloc_interface::swap (other);
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  swap_unequal_no_propagate (small_vector_base& other) {
    assert (get_capacity () <= other.get_capacity ());

    if (get_capacity () < other.get_size ()) {
      // Reallocation required.
      // We should always be able to reuse the allocation of `other`.
      const size_ty new_capacity = unchecked_calculate_new_capacity (other.get_size ());
      const ptr     new_data_ptr = unchecked_allocate (new_capacity, end_ptr ());

      PLUMED_GCH_TRY {
        uninitialized_move (other.begin_ptr (), other.end_ptr (), new_data_ptr);
        PLUMED_GCH_TRY
        {
          destroy_range (
            std::move (begin_ptr (), end_ptr (), other.begin_ptr ()),
            other.end_ptr ());
        }
        PLUMED_GCH_CATCH (...) {
          destroy_range (new_data_ptr, unchecked_next (new_data_ptr, other.get_size ()));
          PLUMED_GCH_THROW;
        }
      }
      PLUMED_GCH_CATCH (...) {
        deallocate (new_data_ptr, new_capacity);
        PLUMED_GCH_THROW;
      }

      destroy_range (begin_ptr (), end_ptr ());
      if (has_allocation ()) {
        deallocate (data_ptr (), get_capacity ());
      }

      set_data_ptr (new_data_ptr);
      set_capacity (new_capacity);
      swap_size (other);
    } else if (get_size () < other.get_size ()) {
      swap_elements (other);
    } else {
      other.swap_elements (*this);
    }

    // This should have no effect.
    alloc_interface::swap (other);
  }

  template <typename A = alloc_ty,
            typename std::enable_if<allocations_are_swappable<A>::value
                                    &&  InlineCapacity == 0>::type * = nullptr>
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  swap (small_vector_base& other) noexcept {
    swap_allocation (other);
    alloc_interface::swap (other);
  }

  template <typename A = alloc_ty,
            typename std::enable_if<allocations_are_swappable<A>::value
                                    &&  InlineCapacity != 0>::type * = nullptr>
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  swap (small_vector_base& other)
  noexcept (std::is_nothrow_move_constructible<value_ty>::value
#ifdef PLUMED_GCH_LIB_IS_SWAPPABLE
            &&  std::is_nothrow_swappable<value_ty>::value
#else
            &&  detail::small_vector_adl::is_nothrow_swappable<value_ty>::value
#endif
                                         ) {
    if (get_capacity () < other.get_capacity ()) {
      swap_default (other);
    } else {
      other.swap_default (*this);
    }
  }

  template <typename A = alloc_ty,
            typename std::enable_if<! allocations_are_swappable<A>::value>::type * = nullptr>
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  swap (small_vector_base& other) {
    if (get_capacity () < other.get_capacity ()) {
      if (other.allocator_ref () == allocator_ref ()) {
        swap_default (other);
      } else {
        swap_unequal_no_propagate (other);
      }
    } else {
      if (other.allocator_ref () == allocator_ref ()) {
        other.swap_default (*this);
      } else {
        other.swap_unequal_no_propagate (*this);
      }
    }
  }

#ifdef __GLIBCXX__

  // These are compatibility fixes for libstdc++ because std::copy doesn't work for
  // `move_iterator`s when constant evaluated.

  template <typename InputIt>
  static PLUMED_GCH_CPP20_CONSTEXPR
  InputIt
  unmove_iterator (InputIt it) {
    return it;
  }

  template <typename InputIt>
  static PLUMED_GCH_CPP20_CONSTEXPR
  auto
  unmove_iterator (std::move_iterator<InputIt> it)
  -> decltype (unmove_iterator (it.base ())) {
    return unmove_iterator (it.base ());
  }

  template <typename InputIt>
  static PLUMED_GCH_CPP20_CONSTEXPR
  auto
  unmove_iterator (std::reverse_iterator<InputIt> it)
  -> std::reverse_iterator<decltype (unmove_iterator (it.base ()))> {
    return std::reverse_iterator<decltype (unmove_iterator (it.base ()))> (
      unmove_iterator (it.base ()));
                                }

#endif

  template <typename InputIt>
  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  copy_range (InputIt first, InputIt last, ptr dest) {
#if defined (PLUMED_GCH_LIB_IS_CONSTANT_EVALUATED) && defined (__GLIBCXX__)
    if (    std::is_constant_evaluated ()
            &&! std::is_same<decltype (unmove_iterator (std::declval<InputIt> ())),
                             InputIt>::value) {
      return std::move (unmove_iterator (first), unmove_iterator (last), dest);
    }
#endif

    return std::copy (first, last, dest);
  }

  template <typename InputIt,
            typename std::enable_if<
              is_memcpyable_iterator<InputIt>::value>::type * = nullptr>
  PLUMED_GCH_CPP20_CONSTEXPR
  InputIt
  copy_n_return_in (InputIt first, size_ty count, ptr dest) noexcept {
#ifdef PLUMED_GCH_LIB_IS_CONSTANT_EVALUATED
    if (std::is_constant_evaluated ()) {
      std::copy_n (first, count, dest);
      return unchecked_next (first, count);
    }
#endif

    if (count != 0) {
      std::memcpy (to_address (dest), to_address (first), count * sizeof (value_ty));
    }
    // Note: The unsafe cast here should be proven to be safe in the caller function.
    return unchecked_next (first, count);
  }

  template <typename InputIt,
            typename std::enable_if<
              is_memcpyable_iterator<InputIt>::value>::type * = nullptr>
  PLUMED_GCH_CPP20_CONSTEXPR
  std::move_iterator<InputIt>
  copy_n_return_in (std::move_iterator<InputIt> first, size_ty count, ptr dest) noexcept {
    return std::move_iterator<InputIt> (copy_n_return_in (first.base (), count, dest));
  }

  template <typename RandomIt,
            typename std::enable_if<
              ! is_memcpyable_iterator<RandomIt>::value
              &&  std::is_base_of<std::random_access_iterator_tag,
                                  typename std::iterator_traits<RandomIt>::iterator_category>::value
              >::type * = nullptr>
  PLUMED_GCH_CPP20_CONSTEXPR
  RandomIt
  copy_n_return_in (RandomIt first, size_ty count, ptr dest) {
#if defined (PLUMED_GCH_LIB_IS_CONSTANT_EVALUATED) && defined (__GLIBCXX__)
    if (    std::is_constant_evaluated ()
            &&! std::is_same<decltype (unmove_iterator (std::declval<RandomIt> ())),
                             RandomIt>::value) {
      auto bfirst = unmove_iterator (first);
      auto blast  = unchecked_next (bfirst, count);
      std::move (bfirst, blast, dest);
      return unchecked_next (first, count);
    }
#endif

    std::copy_n (first, count, dest);
    // Note: This unsafe cast should be proven safe in the caller function.
    return unchecked_next (first, count);
  }

  template <typename InputIt,
            typename std::enable_if<
              ! is_memcpyable_iterator<InputIt>::value
              &&! std::is_base_of<std::random_access_iterator_tag,
                                  typename std::iterator_traits<InputIt>::iterator_category>::value
              >::type * = nullptr>
  PLUMED_GCH_CPP20_CONSTEXPR
  InputIt
  copy_n_return_in (InputIt first, size_ty count, ptr dest) {

    for (; count != 0; --count, static_cast<void> (++dest), static_cast<void> (++first)) {
      *dest = *first;
    }
    return first;
  }

  template <typename V = value_ty,
            typename std::enable_if<is_memcpyable<V>::value>::type * = nullptr>
  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  move_left (ptr first, ptr last, ptr d_first) {
    // Shift initialized elements to the left.

#ifdef PLUMED_GCH_LIB_IS_CONSTANT_EVALUATED
    if (std::is_constant_evaluated ()) {
      return std::move (first, last, d_first);
    }
#endif

    const size_ty num_moved = internal_range_length (first, last);
    if (num_moved != 0) {
      std::memmove (to_address (d_first), to_address (first), num_moved * sizeof (value_ty));
    }
    return unchecked_next (d_first, num_moved);
  }

  template <typename V = value_ty,
            typename std::enable_if<! is_memcpyable<V>::value>::type * = nullptr>
  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  move_left (ptr first, ptr last, ptr d_first) {
    // Shift initialized elements to the left.
    return std::move (first, last, d_first);
  }

  template <typename V = value_ty,
            typename std::enable_if<is_memcpyable<V>::value, bool>::type = true>
  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  move_right (ptr first, ptr last, ptr d_last) {
    // Move initialized elements to the right.

#ifdef PLUMED_GCH_LIB_IS_CONSTANT_EVALUATED
    if (std::is_constant_evaluated ()) {
      return std::move_backward (first, last, d_last);
    }
#endif

    const size_ty num_moved = internal_range_length (first, last);
    const ptr     dest      = unchecked_prev (d_last, num_moved);
    if (num_moved != 0) {
      std::memmove (to_address (dest), to_address (first), num_moved * sizeof (value_ty));
    }
    return dest;
  }

  template <typename V = value_ty,
            typename std::enable_if<! is_memcpyable<V>::value, bool>::type = false>
  PLUMED_GCH_CPP20_CONSTEXPR
  ptr
  move_right (ptr first, ptr last, ptr d_last) {
    // move initialized elements to the right
    // n should not be 0
    return std::move_backward (first, last, d_last);
  }

public:
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  set_default (void) {
    set_to_inline_storage ();
    set_size (0);
  }

  PLUMED_GCH_NODISCARD PLUMED_GCH_CPP14_CONSTEXPR
  ptr
  data_ptr (void) noexcept {
    return m_data.data_ptr ();
  }

  PLUMED_GCH_NODISCARD constexpr
  cptr
  data_ptr (void) const noexcept {
    return m_data.data_ptr ();
  }

  PLUMED_GCH_NODISCARD constexpr
  size_ty
  get_capacity (void) const noexcept {
    return m_data.capacity ();
  }

  PLUMED_GCH_NODISCARD constexpr
  size_ty
  get_size (void) const noexcept {
    return m_data.size ();
  }

  PLUMED_GCH_NODISCARD constexpr
  size_ty
  num_uninitialized (void) const noexcept {
    return get_capacity () - get_size ();
  }

  PLUMED_GCH_NODISCARD PLUMED_GCH_CPP14_CONSTEXPR
  ptr
  begin_ptr (void) noexcept {
    return data_ptr ();
  }

  PLUMED_GCH_NODISCARD
  constexpr
  cptr
  begin_ptr (void) const noexcept {
    return data_ptr ();
  }

  PLUMED_GCH_NODISCARD PLUMED_GCH_CPP14_CONSTEXPR
  ptr
  end_ptr (void) noexcept {
    return unchecked_next (begin_ptr (), get_size ());
  }

  PLUMED_GCH_NODISCARD constexpr
  cptr
  end_ptr (void) const noexcept {
    return unchecked_next (begin_ptr (), get_size ());
  }

  PLUMED_GCH_NODISCARD PLUMED_GCH_CPP14_CONSTEXPR
  ptr
  allocation_end_ptr (void) noexcept {
    return unchecked_next (begin_ptr (), get_capacity ());
  }

  PLUMED_GCH_NODISCARD constexpr
  cptr
  allocation_end_ptr (void) const noexcept {
    return unchecked_next (begin_ptr (), get_capacity ());
  }

  PLUMED_GCH_NODISCARD constexpr
  alloc_ty
  copy_allocator (void) const noexcept {
    return alloc_ty (allocator_ref ());
  }

  PLUMED_GCH_NODISCARD PLUMED_GCH_CPP14_CONSTEXPR
  ptr
  storage_ptr (void) noexcept {
    return m_data.storage ();
  }

  PLUMED_GCH_NODISCARD constexpr
  cptr
  storage_ptr (void) const noexcept {
    return m_data.storage ();
  }

  PLUMED_GCH_NODISCARD constexpr
  bool
  has_allocation (void) const noexcept {
#ifdef PLUMED_GCH_LIB_IS_CONSTANT_EVALUATED
    if (std::is_constant_evaluated ()) {
      return true;
    }
#endif
    return InlineCapacity < get_capacity ();
  }

  PLUMED_GCH_NODISCARD constexpr
  bool
  is_inlinable (void) const noexcept {
    return get_size () <= InlineCapacity;
  }

private:
  small_vector_data<ptr, size_type, value_ty, InlineCapacity> m_data;
};

} // namespace gch::detail

template <typename T, unsigned InlineCapacity, typename Allocator>
#ifdef PLUMED_GCH_LIB_CONCEPTS
requires concepts::small_vector::AllocatorFor<Allocator, T>
#endif
class small_vector
: private detail::small_vector_base<Allocator, InlineCapacity> {
  using base = detail::small_vector_base<Allocator, InlineCapacity>;

public:
  static_assert (std::is_same<T, typename Allocator::value_type>::value,
                 "`Allocator::value_type` must be the same as `T`.");

  template <typename SameT, unsigned DifferentInlineCapacity, typename SameAllocator>
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires concepts::small_vector::AllocatorFor<SameAllocator, SameT>
#endif
  friend class small_vector;

  using value_type             = T;
  using allocator_type         = Allocator;
  using size_type              = typename base::size_type;
  using difference_type        = typename base::difference_type;
  using reference              =       value_type&;
  using const_reference        = const value_type&;
  using pointer                = typename std::allocator_traits<allocator_type>::pointer;
  using const_pointer          = typename std::allocator_traits<allocator_type>::const_pointer;

  using iterator               = small_vector_iterator<pointer, difference_type>;
  using const_iterator         = small_vector_iterator<const_pointer, difference_type>;
  using reverse_iterator       = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  static_assert (InlineCapacity <= (std::numeric_limits<size_type>::max) (),
                 "InlineCapacity must be less than or equal to the maximum value of size_type.");

  static constexpr
  unsigned
  inline_capacity_v = InlineCapacity;

#ifdef PLUMED_GCH_LIB_CONCEPTS

private:
  static constexpr
  bool
  Destructible = concepts::small_vector::Destructible<value_type>;

  static constexpr
  bool
  MoveAssignable = concepts::small_vector::MoveAssignable<value_type>;

  static constexpr
  bool
  CopyAssignable = concepts::small_vector::CopyAssignable<value_type>;

  static constexpr
  bool
  MoveConstructible = concepts::small_vector::MoveConstructible<value_type>;

  static constexpr
  bool
  CopyConstructible = concepts::small_vector::CopyConstructible<value_type>;

  static constexpr
  bool
  Swappable = concepts::small_vector::Swappable<value_type>;

  static constexpr
  bool
  DefaultInsertable = concepts::small_vector::DefaultInsertable<value_type, small_vector,
  allocator_type>;

  static constexpr
  bool
  MoveInsertable = concepts::small_vector::MoveInsertable<value_type, small_vector,
  allocator_type>;

  static constexpr
  bool
  CopyInsertable = concepts::small_vector::CopyInsertable<value_type, small_vector,
  allocator_type>;

  static constexpr
  bool
  Erasable = concepts::small_vector::Erasable<value_type, small_vector, allocator_type>;

  template <typename ...Args>
  struct EmplaceConstructible {
    static constexpr
    bool
    value = concepts::small_vector::EmplaceConstructible<value_type, small_vector,
        allocator_type, Args...>;
                                                        };

public:

#endif

  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector (void)
  noexcept (noexcept (allocator_type ()))
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires concepts::DefaultConstructible<allocator_type>
#endif
    = default;

  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector (const small_vector& other)
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires CopyInsertable
#endif
:
  base (base::bypass, other)
  { }

  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector (small_vector&& other)
  noexcept (std::is_nothrow_move_constructible<value_type>::value || InlineCapacity == 0)
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires MoveInsertable
#endif
:
  base (base::bypass, std::move (other))
  { }

  PLUMED_GCH_CPP20_CONSTEXPR explicit
  small_vector (const allocator_type& alloc) noexcept
    : base (alloc)
  { }

  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector (const small_vector& other, const allocator_type& alloc)
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires CopyInsertable
#endif
:
  base (base::bypass, other, alloc)
  { }

  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector (small_vector&& other, const allocator_type& alloc)
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires MoveInsertable
#endif
:
  base (base::bypass, std::move (other), alloc)
  { }

  PLUMED_GCH_CPP20_CONSTEXPR explicit
  small_vector (size_type count, const allocator_type& alloc = allocator_type ())
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires DefaultInsertable
#endif
:
  base (count, alloc)
  { }

  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector (size_type count, const_reference value,
                const allocator_type& alloc = allocator_type ())
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires CopyInsertable
#endif
:
  base (count, value, alloc)
  { }

#ifdef PLUMED_GCH_LIB_CONCEPTS
  template <typename Generator>
  requires std::invocable<Generator&>
  &&  EmplaceConstructible<std::invoke_result_t<Generator&>>::value
#else
  template <typename Generator,
            typename std::enable_if<
              ! std::is_convertible<Generator, const_reference>::value>::type * = nullptr>
#endif
  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector (size_type count, Generator g, const allocator_type& alloc = allocator_type ())
    : base (count, g, alloc)
  { }

#ifdef PLUMED_GCH_LIB_CONCEPTS
  template <std::input_iterator InputIt>
  requires EmplaceConstructible<std::iter_reference_t<InputIt>>::value
  &&  (std::forward_iterator<InputIt> || MoveInsertable)
#else
  template <typename InputIt,
            typename std::enable_if<
              std::is_base_of<
                std::input_iterator_tag,
                typename std::iterator_traits<InputIt>::iterator_category>::value
              >::type * = nullptr>
#endif
  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector (InputIt first, InputIt last, const allocator_type& alloc = allocator_type ())
    : base (first, last, typename std::iterator_traits<InputIt>::iterator_category { }, alloc)
  { }

  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector (std::initializer_list<value_type> init,
                const allocator_type& alloc = allocator_type ())
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires EmplaceConstructible<const_reference>::value
#endif
:
  small_vector (init.begin (), init.end (), alloc)
  { }

  template <unsigned I>
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires CopyInsertable
#endif
  PLUMED_GCH_CPP20_CONSTEXPR explicit
  small_vector (const small_vector<T, I, Allocator>& other)
    : base (base::bypass, other)
  { }

  template <unsigned I>
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires MoveInsertable
#endif
  PLUMED_GCH_CPP20_CONSTEXPR explicit
  small_vector (small_vector<T, I, Allocator>&& other)
  noexcept (std::is_nothrow_move_constructible<value_type>::value && I < InlineCapacity)
              : base (base::bypass, std::move (other))
  { }

  template <unsigned I>
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires CopyInsertable
#endif
  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector (const small_vector<T, I, Allocator>& other, const allocator_type& alloc)
    : base (base::bypass, other, alloc)
  { }

  template <unsigned I>
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires MoveInsertable
#endif
  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector (small_vector<T, I, Allocator>&& other, const allocator_type& alloc)
    : base (base::bypass, std::move (other), alloc)
  { }

  PLUMED_GCH_CPP20_CONSTEXPR
  ~small_vector (void)
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires Erasable
#endif
    = default;

  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector&
  operator= (const small_vector& other)
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires CopyInsertable && CopyAssignable
#endif
  {
    assign (other);
    return *this;
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector&
  operator= (small_vector&& other)
  noexcept (  (  std::is_same<std::allocator<value_type>, Allocator>::value
                 ||  std::allocator_traits<Allocator>::propagate_on_container_move_assignment::value
#ifdef PLUMED_GCH_LIB_IS_ALWAYS_EQUAL
                 ||  std::allocator_traits<Allocator>::is_always_equal::value
#endif
              )
              &&  (  (  std::is_nothrow_move_assignable<value_type>::value
                        &&  std::is_nothrow_move_constructible<value_type>::value
                     )
                     ||  InlineCapacity == 0
                  )
           )
#ifdef PLUMED_GCH_LIB_CONCEPTS
  // Note: The standard says here that
  // std::allocator_traits<allocator_type>::propagate_on_container_move_assignment == false
  // implies MoveInsertable && MoveAssignable, but since we have inline storage we must always
  // require moves [tab:container.alloc.req].
  requires MoveInsertable && MoveAssignable
#endif
  {
    assign (std::move (other));
    return *this;
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector&
  operator= (std::initializer_list<value_type> ilist)
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires CopyInsertable && CopyAssignable
#endif
  {
    assign (ilist);
    return *this;
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  assign (size_type count, const_reference value)
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires CopyInsertable && CopyAssignable
#endif
  {
    base::assign_with_copies (count, value);
  }

#ifdef PLUMED_GCH_LIB_CONCEPTS
  template <std::input_iterator InputIt>
  requires EmplaceConstructible<std::iter_reference_t<InputIt>>::value
  &&  (std::forward_iterator<InputIt> || MoveInsertable)
#else
  template <typename InputIt,
            typename std::enable_if<std::is_base_of<
                                      std::input_iterator_tag,
                                      typename std::iterator_traits<InputIt>::iterator_category
                                      >::value>::type * = nullptr>
#endif
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  assign (InputIt first, InputIt last) {
    using iterator_cat = typename std::iterator_traits<InputIt>::iterator_category;
    base::assign_with_range (first, last, iterator_cat { });
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  assign (std::initializer_list<value_type> ilist)
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires EmplaceConstructible<const_reference>::value
#endif
  {
    assign (ilist.begin (), ilist.end ());
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  assign (const small_vector& other)
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires CopyInsertable && CopyAssignable
#endif
  {
    if (&other != this) {
      base::copy_assign (other);
    }
  }

  template <unsigned I>
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires CopyInsertable && CopyAssignable
#endif
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  assign (const small_vector<T, I, Allocator>& other) {
    base::copy_assign (other);
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  assign (small_vector&& other)
  noexcept (  (  std::is_same<std::allocator<value_type>, Allocator>::value
                 ||  std::allocator_traits<Allocator>::propagate_on_container_move_assignment::value
#ifdef PLUMED_GCH_LIB_IS_ALWAYS_EQUAL
                 ||  std::allocator_traits<Allocator>::is_always_equal::value
#endif
              )
              &&  (  (  std::is_nothrow_move_assignable<value_type>::value
                        &&  std::is_nothrow_move_constructible<value_type>::value
                     )
                     ||  InlineCapacity == 0
                  )
           )
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires MoveInsertable && MoveAssignable
#endif
  {
    if (&other != this) {
      base::move_assign (std::move (other));
    }
  }

  template <unsigned I>
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires MoveInsertable && MoveAssignable
#endif
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  assign (small_vector<T, I, Allocator>&& other)
  noexcept (  I <= InlineCapacity
                 &&  (  std::is_same<std::allocator<value_type>, Allocator>::value
                        ||  std::allocator_traits<Allocator>::propagate_on_container_move_assignment::value
#ifdef PLUMED_GCH_LIB_IS_ALWAYS_EQUAL
                        ||  std::allocator_traits<Allocator>::is_always_equal::value
#endif
                                                 )
                 &&  std::is_nothrow_move_assignable<value_type>::value
                 &&  std::is_nothrow_move_constructible<value_type>::value
  ) {
    base::move_assign (std::move (other));
  }

#ifndef PLUMED_GCH_LIB_CONCEPTS
  template <typename ValueType = value_type,
            typename std::enable_if<
              (  std::is_move_constructible<ValueType>::value
                 &&  std::is_move_assignable<ValueType>::value
#ifdef PLUMED_GCH_LIB_IS_SWAPPABLE
                 &&  std::is_swappable<ValueType>::value
#endif
                                      )
              ||  (  (  std::is_same<std::allocator<value_type>, Allocator>::value
                        ||  std::allocator_traits<Allocator>::propagate_on_container_swap::value
#ifdef PLUMED_GCH_LIB_IS_ALWAYS_EQUAL
                        ||  std::allocator_traits<Allocator>::is_always_equal::value
#endif
                                                 )
                     &&  InlineCapacity == 0
                                                 )
              >::type * = nullptr>
#endif
  PLUMED_GCH_CPP20_CONSTEXPR
  void
  swap (small_vector& other)
  noexcept (  (  std::is_same<std::allocator<value_type>, Allocator>::value
                 ||  std::allocator_traits<Allocator>::propagate_on_container_swap::value
#ifdef PLUMED_GCH_LIB_IS_ALWAYS_EQUAL
                 ||  std::allocator_traits<Allocator>::is_always_equal::value
#endif
                                          )
              &&  (  (  std::is_nothrow_move_constructible<value_type>::value
                        &&  std::is_nothrow_move_assignable<value_type>::value
#ifdef PLUMED_GCH_LIB_IS_SWAPPABLE
                        &&  std::is_nothrow_swappable<value_type>::value
#else
                        &&  detail::small_vector_adl::is_nothrow_swappable<value_type>::value
#endif
                                                     )
                     ||  InlineCapacity == 0
                                                           )
                                                          )
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires (MoveInsertable && MoveAssignable && Swappable)
  ||  (  (  std::is_same<std::allocator<value_type>, Allocator>::value
            ||  std::allocator_traits<Allocator>::propagate_on_container_swap::value
#ifdef PLUMED_GCH_LIB_IS_ALWAYS_EQUAL
            ||  std::allocator_traits<Allocator>::is_always_equal::value
#endif
                                     )
         &&  InlineCapacity == 0
                                     )
#endif
  {
    base::swap (other);
  }

  PLUMED_GCH_CPP14_CONSTEXPR
  iterator
  begin (void) noexcept {
    return iterator { base::begin_ptr () };
  }

  constexpr
  const_iterator
  begin (void) const noexcept {
    return const_iterator { base::begin_ptr () };
  }

  constexpr
  const_iterator
  cbegin (void) const noexcept {
    return begin ();
  }

  PLUMED_GCH_CPP14_CONSTEXPR
  iterator
  end (void) noexcept {
    return iterator { base::end_ptr () };
  }

  constexpr
  const_iterator
  end (void) const noexcept {
    return const_iterator { base::end_ptr () };
  }

  constexpr
  const_iterator
  cend (void) const noexcept {
    return end ();
  }

  PLUMED_GCH_CPP14_CONSTEXPR
  reverse_iterator
  rbegin (void) noexcept {
    return reverse_iterator { end () };
  }

  constexpr
  const_reverse_iterator
  rbegin (void) const noexcept {
    return const_reverse_iterator { end () };
  }

  constexpr
  const_reverse_iterator
  crbegin (void) const noexcept {
    return rbegin ();
  }

  PLUMED_GCH_CPP14_CONSTEXPR
  reverse_iterator
  rend (void) noexcept {
    return reverse_iterator { begin () };
  }

  constexpr
  const_reverse_iterator
  rend (void) const noexcept {
    return const_reverse_iterator { begin () };
  }

  constexpr
  const_reverse_iterator
  crend (void) const noexcept {
    return rend ();
  }

  PLUMED_GCH_CPP14_CONSTEXPR
  reference
  at (size_type pos) {
    if (size () <= pos) {
      base::throw_index_error ();
    }
    return begin ()[static_cast<difference_type> (pos)];
  }

  PLUMED_GCH_CPP14_CONSTEXPR
  const_reference
  at (size_type pos) const {
    if (size () <= pos) {
      base::throw_index_error ();
    }
    return begin ()[static_cast<difference_type> (pos)];
  }

  PLUMED_GCH_CPP14_CONSTEXPR
  reference
  operator[] (size_type pos) {
#ifdef _GLIBCXX_DEBUG
    if (size () <= pos) {
      base::throw_index_error ();
    }
#endif
    return begin ()[static_cast<difference_type> (pos)];
  }

  constexpr
  const_reference
  operator[] (size_type pos) const {
#ifdef _GLIBCXX_DEBUG
    if (size () <= pos) {
      base::throw_index_error ();
    }
#endif
    return begin ()[static_cast<difference_type> (pos)];
  }

  PLUMED_GCH_CPP14_CONSTEXPR
  reference
  front (void) {
    return (*this)[0];
  }

  constexpr
  const_reference
  front (void) const {
    return (*this)[0];
  }

  PLUMED_GCH_CPP14_CONSTEXPR
  reference
  back (void) {
    return (*this)[size () - 1];
  }

  constexpr
  const_reference
  back (void) const {
    return (*this)[size () - 1];
  }

  PLUMED_GCH_CPP14_CONSTEXPR
  pointer
  data (void) noexcept {
    return base::begin_ptr ();
  }

  constexpr
  const_pointer
  data (void) const noexcept {
    return base::begin_ptr ();
  }

  constexpr
  size_type
  size (void) const noexcept {
    return static_cast<size_type> (base::get_size ());
  }

  PLUMED_GCH_NODISCARD constexpr
  bool
  empty (void) const noexcept {
    return size () == 0;
  }

  PLUMED_GCH_CPP14_CONSTEXPR
  size_type
  max_size (void) const noexcept {
    return static_cast<size_type> (base::get_max_size ());
  }

  constexpr
  size_type
  capacity (void) const noexcept {
    return static_cast<size_type> (base::get_capacity ());
  }

  constexpr
  allocator_type
  get_allocator (void) const noexcept {
    return base::copy_allocator ();
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  iterator
  insert (const_iterator pos, const_reference value)
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires CopyInsertable && CopyAssignable
#endif
  {
    return emplace (pos, value);
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  iterator
  insert (const_iterator pos, value_type&& value)
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires MoveInsertable && MoveAssignable
#endif
  {
    return emplace (pos, std::move (value));
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  iterator
  insert (const_iterator pos, size_type count, const_reference value)
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires CopyInsertable && CopyAssignable
#endif
  {
    return iterator (base::insert_copies (base::ptr_cast (pos), count, value));
  }

  // Note: Unlike std::vector, this does not require MoveConstructible because we
  //       don't use std::rotate (as was the reason for the change in C++17).
  //       Relevant: https://cplusplus.github.io/LWG/issue2266).
#ifdef PLUMED_GCH_LIB_CONCEPTS
  template <std::input_iterator InputIt>
  requires EmplaceConstructible<std::iter_reference_t<InputIt>>::value
  &&  MoveInsertable
  &&  MoveAssignable
#else
  template <typename InputIt,
            typename std::enable_if<std::is_base_of<
                                      std::input_iterator_tag,
                                      typename std::iterator_traits<InputIt>::iterator_category
                                      >::value>::type * = nullptr>
#endif
  PLUMED_GCH_CPP20_CONSTEXPR
  iterator
  insert (const_iterator pos, InputIt first, InputIt last) {
    if (first == last) {
      return iterator (base::ptr_cast (pos));
    }

    using iterator_cat = typename std::iterator_traits<InputIt>::iterator_category;
    return iterator (base::insert_range (base::ptr_cast (pos), first, last, iterator_cat { }));
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  iterator
  insert (const_iterator pos, std::initializer_list<value_type> ilist)
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires EmplaceConstructible<const_reference>::value
  &&  MoveInsertable
  &&  MoveAssignable
#endif
  {
    return insert (pos, ilist.begin (), ilist.end ());
  }

  template <typename ...Args>
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires EmplaceConstructible<Args...>::value
  &&  MoveInsertable
  &&  MoveAssignable
#endif
  PLUMED_GCH_CPP20_CONSTEXPR
  iterator
  emplace (const_iterator pos, Args&&... args) {
    return iterator (base::emplace_at (base::ptr_cast (pos), std::forward<Args> (args)...));
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  iterator
  erase (const_iterator pos)
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires MoveAssignable && Erasable
#endif
  {
    assert (0 <= (pos    - begin ()) && "`pos` is out of bounds (before `begin ()`)."   );
    assert (0 <  (end () - pos)      && "`pos` is out of bounds (at or after `end ()`).");

    return iterator (base::erase_at (base::ptr_cast (pos)));
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  iterator
  erase (const_iterator first, const_iterator last)
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires MoveAssignable && Erasable
#endif
  {
    assert (0 <= (last   - first)    && "Invalid range.");
    assert (0 <= (first  - begin ()) && "`first` is out of bounds (before `begin ()`)."  );
    assert (0 <= (end () - last)     && "`last` is out of bounds (after `end ()`).");

    return iterator (base::erase_range (base::ptr_cast (first), base::ptr_cast (last)));
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  push_back (const_reference value)
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires CopyInsertable
#endif
  {
    emplace_back (value);
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  push_back (value_type&& value)
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires MoveInsertable
#endif
  {
    emplace_back (std::move (value));
  }

  template <typename ...Args>
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires EmplaceConstructible<Args...>::value && MoveInsertable
#endif
  PLUMED_GCH_CPP20_CONSTEXPR
  reference
  emplace_back (Args&&... args) {
    return *base::append_element (std::forward<Args> (args)...);
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  pop_back (void)
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires Erasable
#endif
  {
    assert (! empty () && "`pop_back ()` called on an empty `small_vector`.");
    base::erase_last ();
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  reserve (size_type new_capacity)
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires MoveInsertable
#endif
  {
    base::request_capacity (new_capacity);
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  shrink_to_fit (void)
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires MoveInsertable
#endif
  {
    base::shrink_to_size ();
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  clear (void) noexcept
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires Erasable
#endif
  {
    base::erase_all ();
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  resize (size_type count)
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires MoveInsertable && DefaultInsertable
#endif
  {
    base::resize_with (count);
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  void
  resize (size_type count, const_reference value)
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires CopyInsertable
#endif
  {
    base::resize_with (count, value);
  }

  PLUMED_GCH_NODISCARD constexpr
  bool
  inlined (void) const noexcept {
    return ! base::has_allocation ();
  }

  PLUMED_GCH_NODISCARD constexpr
  bool
  inlinable (void) const noexcept {
    return base::is_inlinable ();
  }

  PLUMED_GCH_NODISCARD
  static PLUMED_GCH_CONSTEVAL
  size_type
  inline_capacity (void) noexcept {
    return static_cast<size_type> (inline_capacity_v);
  }

#ifdef PLUMED_GCH_LIB_CONCEPTS
  template <std::input_iterator InputIt>
  requires EmplaceConstructible<std::iter_reference_t<InputIt>>::value
  &&  MoveInsertable
#else
  template <typename InputIt,
            typename std::enable_if<std::is_base_of<
                                      std::input_iterator_tag,
                                      typename std::iterator_traits<InputIt>::iterator_category
                                      >::value>::type * = nullptr>
#endif
  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector&
  append (InputIt first, InputIt last) {
    using policy = typename base::strong_exception_policy;
    using iterator_cat = typename std::iterator_traits<InputIt>::iterator_category;
    base::template append_range<policy> (first, last, iterator_cat { });
    return *this;
  }

  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector&
  append (std::initializer_list<value_type> ilist)
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires EmplaceConstructible<const_reference>::value
  &&  MoveInsertable
#endif
  {
    return append (ilist.begin (), ilist.end ());
  }

  template <unsigned I>
  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector&
  append (const small_vector<T, I, Allocator>& other)
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires CopyInsertable
#endif
  {
    return append (other.begin (), other.end ());
  }

  template <unsigned I>
  PLUMED_GCH_CPP20_CONSTEXPR
  small_vector&
  append (small_vector<T, I, Allocator>&& other)
#ifdef PLUMED_GCH_LIB_CONCEPTS
  requires MoveInsertable
#endif
  {
    // Provide a strong exception guarantee for `other` as well.
    using move_iter_type = typename std::conditional<
      base::template relocate_with_move<value_type>::value,
      std::move_iterator<iterator>,
      iterator>::type;

    append (move_iter_type { other.begin () }, move_iter_type { other.end () });
    other.clear ();
    return *this;
  }
};

template <typename T, unsigned InlineCapacityLHS, unsigned InlineCapacityRHS, typename Allocator>
inline PLUMED_GCH_CPP20_CONSTEXPR
bool
operator== (const small_vector<T, InlineCapacityLHS, Allocator>& lhs,
            const small_vector<T, InlineCapacityRHS, Allocator>& rhs) {
  return lhs.size () == rhs.size () && std::equal (lhs.begin (), lhs.end (), rhs.begin ());
}

template <typename T, unsigned InlineCapacity, typename Allocator>
inline PLUMED_GCH_CPP20_CONSTEXPR
bool
operator== (const small_vector<T, InlineCapacity, Allocator>& lhs,
            const small_vector<T, InlineCapacity, Allocator>& rhs) {
  return lhs.size () == rhs.size () && std::equal (lhs.begin (), lhs.end (), rhs.begin ());
}

#ifdef PLUMED_GCH_LIB_THREE_WAY_COMPARISON

template <typename T, unsigned InlineCapacityLHS, unsigned InlineCapacityRHS, typename Allocator>
requires std::three_way_comparable<T>
constexpr
auto
operator<=> (const small_vector<T, InlineCapacityLHS, Allocator>& lhs,
             const small_vector<T, InlineCapacityRHS, Allocator>& rhs) {
  return std::lexicographical_compare_three_way (
    lhs.begin (), lhs.end (),
    rhs.begin (), rhs.end (),
    std::compare_three_way { });
}

template <typename T, unsigned InlineCapacity, typename Allocator>
requires std::three_way_comparable<T>
constexpr
auto
operator<=> (const small_vector<T, InlineCapacity, Allocator>& lhs,
             const small_vector<T, InlineCapacity, Allocator>& rhs) {
  return std::lexicographical_compare_three_way (
    lhs.begin (), lhs.end (),
    rhs.begin (), rhs.end (),
    std::compare_three_way { });
}

template <typename T, unsigned InlineCapacityLHS, unsigned InlineCapacityRHS, typename Allocator>
constexpr
auto
operator<=> (const small_vector<T, InlineCapacityLHS, Allocator>& lhs,
             const small_vector<T, InlineCapacityRHS, Allocator>& rhs) {
  constexpr auto comparison = [](const T& l, const T& r) {
    return (l < r) ? std::weak_ordering::less
            : (r < l) ? std::weak_ordering::greater
               : std::weak_ordering::equivalent;
  };

  return std::lexicographical_compare_three_way (
           lhs.begin (), lhs.end (),
           rhs.begin (), rhs.end (),
           comparison);
}

template <typename T, unsigned InlineCapacity, typename Allocator>
constexpr
auto
operator<=> (const small_vector<T, InlineCapacity, Allocator>& lhs,
             const small_vector<T, InlineCapacity, Allocator>& rhs) {
  constexpr auto comparison = [](const T& l, const T& r) {
    return (l < r) ? std::weak_ordering::less
            : (r < l) ? std::weak_ordering::greater
               : std::weak_ordering::equivalent;
  };

  return std::lexicographical_compare_three_way (
           lhs.begin (), lhs.end (),
           rhs.begin (), rhs.end (),
           comparison);
}

#else

template <typename T, unsigned InlineCapacityLHS, unsigned InlineCapacityRHS, typename Allocator>
inline PLUMED_GCH_CPP20_CONSTEXPR
bool
operator!= (const small_vector<T, InlineCapacityLHS, Allocator>& lhs,
            const small_vector<T, InlineCapacityRHS, Allocator>& rhs) {
  return ! (lhs == rhs);
}

template <typename T, unsigned InlineCapacity, typename Allocator>
inline PLUMED_GCH_CPP20_CONSTEXPR
bool
operator!= (const small_vector<T, InlineCapacity, Allocator>& lhs,
            const small_vector<T, InlineCapacity, Allocator>& rhs) {
  return ! (lhs == rhs);
}

template <typename T, unsigned InlineCapacityLHS, unsigned InlineCapacityRHS, typename Allocator>
inline PLUMED_GCH_CPP20_CONSTEXPR
bool
operator<  (const small_vector<T, InlineCapacityLHS, Allocator>& lhs,
const small_vector<T, InlineCapacityRHS, Allocator>& rhs) {
  return std::lexicographical_compare (lhs.begin (), lhs.end (), rhs.begin (), rhs.end ());
}

template <typename T, unsigned InlineCapacity, typename Allocator>
inline PLUMED_GCH_CPP20_CONSTEXPR
bool
operator<  (const small_vector<T, InlineCapacity, Allocator>& lhs,
const small_vector<T, InlineCapacity, Allocator>& rhs) {
  return std::lexicographical_compare (lhs.begin (), lhs.end (), rhs.begin (), rhs.end ());
}

template <typename T, unsigned InlineCapacityLHS, unsigned InlineCapacityRHS, typename Allocator>
inline PLUMED_GCH_CPP20_CONSTEXPR
bool
operator>= (const small_vector<T, InlineCapacityLHS, Allocator>& lhs,
            const small_vector<T, InlineCapacityRHS, Allocator>& rhs) {
  return ! (lhs < rhs);
}

template <typename T, unsigned InlineCapacity, typename Allocator>
inline PLUMED_GCH_CPP20_CONSTEXPR
bool
operator>= (const small_vector<T, InlineCapacity, Allocator>& lhs,
            const small_vector<T, InlineCapacity, Allocator>& rhs) {
  return ! (lhs < rhs);
}

template <typename T, unsigned InlineCapacityLHS, unsigned InlineCapacityRHS, typename Allocator>
inline PLUMED_GCH_CPP20_CONSTEXPR
bool
operator>  (const small_vector<T, InlineCapacityLHS, Allocator>& lhs,
            const small_vector<T, InlineCapacityRHS, Allocator>& rhs) {
  return rhs < lhs;
}

template <typename T, unsigned InlineCapacity, typename Allocator>
inline PLUMED_GCH_CPP20_CONSTEXPR
bool
operator>  (const small_vector<T, InlineCapacity, Allocator>& lhs,
            const small_vector<T, InlineCapacity, Allocator>& rhs) {
  return rhs < lhs;
}

template <typename T, unsigned InlineCapacityLHS, unsigned InlineCapacityRHS, typename Allocator>
inline PLUMED_GCH_CPP20_CONSTEXPR
bool
operator<= (const small_vector<T, InlineCapacityLHS, Allocator>& lhs,
const small_vector<T, InlineCapacityRHS, Allocator>& rhs) {
  return rhs >= lhs;
}

template <typename T, unsigned InlineCapacity, typename Allocator>
inline PLUMED_GCH_CPP20_CONSTEXPR
bool
operator<= (const small_vector<T, InlineCapacity, Allocator>& lhs,
const small_vector<T, InlineCapacity, Allocator>& rhs) {
  return rhs >= lhs;
}

#endif

template <typename T, unsigned InlineCapacity, typename Allocator
#ifndef PLUMED_GCH_LIB_CONCEPTS
          , typename std::enable_if<std::is_move_constructible<T>::value
                                    &&  std::is_move_assignable<T>::value
#ifdef PLUMED_GCH_LIB_IS_SWAPPABLE
                                    &&  std::is_swappable<T>::value
#endif
                                    >::type * = nullptr
#endif
          >
inline PLUMED_GCH_CPP20_CONSTEXPR
void
swap (small_vector<T, InlineCapacity, Allocator>& lhs,
      small_vector<T, InlineCapacity, Allocator>& rhs)
noexcept (noexcept (lhs.swap (rhs)))
#ifdef PLUMED_GCH_LIB_CONCEPTS
requires concepts::MoveInsertable<T, small_vector<T, InlineCapacity, Allocator>, Allocator>
&& concepts::Swappable<T>
#endif
{
  lhs.swap (rhs);
}

template <typename T, unsigned InlineCapacity, typename Allocator, typename U>
inline PLUMED_GCH_CPP20_CONSTEXPR
typename small_vector<T, InlineCapacity, Allocator>::size_type
erase (small_vector<T, InlineCapacity, Allocator>& v, const U& value) {
  const auto original_size = v.size ();
  v.erase (std::remove (v.begin (), v.end (), value), v.end ());
  return original_size - v.size ();
}

template <typename T, unsigned InlineCapacity, typename Allocator, typename Pred>
inline PLUMED_GCH_CPP20_CONSTEXPR
typename small_vector<T, InlineCapacity, Allocator>::size_type
erase_if (small_vector<T, InlineCapacity, Allocator>& v, Pred pred) {
  const auto original_size = v.size ();
  v.erase (std::remove_if (v.begin (), v.end (), pred), v.end ());
  return original_size - v.size ();
}

template <typename T, unsigned InlineCapacity, typename Allocator>
constexpr
typename small_vector<T, InlineCapacity, Allocator>::iterator
begin (small_vector<T, InlineCapacity, Allocator>& v) noexcept {
  return v.begin ();
}

template <typename T, unsigned InlineCapacity, typename Allocator>
constexpr
typename small_vector<T, InlineCapacity, Allocator>::const_iterator
begin (const small_vector<T, InlineCapacity, Allocator>& v) noexcept {
  return v.begin ();
}

template <typename T, unsigned InlineCapacity, typename Allocator>
constexpr
typename small_vector<T, InlineCapacity, Allocator>::const_iterator
cbegin (const small_vector<T, InlineCapacity, Allocator>& v) noexcept {
  return begin (v);
}

template <typename T, unsigned InlineCapacity, typename Allocator>
constexpr
typename small_vector<T, InlineCapacity, Allocator>::iterator
end (small_vector<T, InlineCapacity, Allocator>& v) noexcept {
  return v.end ();
}

template <typename T, unsigned InlineCapacity, typename Allocator>
constexpr
typename small_vector<T, InlineCapacity, Allocator>::const_iterator
end (const small_vector<T, InlineCapacity, Allocator>& v) noexcept {
  return v.end ();
}

template <typename T, unsigned InlineCapacity, typename Allocator>
constexpr
typename small_vector<T, InlineCapacity, Allocator>::const_iterator
cend (const small_vector<T, InlineCapacity, Allocator>& v) noexcept {
  return end (v);
}

template <typename T, unsigned InlineCapacity, typename Allocator>
constexpr
typename small_vector<T, InlineCapacity, Allocator>::reverse_iterator
rbegin (small_vector<T, InlineCapacity, Allocator>& v) noexcept {
  return v.rbegin ();
}

template <typename T, unsigned InlineCapacity, typename Allocator>
constexpr
typename small_vector<T, InlineCapacity, Allocator>::const_reverse_iterator
rbegin (const small_vector<T, InlineCapacity, Allocator>& v) noexcept {
  return v.rbegin ();
}

template <typename T, unsigned InlineCapacity, typename Allocator>
constexpr
typename small_vector<T, InlineCapacity, Allocator>::const_reverse_iterator
crbegin (const small_vector<T, InlineCapacity, Allocator>& v) noexcept {
  return rbegin (v);
}

template <typename T, unsigned InlineCapacity, typename Allocator>
constexpr
typename small_vector<T, InlineCapacity, Allocator>::reverse_iterator
rend (small_vector<T, InlineCapacity, Allocator>& v) noexcept {
  return v.rend ();
}

template <typename T, unsigned InlineCapacity, typename Allocator>
constexpr
typename small_vector<T, InlineCapacity, Allocator>::const_reverse_iterator
rend (const small_vector<T, InlineCapacity, Allocator>& v) noexcept {
  return v.rend ();
}

template <typename T, unsigned InlineCapacity, typename Allocator>
constexpr
typename small_vector<T, InlineCapacity, Allocator>::const_reverse_iterator
crend (const small_vector<T, InlineCapacity, Allocator>& v) noexcept {
  return rend (v);
}

template <typename T, unsigned InlineCapacity, typename Allocator>
constexpr
typename small_vector<T, InlineCapacity, Allocator>::size_type
size (const small_vector<T, InlineCapacity, Allocator>& v) noexcept {
  return v.size ();
}

template <typename T, unsigned InlineCapacity, typename Allocator>
constexpr
typename std::common_type<
  std::ptrdiff_t,
  typename std::make_signed<
    typename small_vector<T, InlineCapacity, Allocator>::size_type>::type>::type
ssize (const small_vector<T, InlineCapacity, Allocator>& v) noexcept {
  using ret_type = typename std::common_type<
    std::ptrdiff_t,
    typename std::make_signed<decltype (v.size ())>::type>::type;
  return static_cast<ret_type> (v.size ());
}

template <typename T, unsigned InlineCapacity, typename Allocator>
PLUMED_GCH_NODISCARD constexpr
bool
empty (const small_vector<T, InlineCapacity, Allocator>& v) noexcept {
  return v.empty ();
}

template <typename T, unsigned InlineCapacity, typename Allocator>
constexpr
typename small_vector<T, InlineCapacity, Allocator>::pointer
data (small_vector<T, InlineCapacity, Allocator>& v) noexcept {
  return v.data ();
}

template <typename T, unsigned InlineCapacity, typename Allocator>
constexpr
typename small_vector<T, InlineCapacity, Allocator>::const_pointer
data (const small_vector<T, InlineCapacity, Allocator>& v) noexcept {
  return v.data ();
}

#ifdef PLUMED_GCH_CTAD_SUPPORT

template <typename InputIt,
          unsigned InlineCapacity = default_buffer_size<
            std::allocator<typename std::iterator_traits<InputIt>::value_type>>::value,
          typename Allocator = std::allocator<typename std::iterator_traits<InputIt>::value_type>>
small_vector (InputIt, InputIt, Allocator = Allocator ())
-> small_vector<typename std::iterator_traits<InputIt>::value_type, InlineCapacity, Allocator>;

#endif

} // namespace gch
} // namespace PLMD

#endif // PLUMED_GCH_SMALL_VECTOR_HPP
#endif
