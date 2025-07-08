#ifndef SZ_DEBUG_HPP
#define SZ_DEBUG_HPP

#include <iostream>

#define DEBUG_PTR(var) std::cout << "[Debug] [" << __LINE__ << "] " << #var << " = 0x" \
    << std::hex << reinterpret_cast<uintptr_t>(var) << std::dec << std::endl;

#define DEBUG_VAR(var) std::cout << "[Debug] " << #var " = " << (var) << std::endl;

template <typename T, typename = void>
struct is_iterable : std::false_type {};

template <typename T>
struct is_iterable<T, std::void_t<
                          decltype(std::begin(std::declval<T>())),
                          decltype(std::end(std::declval<T>()))
                          >> : std::true_type {};

template <typename Container,
          typename = std::enable_if_t<
              is_iterable<Container>::value &&
              !std::is_same<Container, std::string>::value &&
              !std::is_same<Container, const char*>::value
              >>
std::ostream& operator<<(std::ostream& os, const Container& container) {
    os << "[";
    bool first = true;
    for (const auto& element : container) {
        if (!first) {
            os << ", ";
        }
        os << element;
        first = false;
    }
    os << "]";
    return os;
}

#endif