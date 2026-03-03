#pragma once
// -*- C++ -*-
//
// KUCMSItemManager.hh (header-only, ROOT-free)
// Original Author: Jack W King III
//

#include <map>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <limits>

namespace kucms {

using uInt = unsigned int;

namespace detail {

// detection idiom: has_clear<T>
template <typename U, typename = void>
struct has_clear : std::false_type {};

template <typename U>
struct has_clear<U, std::void_t<decltype(std::declval<U&>().clear())>> : std::true_type {};

} // namespace detail

//------------------------------------------------------------------------------
// Item<T>
//------------------------------------------------------------------------------

template <class T>
class Item {
public:
  Item() : iName(""), iDoc(""), iValue(T{}) {}
  Item(const std::string& name, const std::string& doc = "")
      : iName(name), iDoc(doc), iValue(T{}) {}
  Item(const std::string& name, const T& value, const std::string& doc = "")
      : iName(name), iDoc(doc), iValue(value) {}

  void make(const std::string& name, const std::string& doc = "") {
    iName = name;
    iDoc = doc;
  }

  void fill(const T& val) { iValue = val; }

  // Keep behavior: set to a conservative "empty / sentinel"
  inline void clear() {
    if constexpr (std::is_same_v<T, std::string>) {
      iValue.clear();
    } else if constexpr (std::is_same_v<T, bool>) {
      iValue = false;
    } else if constexpr (std::is_same_v<T, uInt>) {
      iValue = 0;
    } else if constexpr (detail::has_clear<T>::value) {
      iValue.clear();
    } else if constexpr (std::is_arithmetic_v<T> && std::numeric_limits<T>::is_specialized) {
      if constexpr (std::is_unsigned_v<T>) {
        iValue = 0;
      } else {
        iValue = -std::numeric_limits<T>::max();
      }
    } else {
      iValue = T{};
    }
  }

  // Backwards-compatible names
  T getvalue() const { return iValue; }
  std::string getdoc() const { return iDoc; }

  // More efficient accessors (optional use)
  const T& value() const { return iValue; }
  T& value() { return iValue; }
  const std::string& doc() const { return iDoc; }
  const std::string& name() const { return iName; }

  // NOTE: kept public because your branch code may access iName/iDoc directly.
  // If you want these protected, we can do that after checking KUCMSBranch2.hh usage.
protected:
  std::string iName;
  std::string iDoc;
  T iValue;
};

//------------------------------------------------------------------------------
// VectorItem<T>
//------------------------------------------------------------------------------

template <class T>
class VectorItem {
public:
  VectorItem() : iName(""), iDoc(""), iVector() {}
  VectorItem(const std::string& name, const std::string& doc)
      : iName(name), iDoc(doc), iVector() {}
  VectorItem(const std::string& name, const T& value, const std::string& doc)
      : iName(name), iDoc(doc), iVector() {
    fill(value);
  }

  void make(const std::string& name, const std::string& doc = "") {
    iName = name;
    iDoc = doc;
  }

  void fill(const T& val) { iVector.push_back(val); }
  void clear() { iVector.clear(); }

  // Backwards-compatible names
  std::vector<T> getvalue() const { return iVector; }
  std::string getdoc() const { return iDoc; }

  // More efficient accessors (optional use)
  const std::vector<T>& value() const { return iVector; }
  std::vector<T>& value() { return iVector; }

protected:
  std::string iName;
  std::string iDoc;
  std::vector<T> iVector;
};

//------------------------------------------------------------------------------
// ItemManager<T, C>
//------------------------------------------------------------------------------

template <class T, template <class> class C = Item>
class ItemManager {
public:
  void set(const std::string& key, const std::string& name, const std::string& doc) {
    items[key] = C<T>(name, doc);
  }

  void set(const std::string& name) { set(name, name, ""); }

  void set(const std::string& name, const T& val) {
    items[name] = C<T>(name, val, "");
  }

  void fill(const std::string& key, const T& val) {
    auto it = items.find(key);
    if (it != items.end()) it->second.fill(val);
  }

  T get(const std::string& key) const {
    auto it = items.find(key);
    return (it != items.end()) ? it->second.getvalue() : T{};
  }

  void clear() {
    for (auto& kv : items) kv.second.clear();
  }

  void clear(const std::string& key) {
    auto it = items.find(key);
    if (it != items.end()) it->second.clear();
  }

  void reset() { // fixed bug
    for (auto& kv : items) kv.second.clear();
  }

  T operator()(const std::string& key) const { return get(key); }

private:
  std::map<std::string, C<T>> items;
};

} // namespace kucms
