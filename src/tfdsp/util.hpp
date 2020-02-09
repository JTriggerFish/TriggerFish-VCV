#pragma once

template <typename Derived>
class enable_down_cast
{
private:
  typedef enable_down_cast Base;

public:
  Derived const *Self() const
  {
    // casting "down" the inheritance hierarchy
    return static_cast<Derived const *>(this);
  }

  Derived *Self()
  {
    return static_cast<Derived *>(this);
  }

protected:
  // disable deletion of Derived* through Base*
  // enable deletion of Base* through Derived*
  ~enable_down_cast() = default; // C++11 only, use ~enable_down_cast() {} in C++98
};
