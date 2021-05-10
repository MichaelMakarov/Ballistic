#include <iostream>
#include <fstream>
#include <concepts>
#include <type_traits>
#include <vector>
#include <array>
#include <string>
#include "general/Times.h"

//template<class T> concept Arithmetic =
//std::is_default_constructible<T>::value &&
//requires (T a, double b) {
//	{ a* b };
//	{ a / b };
//	{ a *= b };
//	{ a /= b };
//	{ b* a };
//} && requires (T a, T b) {
//	{ a += b };
//	{ a -= b };
//	{ a + b };
//	{ a - b };
//};
//
//template<Arithmetic A> A mult_with(const A& a, const A& b, const double n)
//{
//	return a * b * n;
//}

template<class ... Args>
class Base
{
public:
	void(*func)(const bool flag, const Args& ... args);

};

template<class ... Args>
class Der : public Base<Args>
{
public:
	void invoke_func(const Args& ... args)
	{
		func(args);
	}
};

template<typename T>
std::vector<std::vector<T>> create_2d_array(const size_t m, const size_t n)
{
	auto result{ std::vector<std::vector<T>>(m) };
	for (auto& arr : result) arr.resize(n);
	return std::move(result);
}

bool check(std::ostream& out)
{
	out << "test: \n";
	return static_cast<bool>(out);
}

template<size_t n = 4>
void test()
{
	double arr[n];
	for (size_t i = 0; i < n; ++i) {
		std::cout << arr[i] << " ";
	}
	std::cout << std::endl;
}

int main()
{
	test<5>();
	auto jd{ general::time::JD(general::time::DateTime(2000, 1, 1, 12, 0, 0)) };
	std::cout << jd << std::endl;
	return 0;
}