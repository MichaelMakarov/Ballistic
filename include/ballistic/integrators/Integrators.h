#pragma once
#include <functional>
#include <type_traits>
#include <concepts>

namespace ball
{
	template<class T> concept Arithmetic = 
		std::is_default_constructible<T>::value &&
		std::is_copy_constructible<T>::value &&
		std::is_copy_assignable<T>::value &&
		requires (T a, double b) {
			{ a * b };
			{ a / b };
			{ a *= b };
			{ a /= b };
		} && requires (T a, T b) { 
			{ a += b };
			{ a -= b };
			{ a + b };
			{ a - b };
		};
	template<class T> concept Time =
		std::is_default_constructible<T>::value &&
		std::is_copy_constructible<T>::value &&
		std::is_copy_assignable<T>::value &&
		requires (T a, double b) {
			{ a + b };
			{ a - b };
			{ a += b };
			{ a -= b };
		};

	template<class F, class Ret, class ... Args> class invoker {
	public:
		Ret operator() (Args ... args) const {
			return static_cast<const F*>(this)->operator()(std::forward<Args>(args)...);
		}
	};
	template<class Ret, class ... Args>
	class static_invoker : public invoker<static_invoker<Ret, Args...>, Ret, Args...> {
		Ret(*_func)(Args ...);
	public:
		static_invoker(Ret(* const function)(Args ...)) : _func{ function } {}
		Ret operator() (Args&& ... args) const {
			return _func(std::forward<Args>(args)...);
		}
	};
	template<class C, class Ret, class ... Args>
	class member_invoker : public invoker<member_invoker<C, Ret, Args...>, Ret, Args...> {
		Ret(C::* _func)(Args ...);
		C* _obj;
	public:
		member_invoker(Ret(C::* const func)(Args ...), C* const obj) : _func{ func }, _obj{ obj } {}
		member_invoker(Ret(C::* const func)(Args ...), const C* const obj) : _func{ func }, _obj{ obj } {}
		Ret operator() (Args&& ... args) const {
			return (_obj->*_func)(std::forward<Args>(args)...);
		}
	};
	template<class C, class Ret, class ... Args>
	class class_invoker : public invoker<class_invoker<C, Ret, Args...>, Ret, Args...> {
		C _obj;
	public:
		class_invoker(const C obj) : _obj{ obj } {}
		Ret operator() (Args&& ... args) const {
			return _obj(std::forward<Args>(args)...);
		}
	};
	
	/*struct S {
		int func(const int, const double) { return 0; }
		int operator() (const int, const double) { return 0; }
	};
	template<class Inv, class R, class T> void set(const invoker<Inv, R, const R, const T>& inv) {

	}
	void test() {
		S s;
		auto f = [](const int, const double)->int { return 1; };
		static_invoker<int, const int, const double> sinv = nullptr;
		member_invoker<S, int, const int, const double> minv(&S::func, &s);
		class_invoker<S, int, const int, const double> cinv = s;
		set(sinv);
		set(minv);
		set(cinv);
		auto fcinv = make_classfunc<decltype(f), int, const int, const double>(f);
		auto fminv = make_memberfunc(&S::func, &s);
	}*/
	template<class Ret, class ... Args> 
	static_invoker<Ret, Args...> make_staticfunc(Ret(*const func)(Args ...)) {
		return static_invoker<Ret, Args...>(func);
	}
	template<class C, class Ret, class ... Args> 
	member_invoker<C, Ret, Args...> make_memberfunc(Ret(C::*const func)(Args ...), C* const obj) {
		return member_invoker<C, Ret, Args...>(func, obj);
	}
	template<class C, class Ret, class ... Args>
	member_invoker<C, Ret, Args...> make_memberfunc(Ret(C::*const func)(Args ...), const C* const obj) {
		return member_invoker<C, Ret, Args...>(func, obj);
	}
	template<class C, class Ret, class ... Args> 
	class_invoker<C, Ret, Args...> make_classfunc(C obj) {
		return class_invoker<C, Ret, Args...>(obj);
	}

	template<class Int, Arithmetic R, Time T> class singlestep_integrator {
	public:
		template<class Inv>	void integrate(
			const R& x0, const T& t0, 
			const double step, 
			R& xk, T& tk,
			const invoker<Inv, R, const R&, const T&>& func) const {
			static_cast<const Int*>(this)->integrate(x0, t0, step, xk, tk, func);
		}
	};

	template<class Int, Arithmetic R, Time T> class multistep_integrator {
	private:
		const size_t _degree;

	public:
		explicit multistep_integrator(const size_t degree) : _degree(degree) {}

		template<class Inv> void integrate(
			R* const pR, T* const pT,
			const double step,
			R& xk, T& tk,
			const invoker<Inv, R, const R&, const T&>& func) const
		{
			static_cast<const Int*>(this)->integrate(pR, pT, step, xk, tk, func);
		}

		size_t degree() const {	return _degree;	}
	};
}