#pragma once
#include <functional>
#include <type_traits>
#include <concepts>

namespace ball
{
	template<class T> concept Arithmetic = 
		std::is_default_constructible<T>::value &&
		requires (T a, double b) {
			{ a * b };
			{ a / b };
			{ a *= b };
			{ a /= b };
			{ b * a };
		} && requires (T a, T b) { 
			{ a += b };
			{ a -= b };
			{ a + b };
			{ a - b };
		};
	template<class T> concept Time =
		requires (T a, double b) {
			{ a + b };
			{ a - b };
			{ a += b };
			{ a -= b };
		};

	template<class F, Arithmetic R, Time T>
	class Func 
	{
	public:
		R operator () (const R& r, const T& t) const
		{
			return static_cast<const F*>(this)->operator()(r, t);
		}
	};
	template<Arithmetic R, Time T>
	class StaticFunc : public Func<StaticFunc<R, T>, R, T>
	{
	private:
		R(*_pFunc)(const R&, const T&);
	public:
		StaticFunc(R(*pfunc)(const R&, const T&)) : _pFunc{ pfunc } {}
		~StaticFunc() = default;

		R operator () (const R& r, const T& t) const
		{
			return _pFunc(r, t);
		}
	};
	template<Arithmetic R, Time T, class C>
	class MemberFunc : public Func<MemberFunc<R, T, C>, R, T>
	{
	private:
		R(C::*_pFunc)(const R&, const T&);
		C* _pObj;
	public:
		MemberFunc(R(C::*pfunc)(const R&, const T&), C* const obj) : _pFunc{ pfunc }, _pObj{ obj } {}
		MemberFunc(R(C::*pfunc)(const R&, const T&) const, const C* const obj) : _pFunc{ pfunc }, _pObj{ obj } {}
		~MemberFunc() = default;

		R operator () (const R& r, const T& t) const
		{
			return (_pObj->*_pFunc)(r, t);
		}
	};
	template<Arithmetic R, Time T, class C>
	class ClassFunc : public Func<ClassFunc<R, T, C>, R, T>
	{
	private:
		C* _pObj;
	public:
		ClassFunc(C* const obj) : _pObj{ obj } {}
		ClassFunc(const C* const obj) : _pObj{ obj } {}
		~ClassFunc() = default;

		R operator () (const R& r, const T& t) const
		{
			return (*_pObj)(r, t);
		}
	};

	template<Arithmetic R, Time T>
	StaticFunc<R, T> make_staticfunc(R(*pfunc)(const R&, const T&))
	{
		return StaticFunc<R, T>(pfunc);
	}
	template<Arithmetic R, Time T, class C>
	MemberFunc<R, T, C> make_memberfunc(R(C::* pfunc)(const R&, const T&), C* const pobj)
	{
		return MemberFunc<R, T, C>(pfunc, pobj);
	}
	template<Arithmetic R, Time T, class C>
	MemberFunc<R, T, C> make_memberfunc(R(C::* pfunc)(const R&, const T&) const, const C* const pobj)
	{
		return MemberFunc<R, T, C>(pfunc, pobj);
	}
	template<Arithmetic R, Time T, class C>
	ClassFunc<R, T, C> make_classfunc(C* const pobj)
	{
		return ClassFunc<R, T, C>(pobj);
	}
	template<Arithmetic R, Time T, class C>
	ClassFunc<R, T, C> make_classfunc(const C* const pobj)
	{
		return ClassFunc<R, T, C>(pobj);
	}


	template<Arithmetic R, Time T>
	class Integrator
	{
	public:
		Integrator() = default;
		~Integrator() = default;

		//std::function<R(const R&, const T& t)> func;
	};

	template<class Int, Arithmetic R, Time T>
	class SinglestepIntegrator : public Integrator<R, T>
	{
	public:
		SinglestepIntegrator() : Integrator<R, T>() {}
		~SinglestepIntegrator() = default;

		template<class Inv>
		void integrate(
			const R& x0,
			const T& t0,
			const double step,
			R& xk,
			T& tk,
			const Func<Inv, R, T>& func) const
		{
			static_cast<const Int*>(this)->integrate(x0, t0, step, xk, tk, func);
		}
	};

	template<class Int, Arithmetic R, Time T>
	class MultistepIntegrator : public Integrator<R, T>
	{
	private:
		const size_t _degree;

	public:
		explicit MultistepIntegrator(const size_t degree) : Integrator<R, T>(), _degree(degree) {}
		~MultistepIntegrator() = default;

		template<class Inv>
		void integrate(
			R* const pR,
			T* const pT,
			const double step,
			R& xk,
			T& tk,
			const Func<Inv, R, T>& func) const
		{
			static_cast<const Int*>(this)->integrate(pR, pT, step, xk, tk, func);
		}

		size_t degree() const {	return _degree;	}
	};
}