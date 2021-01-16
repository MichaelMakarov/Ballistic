#pragma once
#include <ostream>

namespace ball
{
	namespace time
	{
		// Class implements the functionality of the date.
		// It contains year, month and day.
		class Date
		{
		private:
			size_t _year;
			unsigned short _month;
			unsigned short _day;

			bool Check(
				const size_t year, 
				const unsigned short month, 
				const unsigned short day);

		public:
			Date() : _year{ 0 }, _month{ 0 }, _day{ 0 } {}
			Date(
				const size_t year,
				const unsigned short month,
				const unsigned short day);
			~Date() noexcept {}

			size_t get_year() const { return _year; }
			unsigned short get_month() const { return _month; }
			unsigned short get_day() const { return _day; }

			friend bool operator > (const Date& f, const Date& s);
			friend bool operator < (const Date& f, const Date& s);
			friend bool operator == (const Date& f, const Date& s);

			friend std::ostream& operator << (std::ostream& o, const Date& d);
		};

		// Class implements time of the day;
		// Contains hour, minute, second and millisecond.
		class Time
		{
		private:
			unsigned short _hour, _minute, _second, _millisec;

		public:
			Time() : _hour{ 0 }, _minute{ 0 }, _second{ 0 }, _millisec{ 0 } {}
			Time(
				const unsigned short hour,
				const unsigned short minute,
				const unsigned short second,
				const unsigned short millisecond = 0);
			~Time() {}

			unsigned short get_hour() const { return _hour; }
			unsigned short get_minute() const { return _minute; }
			unsigned short get_second() const { return _second; }
			unsigned short get_millisecond() const { return _millisec; }

			friend bool operator > (const Time& f, const Time& s);
			friend bool operator < (const Time& f, const Time& s);
			friend bool operator == (const Time& f, const Time& s);

			friend std::ostream& operator << (std::ostream& o, const Time& d);
		};

		// Class implements both date and time od the day.
		class DateTime
		{
		private:
			Date _date;
			Time _time;

		public:
			DateTime() : _date{}, _time{} {}
			DateTime(const Date& date, const Time& time) : _date{ date }, _time{ time } {}
			DateTime(
				const size_t year,
				const unsigned short month,
				const unsigned short day,
				const unsigned short hour,
				const unsigned short minute,
				const unsigned short second,
				const unsigned short millisecond = 0) :
				_date{ Date(year, month, day) }, 
				_time{ Time(hour, minute, second, millisecond) }
			{}

			const Date& get_date() const { return _date; }
			const Time& get_time() const { return _time; }

			friend bool operator > (const DateTime& f, const DateTime& s);
			friend bool operator < (const DateTime& f, const DateTime& s);
			friend bool operator == (const DateTime& f, const DateTime& s);

			friend std::ostream& operator << (std::ostream& o, const DateTime& d);
		};

		class JD
		{
		private:
			size_t _day;
			double _time;

		public:
			JD() : _day{ 0 }, _time{ 0 } {}
			explicit JD(const DateTime& datetime);
			JD(const size_t day, const double time) : _day{ day }, _time{ time } {}
			JD(const double jd);
			~JD() {}

			JD& operator = (const double jd);
			
			size_t JDN() const { return _day; }
			size_t T() const { return _time; }
			DateTime to_datetime() const;
			double to_double() const;

			void add_days(const int n);
			void add_hours(const int n);
			void add_minutes(const int n);
			void add_seconds(const int n);

			friend JD operator + (const JD& jd, const double dt);
			friend JD operator - (const JD& jd, const double dt);
			friend double operator - (const JD& f, const JD& s);

			friend bool operator < (const JD& f, const JD& s);
			friend bool operator > (const JD& f, const JD& s);
			friend bool operator == (const JD& f, const JD& s);
			friend bool operator <= (const JD& f, const JD& s);
			friend bool operator >= (const JD& f, const JD& s);

			friend std::ostream& operator << (std::ostream& o, const JD& jd);
		};
	}
}