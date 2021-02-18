#pragma once
#include <vector>
#include <functional>
#include <fstream>
#include "Common.h"

namespace ball
{
	namespace graphics
	{

		class Object
		{
		private:
			std::vector<Vertex> _vertices;

		public:
			Object(
				std::ifstream& filein,
				std::function<std::ifstream&(std::ifstream& filein, std::vector<Vertex>& vertices)>& readfunc
			);
			Object(std::vector<Vertex>& vertices) noexcept : _vertices{ std::move(vertices) } {}
			Object(const Object& o) : _vertices{ o._vertices } {}
			Object(Object&& o) noexcept : _vertices{ std::move(o._vertices) } {}
			~Object() = default;

			Object& operator = (const Object& o);
			Object& operator = (Object&& o) noexcept;

			const std::vector<Vertex>& get_vertices() const { return _vertices; }
		};

	}
}