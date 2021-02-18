#include "Object.h"

namespace ball
{
	namespace graphics
	{
		Object::Object(
			std::ifstream& filein, 
			std::function<std::ifstream& (std::ifstream& filein, std::vector<Vertex>& vertices)>& readfunc)
		{
			try {
				readfunc(filein, _vertices);
			}
			catch (std::exception& ) { return; }
		}
		Object& Object::operator=(const Object& o)
		{
			_vertices = o._vertices;
			return *this;
		}
		Object& Object::operator=(Object&& o) noexcept
		{
			_vertices = std::move(o._vertices);
			return *this;
		}
	}
}