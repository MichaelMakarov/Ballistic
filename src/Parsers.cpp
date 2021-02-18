#include "Parsers.h"

namespace ball
{
	namespace graphics
	{
		// Object represents a single triangle from blender
		struct Face
		{
			size_t Vertind[3];	// indices of the vertices
			size_t Textind[3];	// indices of the texture
			size_t Normind[3];	// indices of the normals
		};

		
		std::ifstream& load_blender_model(std::ifstream& filein, std::vector<Vertex>& vertices)
		{
			using namespace general::math;
			if (!filein.is_open()) {
				throw std::invalid_argument("filein is not opened!");
			}
			size_t vertcount{ 0 }, textcount{ 0 }, normcount{ 0 }, facecount{ 0 };
			char input;
			filein.get(input);
			while (!filein.eof()) {
				if (input == 'v') {
					filein.get(input);
					if (input == ' ') vertcount++;
					else if (input == 't') textcount++;
					else if (input == 'n') normcount++;
				}
				else if (input == 'f') {
					filein.get(input);
					if (input == ' ') facecount++;
				}
				while (input != '\n') filein.get(input);
				try {
					filein.get(input);
				} catch (std::ifstream::failure&) {}
			}
			filein.clear();
			filein.seekg(0, filein.beg);

			auto points{ std::vector<Vec3>(vertcount) };
			auto textures{ std::vector<Vec2>(textcount) };
			auto normals{ std::vector<Vec3>(normcount) };
			auto faces{ std::vector<Face>(facecount) };
			vertcount = textcount = normcount = facecount = 0;

			filein.get(input);
			while (!filein.eof()) {
				if (input == 'v') {
					filein.get(input);
					if (input == ' ') {
						filein >> points[vertcount];
						points[vertcount].Z = -points[vertcount].Z;
						vertcount++;
					}
					else if (input == 't') {
						filein >> textures[textcount++];
					}
					else if (input == 'n') {
						filein >> normals[normcount];
						normals[normcount].Z = -normals[normcount].Z;
						normcount++;
					}
				}
				else if (input == 'f') {
					//filein.get(input);
					if (input == ' ') {
						filein >>
							input >> faces[facecount].Vertind[2] >>
							input >> faces[facecount].Textind[2] >>
							input >> faces[facecount].Normind[2] >>
							input >> faces[facecount].Vertind[1] >>
							input >> faces[facecount].Textind[1] >>
							input >> faces[facecount].Normind[1] >>
							input >> faces[facecount].Vertind[0] >>
							input >> faces[facecount].Textind[0] >>
							input >> faces[facecount].Normind[0];
						facecount++;
					}
				}
				while (input != '\n') filein.get(input);
				try {
					filein.get(input);
				}
				catch (std::ifstream::failure&) {}
			}

			vertices.resize(faces.size() * 3);
			vertcount = textcount = normcount = facecount = 0;
			size_t a{ 0 };
			for (size_t i = 0; i < faces.size(); ++i) {
				vertcount = faces[a].Vertind[0] - 1;
				textcount = faces[a].Textind[0] - 1;
				normcount = faces[a].Normind[0] - 1;
				vertices[i].Pos = points[vertcount];
				vertices[i].Nor = points[normcount];
				++i;
				vertcount = faces[a].Vertind[1] - 1;
				textcount = faces[a].Textind[1] - 1;
				normcount = faces[a].Normind[1] - 1;
				vertices[i].Pos = points[vertcount];
				vertices[i].Nor = points[normcount];
				++i;
				vertcount = faces[a].Vertind[2] - 1;
				textcount = faces[a].Textind[2] - 1;
				normcount = faces[a].Normind[2] - 1;
				vertices[i].Pos = points[vertcount];
				vertices[i].Nor = points[normcount];
			}
			return filein;
		}
	}
}