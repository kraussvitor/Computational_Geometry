#ifndef EDGE_H_INCLUDED
#define EDGE_H_INCLUDED

#include <memory>

class Face;
class Vertex;

class Edge
{
	private: 
		std::shared_ptr<Vertex> endpts[2];
		std::shared_ptr<Face> adjfaces[2];
		bool used;
			   
	public:
		Edge();
		Edge(const Edge& edge);
		void setEndPts(std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2);
		void setAdjFaces(std::shared_ptr<Face> f1, std::shared_ptr<Face> f2);
		void setAdjFace(std::shared_ptr<Face> f, unsigned int i);
		std::shared_ptr<Vertex> getEndPts(unsigned int i);
		std::shared_ptr<Face> getAdjFace(unsigned int i);
		bool getUsed();
		void setUsed(bool used);
};

#endif
