#ifndef FACE_H_INCLUDED
#define FACE_H_INCLUDED

#include <memory>
#include <Eigen/Dense>
#include "Vertex.h"

class Edge;

class Face
{
	private:
		std::shared_ptr<Edge> edges[3];
		std::shared_ptr<Vertex> vertices[3];
		Eigen::Vector3d normal;
		Eigen::Vector3d center;
		bool isVisible;

	public:
		Face();
		Face(const Face& face);
		void setEdges(std::shared_ptr<Edge> e1, std::shared_ptr<Edge> e2, std::shared_ptr<Edge> e3);
		void setEdge(std::shared_ptr<Edge> e, unsigned int i);
		void setVertices(std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2, std::shared_ptr<Vertex> v3);
		void setNormal(Eigen::Vector3d normal);
		void setIsVisible(bool isVisible);
		bool checkVisibility(std::shared_ptr<Vertex> v);
		bool getIsVisible() const;
		Eigen::Vector3d getNormal() const;
		Eigen::Vector3d getCenter() const;
		std::shared_ptr<Edge> getEdge(unsigned int i) const;
		std::shared_ptr<Vertex> getVertex(unsigned int i) const;
};


#endif 