#ifndef HULL_H_INCLUDED
#define HULL_H_INCLUDED

#include <Eigen/Dense>
#include "Face.h"
#include "Vertex.h"
#include "Edge.h"
#include <list>
#include <memory>

class Hull
{
	private:
		std::list<std::shared_ptr<Edge>> edges;
		std::list<std::shared_ptr<Face>> faces;
		std::shared_ptr<Vertex> bottom, top, left, right;
		Eigen::Vector3d center;
		
	public:
		Hull(std::list<std::shared_ptr<Edge>> edges, std::list<std::shared_ptr<Face>> faces, std::shared_ptr<Vertex> bottom, std::shared_ptr<Vertex> top, std::shared_ptr<Vertex> left, std::shared_ptr<Vertex> right, Eigen::Vector3d center);
		std::list<std::shared_ptr<Face>>& getFaces();
		std::list<std::shared_ptr<Edge>>& getEdges();
		std::shared_ptr<Vertex> getBottom() const;
		std::shared_ptr<Vertex> getTop() const;
		std::shared_ptr<Vertex> getLeft() const;
		std::shared_ptr<Vertex> getRight() const;
		Eigen::Vector3d getCenter() const;
		void setBottom(std::shared_ptr<Vertex> v);
		void setTop(std::shared_ptr<Vertex> v);
		void setLeft(std::shared_ptr<Vertex> v);
		void setRight(std::shared_ptr<Vertex> v);
		void setCenter(Eigen::Vector3d center);
		void addFace(std::shared_ptr<Face> face);
		void addEdge(std::shared_ptr<Edge> edge);
		void deleteFace(std::list<std::shared_ptr<Face>>::iterator face);
		void deleteFace(std::list<std::shared_ptr<Edge>>::iterator edge);
		void deleteVisibles();
		void updateExtremities(std::shared_ptr<Vertex> bottom, std::shared_ptr<Vertex> top, std::shared_ptr<Vertex> left, std::shared_ptr<Vertex> right);
};

#endif