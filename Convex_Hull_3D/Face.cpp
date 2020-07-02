#include "Face.h"
#include "Edge.h"
#include <Eigen/Dense>
#include <iostream>

Face::Face()
{
	isVisible = false;
	edges[0] = std::make_shared<Edge>();
	edges[1] = std::make_shared<Edge>();
	edges[2] = std::make_shared<Edge>();
	vertices[0] = std::make_shared<Vertex>();
	vertices[1] = std::make_shared<Vertex>();
	vertices[2] = std::make_shared<Vertex>();
}

Face::Face(const Face& face)
{
	std::cout << "-- COPYING A FACE --" << std::endl;
}


void Face::setEdges(std::shared_ptr<Edge> e1, std::shared_ptr<Edge> e2, std::shared_ptr<Edge> e3)
{
	edges[0] = e1;
	edges[1] = e2;
	edges[2] = e3;
}

void Face::setEdge(std::shared_ptr<Edge> e, unsigned int i)
{
	edges[i] = e;
}

void Face::setVertices(std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2, std::shared_ptr<Vertex> v3)
{
	vertices[0] = v1;
	vertices[1] = v2;
	vertices[2] = v3;
	center = (v1->getLocation() + v2->getLocation() + v3->getLocation()) / 3.0;
}

void Face::setNormal(Eigen::Vector3d normal)
{
	this->normal = normal;
}

void Face::setIsVisible(bool isVisible)
{
	this->isVisible = isVisible;
}

bool Face::checkVisibility(std::shared_ptr<Vertex> v)
{
	Eigen::Vector3d diff = this->center - v->getLocation();
	if (diff.dot(this->normal) > 1e-9)
	{
		return false;
	}
	else
	{
		return true;
	}
}

bool Face::getIsVisible() const
{
	return isVisible;
}

Eigen::Vector3d Face::getNormal() const
{
	return normal;
}

std::shared_ptr<Edge> Face::getEdge(unsigned int i) const
{
	return edges[i];
}

std::shared_ptr<Vertex> Face::getVertex(unsigned int i) const
{
	return vertices[i];
}

Eigen::Vector3d Face::getCenter() const
{
	return center;
}