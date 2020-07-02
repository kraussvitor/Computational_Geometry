#include "Vertex.h"

Vertex::Vertex()
{
	this->location = Eigen::Vector3d(0.0, 0.0, 0.0);
}

Vertex::Vertex(Eigen::Vector3d location)
{
	this->location = location;
}

Eigen::Vector3d Vertex::getLocation() const
{
	return location;
}