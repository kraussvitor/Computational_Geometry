#include "Edge.h"
#include <iostream>

Edge::Edge() 
{
	used = false;
}

Edge::Edge(const Edge& edge)
{
	std::cout << "-- COPYING EDGE --" << std::endl;
}

void Edge::setEndPts(std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2)
{
	endpts[0] = v1;
	endpts[1] = v2;
}

void Edge::setAdjFaces(std::shared_ptr<Face> f1, std::shared_ptr<Face> f2)
{
	adjfaces[0] = f1;
	adjfaces[1] = f2;
}

void Edge::setAdjFace(std::shared_ptr<Face> f, unsigned int i)
{
	adjfaces[i].reset(f.get()); // adjfaces[i] = f;
}

std::shared_ptr<Vertex> Edge::getEndPts(unsigned int i)
{
	return endpts[i];
}

std::shared_ptr<Face> Edge::getAdjFace(unsigned int i)
{
	return adjfaces[i];
}

bool Edge::getUsed()
{
	return used;
}

void Edge::setUsed(bool used)
{
	this->used = used;
}