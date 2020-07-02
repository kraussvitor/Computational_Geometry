#include "Hull.h"
#include <iostream>

Hull::Hull(std::list<std::shared_ptr<Edge>> edges, std::list<std::shared_ptr<Face>> faces, std::shared_ptr<Vertex> bottom, std::shared_ptr<Vertex> top, std::shared_ptr<Vertex> left, std::shared_ptr<Vertex> right, Eigen::Vector3d center)
{
	this->edges = edges;
	this->faces = faces;
	this->bottom = bottom;
	this->top = top;
	this->left = left;
	this->right = right;
	this->center = center;
}

std::list<std::shared_ptr<Face>>& Hull::getFaces()
{
	return faces;
}

std::list<std::shared_ptr<Edge>>& Hull::getEdges()
{
	return edges;
}

std::shared_ptr<Vertex> Hull::getBottom() const
{
	return bottom;
}

std::shared_ptr<Vertex> Hull::getTop() const
{
	return top;
}

std::shared_ptr<Vertex> Hull::getLeft() const
{
	return left;
}

std::shared_ptr<Vertex> Hull::getRight() const
{
	return right;
}

Eigen::Vector3d Hull::getCenter() const
{
	return center;
}

void Hull::setBottom(std::shared_ptr<Vertex> v)
{
	bottom = v;
}

void Hull::setTop(std::shared_ptr<Vertex> v)
{
	top = v;
}

void Hull::setLeft(std::shared_ptr<Vertex> v)
{
	left = v;
}

void Hull::setRight(std::shared_ptr<Vertex> v)
{
	right = v;
}

void Hull::setCenter(Eigen::Vector3d center)
{
	this->center = center;
}

void Hull::addFace(std::shared_ptr<Face> face)
{
	faces.push_back(face);
}

void Hull::addEdge(std::shared_ptr<Edge> edge)
{
	edges.push_back(edge);
}

void Hull::deleteFace(std::list<std::shared_ptr<Face>>::iterator face)
{
	faces.erase(face);
}
void Hull::deleteFace(std::list<std::shared_ptr<Edge>>::iterator edge)
{
	edges.erase(edge);
}

void Hull::deleteVisibles()
{
	unsigned int count;
	std::list<std::shared_ptr<Edge>>::iterator e = edges.begin();
	while(e != edges.end())
	{
		count = 0;
		if ((*e)->getAdjFace(0))
		{
			if ((*e)->getAdjFace(0)->getIsVisible())
			{
				count++;
			}
			else
			{
				count--;
			}
		}
		else
		{
			count--;
		}

		if ((*e)->getAdjFace(1))
		{
			if ((*e)->getAdjFace(1)->getIsVisible())
			{
				count++;
			}
			else
			{
				count--;
			}
		}
		else
		{
			count--;
		}

		if (count == 2)
		{
			e = edges.erase(e);
		}
		else
		{
			e++;
		}
	}
	std::list<std::shared_ptr<Face>>::iterator f = faces.begin();
	while (f != faces.end())
	{
		if ((*f)->getIsVisible())
		{
			f = faces.erase(f);
		}
		else
		{
			f++;
		}
	}
}

void Hull::updateExtremities(std::shared_ptr<Vertex> bottom, std::shared_ptr<Vertex> top, std::shared_ptr<Vertex> left, std::shared_ptr<Vertex> right)
{
	if (bottom->getLocation()(2) < this->bottom->getLocation()(2))
	{
		this->bottom = bottom;
	}
	if (top->getLocation()(2) > this->top->getLocation()(2))
	{
		this->top = top;
	}
	if (left->getLocation()(1) < this->left->getLocation()(1))
	{
		this->left = left;
	}
	if (right->getLocation()(1) > this->right->getLocation()(1))
	{
		this->right = right;
	}
}
