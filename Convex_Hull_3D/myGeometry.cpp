#include "myGeometry.h"
#include <iostream>
#include <random>
#include <cmath>

Eigen::Vector3d myGeometry::outwardNormal(Eigen::Vector3d p1, Eigen::Vector3d p2, Eigen::Vector3d p3, Eigen::Vector3d centroid)
{
	Eigen::Vector3d diff1 = p2 - p1;
	Eigen::Vector3d diff2 = p3 - p1;
	Eigen::Vector3d normal = diff1.cross(diff2);
	Eigen::Vector3d tCenter = 0.3333 * (p1 + p2 + p3);
	if (normal.dot(centroid - tCenter) < -1e-10)
	{
		return normal / normal.norm();
	}
	else
	{
		return (-1.0 / normal.norm()) * normal;
	}
}

std::shared_ptr<Vertex> myGeometry::findBottom(std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2, std::shared_ptr<Vertex> v3)
{
	if (v1->getLocation()(2) < v2->getLocation()(2) && v1->getLocation()(2) < v3->getLocation()(2))
	{
		return v1;
	}
	else if (v2->getLocation()(2) < v1->getLocation()(2) && v2->getLocation()(2) < v3->getLocation()(2))
	{
		return v2;
	}
	else
	{
		return v3;
	}
}

std::shared_ptr<Vertex> myGeometry::findBottom(std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2, std::shared_ptr<Vertex> v3, std::shared_ptr<Vertex> v4)
{
	if (v1->getLocation()(2) < v2->getLocation()(2) && v1->getLocation()(2) < v3->getLocation()(2) && v1->getLocation()(2) < v4->getLocation()(2))
	{
		return v1;
	}
	else if (v2->getLocation()(2) < v1->getLocation()(2) && v2->getLocation()(2) < v3->getLocation()(2) && v2->getLocation()(2) < v4->getLocation()(2))
	{
		return v2;
	}
	else if (v3->getLocation()(2) < v1->getLocation()(2) && v3->getLocation()(2) < v2->getLocation()(2) && v3->getLocation()(2) < v4->getLocation()(2))
	{
		return v3;
	}
	else
	{
		return v4;
	}
}

std::shared_ptr<Vertex> myGeometry::findTop(std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2, std::shared_ptr<Vertex> v3)
{
	if (v1->getLocation()(2) > v2->getLocation()(2) && v1->getLocation()(2) > v3->getLocation()(2))
	{
		return v1;
	}
	else if (v2->getLocation()(2) > v1->getLocation()(2) && v2->getLocation()(2) > v3->getLocation()(2))
	{
		return v2;
	}
	else
	{
		return v3;
	}
}

std::shared_ptr<Vertex> myGeometry::findTop(std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2, std::shared_ptr<Vertex> v3, std::shared_ptr<Vertex> v4)
{
	if (v1->getLocation()(2) > v2->getLocation()(2) && v1->getLocation()(2) > v3->getLocation()(2) && v1->getLocation()(2) > v4->getLocation()(2))
	{
		return v1;
	}
	else if (v2->getLocation()(2) > v1->getLocation()(2) && v2->getLocation()(2) > v3->getLocation()(2) && v2->getLocation()(2) > v4->getLocation()(2))
	{
		return v2;
	}
	else if (v3->getLocation()(2) > v1->getLocation()(2) && v3->getLocation()(2) > v2->getLocation()(2) && v3->getLocation()(2) > v4->getLocation()(2))
	{
		return v3;
	}
	else
	{
		return v4;
	}
}

std::shared_ptr<Vertex> myGeometry::findLeft(std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2, std::shared_ptr<Vertex> v3)
{
	if (v1->getLocation()(1) < v2->getLocation()(1) && v1->getLocation()(1) < v3->getLocation()(1))
	{
		return v1;
	}
	else if (v2->getLocation()(1) < v1->getLocation()(1) && v2->getLocation()(1) < v3->getLocation()(1))
	{
		return v2;
	}
	else
	{
		return v3;
	}
}

std::shared_ptr<Vertex> myGeometry::findLeft(std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2, std::shared_ptr<Vertex> v3, std::shared_ptr<Vertex> v4)
{
	if (v1->getLocation()(1) < v2->getLocation()(1) && v1->getLocation()(1) < v3->getLocation()(1) && v1->getLocation()(1) < v4->getLocation()(1))
	{
		return v1;
	}
	else if (v2->getLocation()(1) < v1->getLocation()(1) && v2->getLocation()(1) < v3->getLocation()(1) && v2->getLocation()(1) < v4->getLocation()(1))
	{
		return v2;
	}
	else if (v3->getLocation()(1) < v1->getLocation()(1) && v3->getLocation()(1) < v2->getLocation()(1) && v3->getLocation()(1) < v4->getLocation()(1))
	{
		return v3;
	}
	else
	{
		return v4;
	}
}

std::shared_ptr<Vertex> myGeometry::findRight(std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2, std::shared_ptr<Vertex> v3)
{
	if (v1->getLocation()(1) > v2->getLocation()(1) && v1->getLocation()(1) > v3->getLocation()(1))
	{
		return v1;
	}
	else if (v2->getLocation()(1) > v1->getLocation()(1) && v2->getLocation()(1) > v3->getLocation()(1))
	{
		return v2;
	}
	else
	{
		return v3;
	}
}

std::shared_ptr<Vertex> myGeometry::findRight(std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2, std::shared_ptr<Vertex> v3, std::shared_ptr<Vertex> v4)
{
	if (v1->getLocation()(1) > v2->getLocation()(1) && v1->getLocation()(1) > v3->getLocation()(1) && v1->getLocation()(1) > v4->getLocation()(1))
	{
		return v1;
	}
	else if (v2->getLocation()(1) > v1->getLocation()(1) && v2->getLocation()(1) > v3->getLocation()(1) && v2->getLocation()(1) > v4->getLocation()(1))
	{
		return v2;
	}
	else if (v3->getLocation()(1) > v1->getLocation()(1) && v3->getLocation()(1) > v2->getLocation()(1) && v3->getLocation()(1) > v4->getLocation()(1))
	{
		return v3;
	}
	else
	{
		return v4;
	}
}

Hull myGeometry::initializeTriangle(const std::vector<Eigen::Vector3d> &points)
{
	std::shared_ptr<Vertex> v1 = std::make_shared<Vertex>(points[0]);
	std::shared_ptr<Vertex> v2 = std::make_shared<Vertex>(points[1]);
	std::shared_ptr<Vertex> v3 = std::make_shared<Vertex>(points[2]);

	std::shared_ptr<Edge> e1 = std::make_shared<Edge>();
	std::shared_ptr<Edge> e2 = std::make_shared<Edge>();
	std::shared_ptr<Edge> e3 = std::make_shared<Edge>();
	e1->setEndPts(v1, v2);
	e2->setEndPts(v2, v3);
	e3->setEndPts(v3, v1);

	std::shared_ptr<Face> f1 = std::make_shared<Face>();
	std::shared_ptr<Face> f2 = std::make_shared<Face>();

	f1->setVertices(v1, v2, v3);
	f2->setVertices(v1, v2, v3);
	f1->setEdges(e1, e2, e3);
	f2->setEdges(e1, e2, e3);
	Eigen::Vector3d diff1 = points[1] - points[0];
	Eigen::Vector3d diff2 = points[2] - points[0];
	Eigen::Vector3d normal = diff1.cross(diff2);
	normal = normal / normal.norm();
	f1->setNormal(normal);
	f2->setNormal(-1.0 * normal);

	e1->setAdjFaces(f1, f2);
	e2->setAdjFaces(f1, f2);
	e3->setAdjFaces(f1, f2);

	std::list<std::shared_ptr<Edge>> edges;
	std::list<std::shared_ptr<Face>> faces;

	edges.push_back(e1);
	edges.push_back(e2);
	edges.push_back(e3);

	faces.push_back(f1);
	faces.push_back(f2);

	Hull hull(
		edges,
		faces,
		myGeometry::findBottom(v1, v2, v3),
		myGeometry::findTop(v1, v2, v3),
		myGeometry::findLeft(v1, v2, v3),
		myGeometry::findRight(v1, v2, v3),
		(points[0] + points[1] + points[2]) / 3.0);

	return hull;
}

Hull myGeometry::initializeTetrahedron(const std::vector<Eigen::Vector3d> &points)
{
	std::shared_ptr<Vertex> v1 = std::make_shared<Vertex>(points[0]);
	std::shared_ptr<Vertex> v2 = std::make_shared<Vertex>(points[1]);
	std::shared_ptr<Vertex> v3 = std::make_shared<Vertex>(points[2]);
	std::shared_ptr<Vertex> v4 = std::make_shared<Vertex>(points[3]);

	std::shared_ptr<Edge> e1 = std::make_shared<Edge>();
	std::shared_ptr<Edge> e2 = std::make_shared<Edge>();
	std::shared_ptr<Edge> e3 = std::make_shared<Edge>();
	std::shared_ptr<Edge> e4 = std::make_shared<Edge>();
	std::shared_ptr<Edge> e5 = std::make_shared<Edge>();
	std::shared_ptr<Edge> e6 = std::make_shared<Edge>();

	e1->setEndPts(v1, v2);
	e2->setEndPts(v2, v3);
	e3->setEndPts(v3, v1);
	e4->setEndPts(v1, v4);
	e5->setEndPts(v2, v4);
	e6->setEndPts(v3, v4);

	std::shared_ptr<Face> f1 = std::make_shared<Face>();
	std::shared_ptr<Face> f2 = std::make_shared<Face>();
	std::shared_ptr<Face> f3 = std::make_shared<Face>();
	std::shared_ptr<Face> f4 = std::make_shared<Face>();

	f1->setVertices(v1, v2, v3);
	f2->setVertices(v1, v2, v4);
	f3->setVertices(v2, v3, v4);
	f4->setVertices(v1, v3, v4);
	f1->setEdges(e1, e2, e3);
	f2->setEdges(e1, e4, e5);
	f3->setEdges(e2, e5, e6);
	f4->setEdges(e3, e4, e6);
	
	Eigen::Vector3d centroid = 0.25 * (points[0] + points[1] + points[2] + points[3]);

	f1->setNormal(myGeometry::outwardNormal(points[0], points[1], points[2], centroid));
	f2->setNormal(myGeometry::outwardNormal(points[0], points[1], points[3], centroid));
	f3->setNormal(myGeometry::outwardNormal(points[1], points[2], points[3], centroid));
	f4->setNormal(myGeometry::outwardNormal(points[0], points[2], points[3], centroid));

	e1->setAdjFaces(f1, f2);
	e2->setAdjFaces(f1, f3);
	e3->setAdjFaces(f1, f4);
	e4->setAdjFaces(f2, f4);
	e5->setAdjFaces(f2, f3);
	e6->setAdjFaces(f3, f4);

	std::list<std::shared_ptr<Edge>> edges;
	std::list<std::shared_ptr<Face>> faces;

	edges.push_back(e1);
	edges.push_back(e2);
	edges.push_back(e3);
	edges.push_back(e4);
	edges.push_back(e5);
	edges.push_back(e6);

	faces.push_back(f1);
	faces.push_back(f2);
	faces.push_back(f3);
	faces.push_back(f4);

	Hull hull(
		edges,
		faces,
		myGeometry::findBottom(v1, v2, v3,v4),
		myGeometry::findTop(v1, v2, v3, v4),
		myGeometry::findLeft(v1, v2, v3, v4),
		myGeometry::findRight(v1, v2, v3, v4),
		(points[0] + points[1] + points[2] + points[3]) / 4.0);
	return hull;
}

// goes through the list of faces and determines which are visible or not
void myGeometry::checkFacesVisibility(const std::list<std::shared_ptr<Face>>& faces, std::shared_ptr<Vertex> bottom, std::shared_ptr<Vertex> top, std::shared_ptr<Vertex> left, std::shared_ptr<Vertex> right)
{
	std::list<std::shared_ptr<Face>>::const_iterator f = faces.begin();
	bool isVisible;
	while (f != faces.end())
	{
		isVisible = (*f)->checkVisibility(bottom);
		if (isVisible)
		{
			(*f)->setIsVisible(true);
		}
		else
		{
			isVisible = (*f)->checkVisibility(top);
			if (isVisible)
			{
				(*f)->setIsVisible(true);
			}
			else
			{
				isVisible = (*f)->checkVisibility(left);
				if (isVisible)
				{
					(*f)->setIsVisible(true);
				}
				else
				{
					isVisible = (*f)->checkVisibility(right);
					(*f)->setIsVisible(isVisible);
				}
			}
		}
		f++;
	}
}

std::list<std::shared_ptr<Edge>> myGeometry::findBridgeEdges(std::list<std::shared_ptr<Edge>>& edges)
{
	std::list<std::shared_ptr<Edge>> bridgeEdges;
	std::list<std::shared_ptr<Edge>>::iterator e = edges.begin();
	int count;
	while (e != edges.end())
	{
		count = 0;
		if ((*e)->getAdjFace(0)->getIsVisible())
		{
			count++;
		}
		else
		{
			count--;
		}

		if ((*e)->getAdjFace(1)->getIsVisible())
		{
			count++;
		}
		else
		{
			count--;
		}

		if (count == 0)
		{
			(*e)->setUsed(false);
			bridgeEdges.push_back(*e);
		}

		e++;
	}
	return bridgeEdges;
}

unsigned int myGeometry::isBelowPlane(std::shared_ptr<Vertex> v, Eigen::Vector3d onPlane, Eigen::Vector3d normal)
{
	Eigen::Vector3d diff = v->getLocation() - onPlane;
	double dot = normal.dot(diff);
	if (dot == 0.0)
	{
		return 0;
	}
	else if(dot < 0.0)
	{
		return -1;
	}
	else
	{
		return 1;
	}
}

Eigen::Vector3d myGeometry::projectIntoYZ(Eigen::Vector3d pos)
{
	Eigen::Vector3d newPos;
	newPos(0) = 0.0;
	newPos(1) = pos(1);
	newPos(2) = pos(2);
	return newPos;
}

bool myGeometry::areCounterClockwise(Eigen::Vector3d p1, Eigen::Vector3d p2, Eigen::Vector3d p3)
{
	double cross = p1(1) * (p2(2) - p3(2)) - p1(2) * (p2(1) - p3(1)) + (p2(1) * p3(2) - p2(2) * p3(1));
	double epsilon = 1e-10;
	if (cross > epsilon)
	{
		return true;
	}
	else
	{
		return false;
	}
}


/*
	I have already debugged this function. If there are three points colinear when projected to the YZ-plane,
	it will fail. But I believe that for random data it will be fine. But keep in mind this is an usual suspect.
*/
std::list<std::shared_ptr<Edge>> myGeometry::sortEdges(const std::list<std::shared_ptr<Edge>>& edges, Eigen::Vector3d center)
{
	std::list<std::shared_ptr<Edge>>::const_iterator e = edges.begin();
	Eigen::Vector3d yAxis(0.0, 1.0, 0.0);
	Eigen::Vector3d projVert1;
	Eigen::Vector3d projVert2;
	Eigen::Vector3d projEdge;
	std::list<std::shared_ptr<Edge>> upEdges;
	std::list<std::shared_ptr<Edge>> downEdges;
	std::list<double> upAngles;
	std::list<double> downAngles;
	std::list<double>::iterator angleIt;
	std::list<std::shared_ptr<Edge>>::iterator edgeIt;
	std::shared_ptr<Vertex> temp1, temp2;
	double angle;
	unsigned int i;
	while (e != edges.end())
	{
		projVert1 = myGeometry::projectIntoYZ((*e)->getEndPts(0)->getLocation());
		projVert2 = myGeometry::projectIntoYZ((*e)->getEndPts(1)->getLocation());
		if (myGeometry::areCounterClockwise(projVert1, projVert2, center)) // we want the edges to be clockwise oriented 
		{
			temp1 = (*e)->getEndPts(0);
			temp2 = (*e)->getEndPts(1);
			(*e)->setEndPts(temp2, temp1);
			projEdge = projVert1 - projVert2;
		}
		else
		{
			projEdge = projVert2 - projVert1;
		}
		angle = yAxis.dot(projEdge) / projEdge.norm();
		i = 0;
		if (projEdge(2) < 0)
		{
			angleIt = downAngles.begin();
			edgeIt = downEdges.begin();
			for (; i < downAngles.size(); i++)
			{
				if (angle > *angleIt)
				{
					break;
				}
				angleIt++;
				edgeIt++;
			}
			downAngles.insert(angleIt, angle);
			downEdges.insert(edgeIt, *e);
		}
		else
		{
			angleIt = upAngles.begin();
			edgeIt = upEdges.begin();
			for (; i < upAngles.size(); i++)
			{
				if (angle < *angleIt)
				{
					break;
				}
				angleIt++;
				edgeIt++;
			}
			upAngles.insert(angleIt, angle);
			upEdges.insert(edgeIt, *e);
		}
		e++;
	}
	upEdges.insert(upEdges.end(), downEdges.begin(), downEdges.end());
	return upEdges;
}

bool myGeometry::isPivot(std::shared_ptr<Vertex> v, Eigen::Vector3d normal, std::shared_ptr<Edge> current, std::shared_ptr<Edge> prev, std::shared_ptr<Edge> next, int &side)
{
	Eigen::Vector3d vLocation = v->getLocation();
	std::vector<std::shared_ptr<Vertex>> vertices;
	vertices.push_back(current->getEndPts(0));
	vertices.push_back(current->getEndPts(1));
	vertices.push_back(prev->getEndPts(0));
	vertices.push_back(prev->getEndPts(1));
	vertices.push_back(next->getEndPts(0));
	vertices.push_back(next->getEndPts(1));
	std::shared_ptr<Face> f = current->getAdjFace(0);
	vertices.push_back(f->getVertex(0));
	vertices.push_back(f->getVertex(1));
	vertices.push_back(f->getVertex(2));
	f = current->getAdjFace(1);
	vertices.push_back(f->getVertex(0));
	vertices.push_back(f->getVertex(1));
	vertices.push_back(f->getVertex(2));
	int k, lastK;
	k = myGeometry::isBelowPlane(vertices[0], vLocation, normal);
	for (unsigned int i = 1; i < vertices.size(); i++)
	{
		lastK = k;
		k = myGeometry::isBelowPlane(vertices[i], vLocation, normal);
		if (k * lastK < 0)
		{
			return false;
		}
	}
	side = k;
	return true;
}


void myGeometry::findFirstBridge(const std::list<std::shared_ptr<Edge>>& bridgeEdges1, const std::list<std::shared_ptr<Edge>>& bridgeEdges2, std::list<std::shared_ptr<Edge>>::const_iterator &pivotEdge, unsigned int &pivotIndex, Eigen::Vector3d &normal)
{
	std::shared_ptr<Edge> prev1, current1, next1;
	prev1 = bridgeEdges1.back();
	current1 = bridgeEdges1.front();
	next1 = *std::next(bridgeEdges1.begin());
	std::list<std::shared_ptr<Edge>>::const_iterator current2, prev2, next2;
	current2 = bridgeEdges2.begin();
	std::shared_ptr<Vertex> u1 = current1->getEndPts(0);
	std::shared_ptr<Vertex> u2 = current1->getEndPts(1);
	std::shared_ptr<Vertex> v;
	Eigen::Vector3d diff = u2->getLocation() - u1->getLocation();
	diff = diff / diff.norm();
	bool foundPivot = false;
	unsigned int i = 0;
	bool isPivot1, isPivot2;
	int side1, side2;
	while (!foundPivot)
	{
		if (current2 != bridgeEdges2.begin())
		{
			prev2 = std::prev(current2);
		}
		else
		{
			prev2 = std::prev(bridgeEdges2.end());
		}
		if (current2 != std::prev(bridgeEdges2.end()))
		{
			next2 = std::next(current2);
		}
		else
		{
			next2 = bridgeEdges2.begin();
		}
		v = (*current2)->getEndPts(i);
		normal = diff.cross(v->getLocation() - u1->getLocation());
		isPivot1 = myGeometry::isPivot(u1, normal, current1, prev1, next1, side1);
		isPivot2 = myGeometry::isPivot(u1, normal, *current2, *prev2, *next2, side2);
		if (isPivot1 && isPivot2 && (side1 == side2))
		{
			foundPivot = true;
		}
		else
		{
			if (i == 0)
			{
				i++;
			}
			else
			{
				i = 0;
				current2 = next2;
			}
		}
	}
	pivotEdge = current2;
	pivotIndex = i;
	if (side1 > 0)
	{
		normal = -1.0 * normal;
	}
}

Hull myGeometry::mergeHulls(Hull hull1, Hull hull2)
{
	myGeometry::checkFacesVisibility(hull1.getFaces(), hull2.getBottom(), hull2.getTop(), hull2.getLeft(), hull2.getRight());
	myGeometry::checkFacesVisibility(hull2.getFaces(), hull1.getBottom(), hull1.getTop(), hull1.getLeft(), hull1.getRight());
	std::list<std::shared_ptr<Edge>> bridgeEdges1 = myGeometry::findBridgeEdges(hull1.getEdges());
	std::list<std::shared_ptr<Edge>> bridgeEdges2 = myGeometry::findBridgeEdges(hull2.getEdges());
	bridgeEdges1 = myGeometry::sortEdges(bridgeEdges1, hull1.getCenter());
	bridgeEdges2 = myGeometry::sortEdges(bridgeEdges2, hull2.getCenter());
	std::list<std::shared_ptr<Edge>>::const_iterator brid = bridgeEdges1.begin();
	std::shared_ptr<Vertex> firstU = bridgeEdges1.front()->getEndPts(0);
	std::list<std::shared_ptr<Edge>>::const_iterator pivotEdge;
	unsigned int pivotIndex;
	Eigen::Vector3d normal, normal1, normal2, diff, newNormal, edge1, edge2;
	Eigen::Vector3d hullCenter = 0.5 * (hull1.getCenter() + hull2.getCenter());
	myGeometry::findFirstBridge(bridgeEdges1, bridgeEdges2, pivotEdge, pivotIndex, normal);
	std::shared_ptr<Vertex> firstV = (*pivotEdge)->getEndPts(pivotIndex);
	std::shared_ptr<Vertex> u = firstU;
	std::shared_ptr<Vertex> v = firstV;
	std::shared_ptr<Vertex> uNext, vNext;
	std::list<std::shared_ptr<Edge>>::const_iterator current1 = bridgeEdges1.begin();
	std::list<std::shared_ptr<Edge>>::const_iterator current2;
	double cos1, cos2;
	std::shared_ptr<Edge> newEdge, lastEdge;
	lastEdge.reset();
	std::shared_ptr<Face> newFace, lastFace;
	if (pivotIndex) // means pivotIndex is 1, so we move to the next edge;
	{
		if (pivotEdge != std::prev(bridgeEdges2.end()))
		{
			current2 = std::next(pivotEdge);
		}
		else
		{
			current2 = bridgeEdges2.begin();
		}
	}
	else
	{
		current2 = pivotEdge;
	}
	uNext = (*current1)->getEndPts(1);
	vNext = (*current2)->getEndPts(1);
	edge1 = uNext->getLocation() - u->getLocation();
	edge2 = vNext->getLocation() - v->getLocation();
	unsigned int nF = hull1.getFaces().size();
	unsigned int nE = hull1.getEdges().size();
	unsigned int i = 0;
	do
	{
		lastEdge = newEdge;
		lastFace = newFace;	
		normal1 = myGeometry::outwardNormal(u->getLocation(), v->getLocation(), uNext->getLocation(), hullCenter);
		normal2 = myGeometry::outwardNormal(u->getLocation(), v->getLocation(), vNext->getLocation(), hullCenter);
		cos1 = normal.dot(normal1) / (normal.norm() * normal1.norm());
		cos2 = normal.dot(normal2) / (normal.norm() * normal2.norm());
		std::cout << "cos1 = " << cos1 << "  , cos2 = " << cos2 << std::endl;
		if ((*current1)->getUsed())
		{
			cos1 = -1.0;
		}
		if ((*current2)->getUsed())
		{
			cos2 = -1.0;
		}
		if (cos1 >= cos2) // the neighbour of u is the winner
		{
			newEdge = std::make_shared<Edge>();
			newFace = std::make_shared<Face>();
			newEdge->setEndPts(uNext, v);
			newFace->setVertices(u, v, uNext);
			newFace->setEdges(*current1, lastEdge, newEdge);
			normal = normal1;
			newFace->setNormal(normal);
			if ((*current1)->getAdjFace(0)->getIsVisible())
			{
				(*current1)->setAdjFaces((*current1)->getAdjFace(1), newFace);
			}
			else
			{
				(*current1)->setAdjFaces((*current1)->getAdjFace(0), newFace);
			}
			(*current1)->setUsed(true);
			if (current1 != std::prev(bridgeEdges1.end()))
			{
				current1++;
			}
			else
			{
				current1 = bridgeEdges1.begin();
			}
			u = (*current1)->getEndPts(0);
			uNext = (*current1)->getEndPts(1);
		}
		else // the neighbour of v is the winner
		{
			newEdge = std::make_shared<Edge>();
			newFace = std::make_shared<Face>();
			newEdge->setEndPts(u, vNext);
			newFace->setVertices(u, v, vNext);
			newFace->setEdges(*current2, lastEdge, newEdge);
			normal = normal2;
			newFace->setNormal(normal);
			if ((*current2)->getAdjFace(0)->getIsVisible())
			{
				(*current2)->setAdjFaces((*current2)->getAdjFace(1), newFace);
			}
			else
			{
				(*current2)->setAdjFaces((*current2)->getAdjFace(0), newFace);
			}
			(*current2)->setUsed(true);
			if (current2 != std::prev(bridgeEdges2.end()))
			{
				current2++;
			}
			else
			{
				current2 = bridgeEdges2.begin();
			}
			v = (*current2)->getEndPts(0);
			vNext = (*current2)->getEndPts(1);
		}
		if (lastEdge.get())
		{
			lastEdge->setAdjFaces(lastFace, newFace);
		}
		hull1.addFace(newFace);
		hull1.addEdge(newEdge);
		i++;
		if (i > 20)
		{
			break;
		}
	} while ((u != firstU || v != firstV));
	std::cout << "exiting merging loop" << std::endl;
	newEdge->setAdjFace(newFace, 0);
	newEdge->setAdjFace(*std::next(hull1.getFaces().begin(), nF), 1);
	(*std::next(hull1.getFaces().begin(), nF))->setEdge(newEdge, 1);
	hull1.deleteVisibles();
	hull2.deleteVisibles();
	current2 = hull2.getEdges().begin();
	while (current2 != hull2.getEdges().end())
	{
		hull1.addEdge(*current2);
		current2++;
	}
	std::list<std::shared_ptr<Face>>::iterator f = hull2.getFaces().begin();
	while (f != hull2.getFaces().end())
	{
		hull1.addFace(*f);
		f++;
	}
	hull1.updateExtremities(hull2.getBottom(), hull2.getTop(), hull2.getLeft(), hull2.getRight());
	hull1.setCenter(hullCenter);
	return hull1;
}

void myGeometry::hullToEigen(Hull& hull, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
	unsigned int n = hull.getFaces().size();
	Eigen::Vector3d diff1;
	Eigen::Vector3d diff2;
	Eigen::MatrixXd MV(3 * n, 3);
	Eigen::MatrixXi MF(n, 3);
	std::list<std::shared_ptr<Face>>::const_iterator f = hull.getFaces().begin();
	unsigned int i = 0;
	while (f != hull.getFaces().end())
	{
		MF(i, 0) = 3*i;
		MF(i, 1) = 3*i + 1;
		MF(i, 2) = 3*i + 2;
		diff1 = (*f)->getVertex(1)->getLocation() - (*f)->getVertex(0)->getLocation();
		diff2 = (*f)->getVertex(2)->getLocation() - (*f)->getVertex(0)->getLocation();
		diff1 = diff1.cross(diff2);
		if (diff1.dot((*f)->getNormal()) > 1e-9)
		{
			MV.row(3 * i) = (*f)->getVertex(0)->getLocation();
			MV.row(3 * i + 1) = (*f)->getVertex(1)->getLocation();
			MV.row(3 * i + 2) = (*f)->getVertex(2)->getLocation();
		}
		else
		{
			MV.row(3 * i) = (*f)->getVertex(0)->getLocation();
			MV.row(3 * i + 1) = (*f)->getVertex(2)->getLocation();
			MV.row(3 * i + 2) = (*f)->getVertex(1)->getLocation();
		}
		i++;
		f++;
	}
	V = MV;
	F = MF;
}

std::vector<Eigen::Vector3d> myGeometry::randomSphere(Eigen::Vector3d center, double radius, unsigned int nHull, unsigned int nInside)
{
	std::default_random_engine generator1;
	std::default_random_engine generator2;
	std::default_random_engine generator3;
	double pi = std::atan(1.0) * 4.0;
	std::uniform_real_distribution<double> rDist(0.0, radius);
	std::uniform_real_distribution<double> thetaDist(0.0, 2 * pi);
	std::uniform_real_distribution<double> phiDist(0.0, pi);
	std::vector<Eigen::Vector3d> points(nHull + nInside);
	double r, theta, phi;
	Eigen::Vector3d vec;
	unsigned int i = 0;
	for (; i < nHull; i++)
	{
		theta = thetaDist(generator2);
		phi = phiDist(generator3);
		vec(0) = center(0) + radius * std::cos(theta) * std::sin(phi);
		vec(1) = center(1) + radius * std::sin(theta) * std::sin(phi);
		vec(2) = center(2) + radius * std::cos(phi);
		points[i] = vec;
	}
	for (; i < nHull + nInside; i++)
	{
		r = rDist(generator1);
		theta = thetaDist(generator2);
		phi = phiDist(generator3);
		vec(0) = center(0) + r * std::cos(theta) * std::sin(phi);
		vec(1) = center(1) + r * std::sin(theta) * std::sin(phi);
		vec(2) = center(2) + r * std::cos(phi);
		points[i] = vec;
	}
	return points;
}

bool myGeometry::compareX(Eigen::Vector3d p1, Eigen::Vector3d p2)
{
	return p1(0) < p2(0);
}

void myGeometry::sortByX(std::vector<Eigen::Vector3d>& points)
{
	std::sort(points.begin(), points.end(), myGeometry::compareX);
}