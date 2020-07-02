#ifndef MY_GEOMETRY_H_INCLUDED
#define MY_GEOMETRY_H_INCLUDED

#include <Eigen/Dense>
#include "Edge.h"
#include "Face.h"
#include "Vertex.h"
#include "Hull.h"
#include <vector>
#include <list>
#include <memory>

namespace myGeometry
{
	Eigen::Vector3d outwardNormal(Eigen::Vector3d p1, Eigen::Vector3d p2, Eigen::Vector3d p3, Eigen::Vector3d centroid);
	std::shared_ptr<Vertex> findBottom(std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2, std::shared_ptr<Vertex> v3);
	std::shared_ptr<Vertex> findBottom(std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2, std::shared_ptr<Vertex> v3, std::shared_ptr<Vertex> v4);
	std::shared_ptr<Vertex> findTop(std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2, std::shared_ptr<Vertex> v3);
	std::shared_ptr<Vertex> findTop(std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2, std::shared_ptr<Vertex> v3, std::shared_ptr<Vertex> v4);
	std::shared_ptr<Vertex> findLeft(std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2, std::shared_ptr<Vertex> v3);
	std::shared_ptr<Vertex> findLeft(std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2, std::shared_ptr<Vertex> v3, std::shared_ptr<Vertex> v4);
	std::shared_ptr<Vertex> findRight(std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2, std::shared_ptr<Vertex> v3);
	std::shared_ptr<Vertex> findRight(std::shared_ptr<Vertex> v1, std::shared_ptr<Vertex> v2, std::shared_ptr<Vertex> v3, std::shared_ptr<Vertex> v4);
	Hull initializeTriangle(const std::vector<Eigen::Vector3d> &points);
	Hull initializeTetrahedron(const std::vector<Eigen::Vector3d> &points);
	void checkFacesVisibility(const std::list<std::shared_ptr<Face>>& faces, std::shared_ptr<Vertex> bottom, std::shared_ptr<Vertex> top, std::shared_ptr<Vertex> left, std::shared_ptr<Vertex> right);
	std::list<std::shared_ptr<Edge>> findBridgeEdges(std::list<std::shared_ptr<Edge>>& edges);
	unsigned int isBelowPlane(std::shared_ptr<Vertex> v, Eigen::Vector3d onPlane, Eigen::Vector3d normal);
	Eigen::Vector3d projectIntoYZ(Eigen::Vector3d pos);
	bool areCounterClockwise(Eigen::Vector3d p1, Eigen::Vector3d p2, Eigen::Vector3d p3);
	std::list<std::shared_ptr<Edge>> sortEdges(const std::list<std::shared_ptr<Edge>>& edges, Eigen::Vector3d center);
	bool isPivot(std::shared_ptr<Vertex> v, Eigen::Vector3d normal, std::shared_ptr<Edge> current, std::shared_ptr<Edge> prev, std::shared_ptr<Edge> next, int& side);
	void findFirstBridge(const std::list<std::shared_ptr<Edge>>& bridgeEdges1, const std::list<std::shared_ptr<Edge>>& bridgeEdges2, std::list<std::shared_ptr<Edge>>::const_iterator& pivotEdge, unsigned int& pivotIndex, Eigen::Vector3d& normal);
	Hull mergeHulls(Hull hull1, Hull hull2);
	void hullToEigen(Hull& hull, Eigen::MatrixXd& V, Eigen::MatrixXi& F);
	std::vector<Eigen::Vector3d> randomSphere(Eigen::Vector3d center, double radius, unsigned int nHull, unsigned int nInside);
	bool compareX(Eigen::Vector3d p1, Eigen::Vector3d p2);
	void sortByX(std::vector<Eigen::Vector3d>& points);
}

#endif // !MY_GEOMETRY_H_INCLUDED

