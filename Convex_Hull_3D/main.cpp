#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <list>
#include "myGeometry.h"
#include "Edge.h"
#include "Face.h"
#include "Vertex.h"
#include "Hull.h"
#include <memory>
#include <GLFW/glfw3.h>
#include <igl/opengl/glfw/Viewer.h>

int main()
{
	std::vector<Eigen::Vector3d> points1;
	points1.push_back(Eigen::Vector3d(-10, 0.0, 25.0));
	points1.push_back(Eigen::Vector3d(0.0, 50.0, 0.0));
	points1.push_back(Eigen::Vector3d(0.0, -50.0, 0.0));
	points1.push_back(Eigen::Vector3d(0.0, 0.0, 50.0));

	std::vector<Eigen::Vector3d> points2;
	points2.push_back(Eigen::Vector3d(50.0, 40.0, 0.0));
	points2.push_back(Eigen::Vector3d(50.0, -40.0, 0.0));
	points2.push_back(Eigen::Vector3d(50.0, 0.0, 30.0));
	points2.push_back(Eigen::Vector3d(55.0, 0.0, 15.0));

	std::vector<Eigen::Vector3d> points3;
	points3.push_back(Eigen::Vector3d(65.0, 15.0, 1.0));
	points3.push_back(Eigen::Vector3d(65.0, -15.0, 1.0));
	points3.push_back(Eigen::Vector3d(65.0, 0.0, 16.0));


	Hull hull1 = myGeometry::initializeTetrahedron(points1);
	Hull hull2 = myGeometry::initializeTetrahedron(points2);
	Hull hull3 = myGeometry::initializeTriangle(points3);

	std::cout << "--------- first call -------------" << std::endl;
	Hull merged = myGeometry::mergeHulls(hull1, hull2);
	std::cout << "--------- second call -------------" << std::endl;
	Hull merged2 = myGeometry::mergeHulls(merged, hull3);
	
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	myGeometry::hullToEigen(merged2, V, F);

	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(V, F);
	viewer.launch();
	
	return 0;
}