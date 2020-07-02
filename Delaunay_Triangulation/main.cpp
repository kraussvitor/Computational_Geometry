#include <iostream>
#include <glm/glm.hpp>
#include <Eigen/Dense>
#include <GLFW/glfw3.h>
#include <igl\readOBJ.h>
#include <igl\readOFF.h>
#include <igl\opengl\glfw\Viewer.h>
#include "geometry.h"

int main()
{
   
    unsigned int n = 500;
    std::vector<glm::dvec2> points = myGeometry::randomPoints(n);
    
    std::vector<glm::uvec3> triangles = myGeometry::incrementalTriangulation(points);
    std::list<glm::dvec2> hullPoints = myGeometry::incrementalConvexHull(points);

    std::vector<glm::uvec3> delaunay = myGeometry::makeDelaunay(points, triangles);
    std::vector<glm::dvec3> terrain = myGeometry::generateTerrain(points, triangles);
    Eigen::MatrixXd V = myGeometry::glmToEigen(terrain);
    Eigen::MatrixXi F = myGeometry::glmToEigen(delaunay);

    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);
    viewer.launch();

    return 0;
}