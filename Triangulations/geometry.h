#ifndef MY_GEOMETRY_H_INCLUDED
#define MY_GEOMETRY_H_INCLUDED

#include <glm/glm.hpp>
#include <vector>
#include <list>
#include <stack>
#include <Eigen/Dense>

namespace myGeometry
{
    void coutVector(const std::vector<glm::dvec2>& vector);
    void coutList(const std::list<glm::dvec2>::iterator start, const std::list<glm::dvec2>::iterator end);
    std::vector<glm::dvec2> randomPoints(unsigned int numPoints);
    std::vector<glm::dvec2> pointGrid(unsigned int numPoints);
    std::vector<glm::dvec2> randomParabola(unsigned int numPoints);
    std::vector<glm::dvec2> createParabola(unsigned int numPoints);
    bool compareX(const glm::dvec2 p1, const glm::dvec2 p2);
    void sortByX(std::vector<glm::dvec2>& vector);
    bool areCounterClockwise(glm::dvec2 p1, glm::dvec2 p2, glm::dvec2 p3);
    void incrementalStep(std::list<glm::dvec2>& points, glm::dvec2& newPoint, bool& isVisible, bool& lastIsVisible, std::list<glm::dvec2>::iterator& it);
    std::list<glm::dvec2> incrementalConvexHull(std::vector<glm::dvec2>& points);
    bool edgeIsValid(glm::dvec2 a, glm::dvec2 b, glm::dvec2 c, glm::dvec2 d, double epsilon);
    unsigned int findTriangle(std::vector<glm::uvec3> triangles, unsigned int v1, unsigned int v2, unsigned int& anchor);
    unsigned int findTriangle(std::vector<glm::uvec3> triangles, unsigned int v1, unsigned int v2, unsigned int& anchor, unsigned int current, bool& success);
    void incrementalTriangulationStep(std::vector<glm::dvec2> points, std::list<unsigned int>& hullIndices, std::vector<glm::uvec3>& triangles, unsigned int i);
    std::vector<glm::uvec3> incrementalTriangulation(std::vector<glm::dvec2>& points);
    void compareTriangulations(std::vector<glm::uvec3> triangles1, std::vector<glm::uvec3> triangles2);
    void addEdges(std::stack<glm::uvec2>& edges, glm::uvec3 triangle);
    std::vector<glm::uvec3> makeDelaunay(std::vector<glm::dvec2> points, std::vector<glm::uvec3> triangles);
    std::vector<glm::dvec3> generateTerrain(std::vector<glm::dvec2> points, std::vector<glm::uvec3> triangles);
    Eigen::MatrixXd glmToEigen(std::vector<glm::dvec3> points);
    Eigen::MatrixXi glmToEigen(std::vector<glm::uvec3> triangles);
}

#endif
