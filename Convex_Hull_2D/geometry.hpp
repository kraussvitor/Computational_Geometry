#ifndef MY_GEOMETRY_H_INCLUDED
#define MY_GEOMETRY_H_INCLUDED

#include "../lib/glm/glm/glm.hpp"
#include <vector>
#include <list>

namespace myGeometry
{
    void coutVector(const std::vector<glm::dvec2>& vector);
    void coutList(const std::list<glm::dvec2>::iterator start, const std::list<glm::dvec2>::iterator end);
    std::vector<glm::dvec2> randomPoints(unsigned int numPoints);
    std::vector<glm::dvec2> randomParabola(unsigned int numPoints);
    std::vector<glm::dvec2> createParabola(unsigned int numPoints);
    bool compareX(const glm::dvec2 p1, const glm::dvec2 p2);
    void sortByX(std::vector<glm::dvec2>& vector);
    bool areCounterClockwise(glm::dvec2 p1, glm::dvec2 p2, glm::dvec2 p3);
    void incrementalStep(std::list<glm::dvec2>& points, glm::dvec2& newPoint, bool& isVisible, bool& lastIsVisible, std::list<glm::dvec2>::iterator& it);
    std::list<glm::dvec2> incrementalConvexHull(std::vector<glm::dvec2>& points);
    bool compareY(const glm::dvec2 p1, const glm::dvec2 p2);
    void sortByY(std::vector<glm::dvec2>& vector);
    double angleBetVec(glm::dvec2 vector1, glm::dvec2 vector2);
    bool compareAngleVec(const glm::dvec2 vector, const glm::dvec2 vertex, const glm::dvec2 p1, const glm::dvec2 p2);
    void sortByAngleVec(std::vector<glm::dvec2>& points, glm::dvec2 vector, glm::dvec2 vertex);
    glm::dvec2 findLargestAngleVec(std::vector<glm::dvec2> points, glm::dvec2 vector, glm::dvec2 vertex);
    std::vector<glm::dvec2> giftWrappingConvexHull(std::vector<glm::dvec2>& points);
    std::vector<glm::dvec2> grahamScanConvexHull(std::vector<glm::dvec2>& points);
    std::vector<glm::dvec2> divideConquerStep(std::vector<glm::dvec2> points);
    std::vector<glm::dvec2> divideConquerConvexHull(std::vector<glm::dvec2> points);
    std::vector<glm::dvec2> mergeHulls(std::vector<glm::dvec2> points1, std::vector<glm::dvec2> points2);
}

#endif