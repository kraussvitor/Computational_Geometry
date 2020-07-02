#include "geometry.h"
#include <random>
#include <iostream>
#include <list>
#include <iterator>
#include <chrono>
#include <functional>
#include <cfloat>
#include <algorithm>

void myGeometry::coutVector(const std::vector<glm::dvec2>& vector)
{
    for (unsigned int i = 0; i < vector.size(); i++)
    {
        std::cout << "x,y = " << vector[i].x << " , " << vector[i].y << std::endl;
    }
}

void myGeometry::coutList(const std::list<glm::dvec2>::iterator start, const std::list<glm::dvec2>::iterator end)
{
    std::list<glm::dvec2>::iterator it = start;
    std::cout << "\n ---------- The List ------------ \n" << std::endl;
    while (it != end)
    {
        std::cout << "  x,y = " << it->x << " , " << it->y << std::endl;
        it++;
    }
}

std::vector<glm::dvec2> myGeometry::pointGrid(unsigned int numPoints)
{
    std::vector<glm::dvec2> points(numPoints);
    for (unsigned int i = 0; i < numPoints; i++)
    {
        points[i].x = i;
        points[i].y = i + 1;
    }
    return points;
}


std::vector<glm::dvec2> myGeometry::randomPoints(unsigned int numPoints)
{
    std::vector<glm::dvec2> vector(numPoints);
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(-20.0, 20.0);
    for (unsigned int i = 0; i < numPoints; i++)
    {
        vector[i].x = distribution(generator);
        vector[i].y = distribution(generator);
    }
    return vector;
}

std::vector<glm::dvec2> myGeometry::randomParabola(unsigned int numPoints)
{
    std::vector<glm::dvec2> vector(2 * numPoints - 2);
    std::mt19937 generator;
    std::normal_distribution<double> distribution(0.0, 100.0);
    for (unsigned int i = 0; i < numPoints; i++)
    {
        vector[i].x = abs(distribution(generator));
        vector[i].y = pow(vector[i].x, 2.0);
    }
    for (unsigned int i = 1; i < numPoints - 1; i++)
    {
        vector[numPoints + i - 1] = 0.3 * (vector[i - 1] + vector[i] + vector[i + 1]);
    }
    return vector;
}

std::vector<glm::dvec2> myGeometry::createParabola(unsigned int numPoints)
{
    std::vector<glm::dvec2> vector(2 * numPoints);
    for (unsigned int i = 0; i < numPoints; i++)
    {
        vector[i].x = (double)i;
        vector[i].y = pow(vector[i].x, 2.0);
        vector[numPoints + i].x = -1.0 * (double)i;
        vector[numPoints + i].y = pow(vector[i].x, 2.0);
    }
    return vector;
}

bool myGeometry::compareX(const glm::dvec2 p1, const glm::dvec2 p2)
{
    if (p1.x != p2.x)
    {
        return p1.x < p2.x;
    }
    else
    {
        return p1.y < p2.y;
    }

}

void myGeometry::sortByX(std::vector<glm::dvec2>& vector)
{
    std::sort(vector.begin(), vector.end(), &myGeometry::compareX);
}

bool myGeometry::areCounterClockwise(glm::dvec2 p1, glm::dvec2 p2, glm::dvec2 p3)
{
    double epsilon = 1e-8;
    glm::dvec3 cross = glm::cross(glm::dvec3(p2 - p1, 0.0), glm::dvec3(p3 - p1, 0.0));
    if (cross.z > epsilon)
    {
        return true;
    }
    else if (cross.z < epsilon)
    {
        return false;
    }
    else
    {
        std::cout << "Points are colinear" << std::endl;
        return true;
    }

}

bool myGeometry::edgeIsValid(glm::dvec2 a, glm::dvec2 b, glm::dvec2 c, glm::dvec2 d, double epsilon)
{
    glm::dmat4x4 mat
    {
        a.x, a.y, std::pow(a.x, 2.0) + std::pow(a.y, 2), 1.0,
        b.x, b.y, std::pow(b.x, 2.0) + std::pow(b.y, 2), 1.0,
        c.x, c.y, std::pow(c.x, 2.0) + std::pow(c.y, 2), 1.0,
        d.x, d.y, std::pow(d.x, 2.0) + std::pow(d.y, 2), 1.0,
    };
    return !(glm::determinant(mat) > epsilon);
}

unsigned int myGeometry::findTriangle(std::vector<glm::uvec3> triangles, unsigned int v1, unsigned int v2, unsigned int& anchor)
{
    unsigned int vCount;
    for (int i = triangles.size() - 1; i >= 0; i--)
    {
        vCount = 0;
        if (triangles[i].x == v1)
        {
            vCount++;
        }
        else if (triangles[i].y == v1)
        {
            vCount++;
        }
        else if (triangles[i].z == v1)
        {
            vCount++;
        }

        if (triangles[i].x == v2)
        {
            vCount++;
        }
        else if (triangles[i].y == v2)
        {
            vCount++;
        }
        else if (triangles[i].z == v2)
        {
            vCount++;
        }

        if (vCount == 2)
        {
            if (triangles[i].x != v1 && triangles[i].x != v2)
            {
                anchor = triangles[i].x;
            }
            else if (triangles[i].y != v1 && triangles[i].y != v2)
            {
                anchor = triangles[i].y;
            }
            else
            {
                anchor = triangles[i].z;
            }
            return (unsigned int)i;
        }
    }
    return 0;
}

void myGeometry::incrementalStep(std::list<glm::dvec2>& points, glm::dvec2& newPoint, bool& isVisible, bool& lastIsVisible, std::list<glm::dvec2>::iterator& it)
{
    std::vector<std::list<glm::dvec2>::iterator> pivots;
    isVisible = !myGeometry::areCounterClockwise(points.back(), points.front(), newPoint);
    it = points.begin();
    while (it != std::prev(points.end()))
    {
        lastIsVisible = isVisible;
        isVisible = !myGeometry::areCounterClockwise(*it, *std::next(it), newPoint);
        if (isVisible != lastIsVisible)
        {
            pivots.push_back(it);
        }
        it++;
    }
    if (pivots.size() < 2)
    {
        pivots.push_back(std::prev(points.end()));
    }
    if (pivots[0] == points.begin())
    {
        if (newPoint.y < points.front().y)
        {
            points.erase(std::next(pivots[0]), pivots[1]);
            points.insert(pivots[1], newPoint);
        }
        else
        {
            points.erase(std::next(pivots[1]), points.end());
            points.insert(points.end(), newPoint);
        }
    }
    else
    {
        points.erase(std::next(pivots[0]), pivots[1]);
        points.insert(pivots[1], newPoint);
    }
}

std::list<glm::dvec2> myGeometry::incrementalConvexHull(std::vector<glm::dvec2>& points)
{
    myGeometry::sortByX(points);
    auto timeStart = std::chrono::high_resolution_clock::now();
    if (!myGeometry::areCounterClockwise(points[0], points[1], points[2]))
    {
        glm::dvec2 temp = points[1];
        points[1] = points[2];
        points[2] = temp;
    }
    std::list<glm::dvec2> pointList;
    pointList.push_front(points[2]);
    pointList.push_front(points[1]);
    pointList.push_front(points[0]);
    bool isVisible, lastIsVisible;
    std::list<glm::dvec2>::iterator it;
    for (unsigned int i = 3; i < points.size(); i++)
    {
        incrementalStep(pointList, points[i], isVisible, lastIsVisible, it);
    }
    auto timeStop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(timeStop - timeStart);
    long long time = duration.count();
    std::cout << "Time elapsed " << time << " microseconds " << std::endl;
    std::cout << "Number of Hull Points = " << pointList.size() << std::endl;
    return pointList;
}

void myGeometry::incrementalTriangulationStep(std::vector<glm::dvec2> points, std::list<unsigned int>& hullIndices, std::vector<glm::uvec3>& triangles, unsigned int i)
{
    std::vector<std::list<unsigned int>::iterator> pivots;
    bool isVisible, lastIsVisible;
    isVisible = !myGeometry::areCounterClockwise(points[hullIndices.back()], points[hullIndices.front()], points[i]);
    std::list<unsigned int>::iterator current = hullIndices.begin();
    std::list<unsigned int>::iterator next;
    while (current != hullIndices.end())
    {
        if (current == std::prev(hullIndices.end()))
        {
            next = hullIndices.begin();
        }
        else
        {
            next = std::next(current);
        }
        lastIsVisible = isVisible;
        isVisible = !myGeometry::areCounterClockwise(points[*current], points[*next], points[i]);
        if (isVisible)
        {
            triangles.push_back(glm::uvec3(*current, i, *next));
        }
        if (isVisible != lastIsVisible)
        {
            pivots.push_back(current);
        }
        current++;
    }
    if (pivots.size() < 2)
    {
        pivots.push_back(std::prev(hullIndices.end()));
    }
    if (pivots[0] == hullIndices.begin())
    {
        if (points[i].y <= points[hullIndices.front()].y)
        {
            hullIndices.erase(std::next(pivots[0]), pivots[1]);
            hullIndices.insert(pivots[1], i);
        }
        else
        {
            hullIndices.erase(std::next(pivots[1]), hullIndices.end());
            hullIndices.insert(hullIndices.end(), i);
        }
    }
    else
    {
        hullIndices.erase(std::next(pivots[0]), pivots[1]);
        hullIndices.insert(pivots[1], i);
    }
}

std::vector<glm::uvec3> myGeometry::incrementalTriangulation(std::vector<glm::dvec2>& points)
{
    myGeometry::sortByX(points);
    auto timeStart = std::chrono::high_resolution_clock::now();
    std::vector<glm::uvec3> triangles;
    std::list<unsigned int> hullIndices;
    if (!myGeometry::areCounterClockwise(points[0], points[1], points[2]))
    {
        glm::dvec2 temp = points[1];
        points[1] = points[2];
        points[2] = temp;
    }
    hullIndices.push_back(0);
    hullIndices.push_back(1);
    hullIndices.push_back(2);
    triangles.push_back(glm::uvec3(0, 1, 2));
    for (unsigned int i = 3; i < points.size(); i++)
    {
        myGeometry::incrementalTriangulationStep(points, hullIndices, triangles, i);
    }
    auto timeStop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(timeStop - timeStart);
    long long time = duration.count();
    std::cout << "-- incrementalTriangulation --\n Time elapsed " << time << " microseconds " << std::endl;
    unsigned int k = points.size() - hullIndices.size();
    std::cout << "Correct number of triangles = " << 2 * k + hullIndices.size() - 2 << std::endl;
    std::cout << "Number of computed triangles = " << triangles.size() << std::endl;
    std::cout << "Number of Hull Points = " << hullIndices.size() << std::endl;
    return triangles;
}

unsigned int myGeometry::findTriangle(std::vector<glm::uvec3> triangles, unsigned int v1, unsigned int v2, unsigned int& anchor, unsigned int current, bool& success)
{
    unsigned int vCount;
    success = false;
    for (int i = triangles.size() - 1; i >= 0; i--)
    {
        if (i != current)
        {
            vCount = 0;
            if (triangles[i].x == v1)
            {
                vCount++;
            }
            else if (triangles[i].y == v1)
            {
                vCount++;
            }
            else if (triangles[i].z == v1)
            {
                vCount++;
            }

            if (triangles[i].x == v2)
            {
                vCount++;
            }
            else if (triangles[i].y == v2)
            {
                vCount++;
            }
            else if (triangles[i].z == v2)
            {
                vCount++;
            }

            if (vCount == 2)
            {
                success = true;
                if (triangles[i].x != v1 && triangles[i].x != v2)
                {
                    anchor = triangles[i].x;
                }
                else if (triangles[i].y != v1 && triangles[i].y != v2)
                {
                    anchor = triangles[i].y;
                }
                else
                {
                    anchor = triangles[i].z;
                }
                return (unsigned int)i;
            }
        }
    }
    return 0;
}

void myGeometry::addEdges(std::stack<glm::uvec2>& edges, glm::uvec3 triangle)
{
    edges.push(glm::uvec2(triangle.x, triangle.y));
    edges.push(glm::uvec2(triangle.y, triangle.z));
    edges.push(glm::uvec2(triangle.z, triangle.x));
}

std::vector<glm::uvec3> myGeometry::makeDelaunay(std::vector<glm::dvec2> points, std::vector<glm::uvec3> triangles)
{
    auto timeStart = std::chrono::high_resolution_clock::now();
    std::stack<glm::uvec2> edges;
    std::vector<glm::uvec2> newTriangles;
    for (unsigned int i = 0; i < triangles.size(); i++)
    {
        myGeometry::addEdges(edges, triangles[i]);
    }
    unsigned int i, j, anchor1, anchor2;
    bool success1, success2;
    glm::uvec2 edge;
    glm::uvec3 t;
    double epsilon = 1e-5;
    while (!edges.empty())
    {
        edge = edges.top();
        i = myGeometry::findTriangle(triangles, edge.x, edge.y, anchor1, triangles.size() + 1, success1);
        j = myGeometry::findTriangle(triangles, edge.x, edge.y, anchor2, i, success2);
        t = triangles[i];
        if (success1 && success2)
        {
            if (myGeometry::edgeIsValid(points[t.x], points[t.y], points[t.z], points[anchor2], epsilon))
            {
                edges.pop();
            }
            else
            {
                if (myGeometry::areCounterClockwise(points[anchor1], points[anchor2], points[edge.x]))
                {
                    triangles[i] = glm::uvec3(anchor1, anchor2, edge.x);
                }
                else
                {
                    triangles[i] = glm::uvec3(anchor1, edge.x, anchor2);
                }

                if (myGeometry::areCounterClockwise(points[anchor1], points[anchor2], points[edge.y]))
                {
                    triangles[j] = glm::uvec3(anchor1, anchor2, edge.y);
                }
                else
                {
                    triangles[j] = glm::uvec3(anchor1, edge.y, anchor2);
                }
                myGeometry::addEdges(edges, triangles[i]);
                myGeometry::addEdges(edges, triangles[j]);
            }
        }
        else
        {
            edges.pop();
        }

    }
    auto timeStop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(timeStop - timeStart);
    long long time = duration.count();
    std::cout << "-- makeDelaunay -- \n Time elapsed " << time << " microseconds " << std::endl;
    return triangles;
}

void myGeometry::compareTriangulations(std::vector<glm::uvec3> triangles1, std::vector<glm::uvec3> triangles2)
{
    unsigned int inBoth = 0;
    glm::uvec3 t;
    for (unsigned int i = 0; i < triangles1.size(); i++)
    {
        t = triangles1[i];
        for (unsigned int j = 0; j < triangles2.size(); j++)
        {
            if (triangles2[j] == glm::uvec3(t.x, t.y, t.z) || triangles2[j] == glm::uvec3(t.z, t.x, t.y) || triangles2[j] == glm::uvec3(t.y, t.z, t.x))
            {
                inBoth++;
            }
        }
    }
    std::cout << "They have " << inBoth << " triangles in common, from a total of  " << triangles1.size() << std::endl;
}

std::vector<glm::dvec3> myGeometry::generateTerrain(std::vector<glm::dvec2> points, std::vector<glm::uvec3> triangles)
{
    std::vector<glm::dvec3> newPoints(points.size());
    double z;
    glm::uvec3 t;
    for (unsigned int i = 0; i < points.size(); i++)
    {
        z = 2.0 * std::sin(points[i].x)* std::cos(points[i].y);
        newPoints[i] = glm::dvec3(points[i].x, points[i].y, z);
    }
    return newPoints;
}

Eigen::MatrixXd myGeometry::glmToEigen(std::vector<glm::dvec3> points)
{
    unsigned int n = points.size();
    Eigen::MatrixXd M(n, 3);
    for (unsigned int i = 0; i < n; i++)
    {
        M(i, 0) = points[i].x;
        M(i, 1) = points[i].y;
        M(i, 2) = points[i].z;
    }
    return M;
}

Eigen::MatrixXi myGeometry::glmToEigen(std::vector<glm::uvec3> triangles)
{
    unsigned int n = triangles.size();
    Eigen::MatrixXi M(n, 3);
    for (unsigned int i = 0; i < n; i++)
    {
        M(i, 0) = triangles[i].x;
        M(i, 1) = triangles[i].y;
        M(i, 2) = triangles[i].z;
    }
    return M;
}


