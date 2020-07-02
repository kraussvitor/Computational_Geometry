#include "../include/geometry.hpp"
#include <random>
#include <iostream>
#include <list>
#include <iterator>
#include <chrono>
#include <functional>
#include <cfloat>

void myGeometry::coutVector(const std::vector<glm::dvec2>& vector)
{
    for(unsigned int i = 0; i < vector.size(); i++)
    {
        std::cout << "x,y = " << vector[i].x << " , " << vector[i].y << std::endl;
    }
}

void myGeometry::coutList(const std::list<glm::dvec2>::iterator start, const std::list<glm::dvec2>::iterator end)
{
    std::list<glm::dvec2>::iterator it = start;
    std::cout << "\n ---------- The List ------------ \n" << std::endl;
    while(it != end)
    {
        std::cout << "  x,y = " << it->x << " , " << it->y << std::endl;
        it++;
    }
}


std::vector<glm::dvec2> myGeometry::randomPoints(unsigned int numPoints)
{
    std::vector<glm::dvec2> vector(numPoints);
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(-1.0, 1.0);
    for(unsigned int i = 0; i < numPoints; i++)
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
    for(unsigned int i = 0; i < numPoints; i++)
    {
        vector[i].x = abs(distribution(generator));
        vector[i].y = pow(vector[i].x, 2.0);
    }
    for(unsigned int i = 1; i < numPoints-1; i++)
    {
        vector[numPoints + i - 1] = 0.3 * (vector[i-1] + vector[i] + vector[i+1]);
    }
    return vector;
}

std::vector<glm::dvec2> myGeometry::createParabola(unsigned int numPoints)
{
    std::vector<glm::dvec2> vector(numPoints);
    for(unsigned int i = 0; i < numPoints; i++)
    {
        vector[i].x = (double) i;
        vector[i].y = pow(vector[i].x, 2.0);
    }
    return vector;
}

bool myGeometry::compareX(const glm::dvec2 p1, const glm::dvec2 p2)
{
    if(p1.x != p2.x)
    {
        return p1.x < p2.x;
    }
    else
    {
        return p1.y < p2.y;
    }
    
}

bool myGeometry::compareY(const glm::dvec2 p1, const glm::dvec2 p2)
{
    if(p1.y != p2.y)
    {
        return p1.y < p2.y; 
    }
    else
    {
        return p1.x > p2.x;
    }
}

void myGeometry::sortByX(std::vector<glm::dvec2>& vector)
{
    std::sort(vector.begin(), vector.end(), &myGeometry::compareX);
}

void myGeometry::sortByY(std::vector<glm::dvec2>& vector)
{
    std::sort(vector.begin(), vector.end(), &myGeometry::compareY);
}


bool myGeometry::areCounterClockwise(glm::dvec2 p1, glm::dvec2 p2, glm::dvec2 p3)
{
    double epsilon = 1e-7;
    glm::dvec3 cross = glm::cross(glm::dvec3(p2-p1, 0.0), glm::dvec3(p3-p1, 0.0));
    if(cross.z > epsilon)
    {
        return true;
    }
    else if(cross.z < epsilon)
    {
        return false;
    }
    else
    {
        std::cout << "Points are colinear" << std::endl;
        return true;
    }
    
}

void myGeometry::incrementalStep(std::list<glm::dvec2>& points, glm::dvec2& newPoint, bool& isVisible, bool& lastIsVisible, std::list<glm::dvec2>::iterator& it)
{
    std::vector<std::list<glm::dvec2>::iterator> pivots;
    isVisible = !myGeometry::areCounterClockwise(points.back(), points.front(), newPoint);
    it = points.begin();
    while(it != std::prev(points.end()))
    {
        lastIsVisible = isVisible;
        isVisible = !myGeometry::areCounterClockwise(*it, *std::next(it), newPoint);
        if(isVisible != lastIsVisible)
        {
            pivots.push_back(it);
        }
        it++;
    }
    if(pivots.size() < 2)
    {
        pivots.push_back(std::prev(points.end()));
    }
    if(pivots[0] == points.begin())
    {
        if(newPoint.y < points.front().y)
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
    if(!myGeometry::areCounterClockwise(points[0], points[1], points[2]))
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
    for(unsigned int i = 3; i < points.size(); i++)
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

double myGeometry::angleBetVec(glm::dvec2 vector1, glm::dvec2 vector2)
{
    if(vector2.x == 0.0 && vector2.y == 0.0)
    {
        return 1.1;
    }
    else
    {
        return glm::dot(glm::normalize(vector1), glm::normalize(vector2));
    }
}

bool myGeometry::compareAngleVec(const glm::dvec2 vector, const glm::dvec2 vertex, const glm::dvec2 p1, const glm::dvec2 p2)
{
    return myGeometry::angleBetVec(vector, p1-vertex) > myGeometry::angleBetVec(vector, p2-vertex);
}

void myGeometry::sortByAngleVec(std::vector<glm::dvec2>& points, glm::dvec2 vector, glm::dvec2 vertex)
{
    auto compFun = std::bind(&myGeometry::compareAngleVec, vector, vertex, std::placeholders::_1, std::placeholders::_2);
    std::sort(points.rbegin(), points.rend(), compFun);
}

glm::dvec2 myGeometry::findLargestAngleVec(std::vector<glm::dvec2> points, glm::dvec2 vector, glm::dvec2 vertex)
{
    auto compFun = std::bind(&myGeometry::compareAngleVec, vector, vertex, std::placeholders::_1, std::placeholders::_2);
    //return *std::max_element(points.rbegin(), points.rend(), compFun);
    return *std::max_element(points.begin(), points.end(), compFun);
}

std::vector<glm::dvec2> myGeometry::giftWrappingConvexHull(std::vector<glm::dvec2>& points)
{
    glm::dvec2 bottom = *std::min_element(points.begin(), points.end(), &myGeometry::compareY);
    auto timeStart = std::chrono::high_resolution_clock::now();
    std::vector<glm::dvec2> hullPoints;
    glm::dvec2 lastEdge(-1.0, 0.0);
    glm::dvec2 newPoint;
    hullPoints.push_back(bottom);
    for(unsigned int i = 0; i < points.size(); i++)
    {
        newPoint = myGeometry::findLargestAngleVec(points, lastEdge, hullPoints.back());
        if(newPoint != hullPoints.front())
        {
            lastEdge = hullPoints.back() - newPoint;
            hullPoints.push_back(newPoint);
        }
        else
        {
            break;
        }
    }
    auto timeStop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(timeStop - timeStart);
    long long time = duration.count();
    std::cout << "\n \nTime elapsed " << time << " microseconds " << std::endl;
    std::cout << "Number of Hull points = " << hullPoints.size() << std::endl;
    return hullPoints;
}

std::vector<glm::dvec2> myGeometry::grahamScanConvexHull(std::vector<glm::dvec2>& points)
{
    glm::dvec2 bottom = *std::min_element(points.begin(), points.end(), &myGeometry::compareY);
    auto timeStart = std::chrono::high_resolution_clock::now();
    std::vector<glm::dvec2> hullPoints;
    hullPoints.push_back(bottom);
    myGeometry::sortByAngleVec(points, glm::dvec2(-1.0, 0.0), bottom);
    hullPoints.push_back(points.front());
    unsigned int i = 1; 
    while(i < points.size()-1)
    { 
        if(myGeometry::areCounterClockwise(*std::next(hullPoints.rbegin()), *hullPoints.rbegin(), points[i]))
        {
            hullPoints.push_back(points[i]);
            i++;
        }
        else
        {   
            hullPoints.pop_back();
        }
    }
    auto timeStop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(timeStop - timeStart);
    long long time = duration.count();
    std::cout << "Time elapsed " << time << " microseconds " << std::endl;
    std::cout << "Number of Hull points = " << hullPoints.size() << std::endl;
    return hullPoints;        
}

std::vector<glm::dvec2> myGeometry::mergeHulls(std::vector<glm::dvec2> points1, std::vector<glm::dvec2> points2)
{
    unsigned int numPoints1 = points1.size();
    unsigned int numPoints2 = points2.size();
    if(numPoints1 == 2)
    {
        if(points1[0].y < points1[1].y)
        {
            glm::dvec2 temp = points1[0];
            points1[0] = points1[1];
            points1[1] = temp;
        }
    }
    if(numPoints2 == 2)
    {
        if(points2[0].y > points2[1].y)
        {
            glm::dvec2 temp = points2[0];
            points2[0] = points2[1];
            points2[1] = temp;
        }
    }
    const std::vector<glm::dvec2>::iterator rightmost = std::max_element(points1.begin(), points1.end(), &myGeometry::compareX);
    const std::vector<glm::dvec2>::iterator leftmost = std::min_element(points2.begin(), points2.end(), &myGeometry::compareX);
    std::vector<glm::dvec2>::iterator alpha = rightmost;
    std::vector<glm::dvec2>::iterator beta = leftmost;
    std::vector<glm::dvec2>::iterator alphaPrev;
    std::vector<glm::dvec2>::iterator alphaNext;
    std::vector<glm::dvec2>::iterator betaPrev;
    std::vector<glm::dvec2>::iterator betaNext;
    bool traverseLeft = false;
    bool traverseRight = true;
    bool foundLeftTg = false;
    bool foundRightTg = false;
    unsigned int i = 0;
    while(!(foundRightTg && foundLeftTg) && i < numPoints1 + numPoints2)
    {
        if(std::next(beta) != points2.end())
        {
            betaNext = std::next(beta);
        }
        else
        {
            betaNext = points2.begin();
        }

        if(beta == points2.begin())
        {
            betaPrev = std::prev(points2.end());
        }
        else
        {
            betaPrev = std::prev(beta);
        }

        if(std::next(alpha) != points1.end())
        {
            alphaNext = std::next(alpha);
        }
        else
        {
            alphaNext = points1.begin();
        }

        if(alpha == points1.begin())
        {
            alphaPrev = std::prev(points1.end());
        }
        else
        {
            alphaPrev = std::prev(alpha);
        }
        foundLeftTg = myGeometry::areCounterClockwise(*alphaPrev, *alpha, *beta) && myGeometry::areCounterClockwise(*alphaNext, *alpha, *beta);
        foundRightTg = myGeometry::areCounterClockwise(*alpha, *beta, *betaPrev) && myGeometry::areCounterClockwise(*alpha, *beta, *betaNext);
        if(foundLeftTg && foundRightTg)
        {
            break;
        }
        else if(traverseRight)
        {
            if(foundRightTg)
            {
                traverseRight = false;
                traverseLeft = true;
            }
            else
            {
                beta = betaNext; // traverses counter clockwise
            }
        }
        else if(traverseLeft)
        {
            if(foundLeftTg)
            {
                traverseRight = true;
                traverseLeft = false;
            }
            else
            {
                alpha = alphaPrev; // traverses clockwise
            }            
        }
        i++;
    }
    // ---- end of 1st loop ----
    std::vector<glm::dvec2>::iterator lowerAlpha = alpha;
    std::vector<glm::dvec2>::iterator lowerBeta = beta;
    alpha = rightmost;
    beta = leftmost;
    foundLeftTg = false;
    foundRightTg = false;
    i = 0;
    // ---- start of 2nd loop ----
    while(!(foundRightTg && foundLeftTg) && i < numPoints1 + numPoints2)
    {
        if(std::next(beta) != points2.end())
        {
            betaNext = std::next(beta);
        }
        else
        {
            betaNext = points2.begin();
        }

        if(beta == points2.begin())
        {
            betaPrev = std::prev(points2.end());
        }
        else
        {
            betaPrev = std::prev(beta);
        }

        if(std::next(alpha) != points1.end())
        {
            alphaNext = std::next(alpha);
        }
        else
        {
            alphaNext = points1.begin();
        }

        if(alpha == points1.begin())
        {
            alphaPrev = std::prev(points1.end());
        }
        else
        {
            alphaPrev = std::prev(alpha);
        }
        foundLeftTg = myGeometry::areCounterClockwise(*alpha, *alphaPrev, *beta) && myGeometry::areCounterClockwise(*alpha, *alphaNext, *beta);
        foundRightTg = myGeometry::areCounterClockwise(*beta, *alpha, *betaPrev) && myGeometry::areCounterClockwise(*beta, *alpha, *betaNext);
        if(foundLeftTg && foundRightTg)
        {
            break;
        }
        else if(traverseRight)
        {
            if(foundRightTg)
            {
                traverseRight = false;
                traverseLeft = true;
            }
            else
            {
                beta = betaPrev; // traverses clockwise
            }
        }
        else if(traverseLeft)
        {
            if(foundLeftTg)
            {
                traverseRight = true;
                traverseLeft = false;
            }
            else
            {
                alpha = alphaNext; // traverses counterclockwise
            }            
        }
        i++;
    } // ---- end of 2nd loop
    std::vector<glm::dvec2>::iterator upperAlpha = alpha;
    std::vector<glm::dvec2>::iterator upperBeta = beta;
    std::vector<glm::dvec2> hullPoints;
    unsigned int uAlpha = std::distance(points1.begin(), upperAlpha);
    unsigned int lAlpha = std::distance(points1.begin(), lowerAlpha);
    unsigned int uBeta = std::distance(points2.begin(), upperBeta);
    unsigned int lBeta = std::distance(points2.begin(), lowerBeta);
    if(uAlpha < lAlpha)
    {
        hullPoints.insert(hullPoints.end(), upperAlpha, std::next(lowerAlpha));
    }
    else if(lAlpha < uAlpha)
    {
        hullPoints.insert(hullPoints.end(), upperAlpha, points1.end());
        hullPoints.insert(hullPoints.end(), points1.begin(), std::next(lowerAlpha));
    }
    else
    {
        hullPoints.insert(hullPoints.end(), *lowerAlpha);
    }
    if(lBeta < uBeta)
    {
        hullPoints.insert(hullPoints.end(), lowerBeta, std::next(upperBeta));
    }
    else if(uBeta < lBeta)
    {
        hullPoints.insert(hullPoints.end(), lowerBeta, points2.end());
        hullPoints.insert(hullPoints.end(), points2.begin(), std::next(upperBeta));
    }
    else
    {
        hullPoints.insert(hullPoints.end(), *lowerBeta);
    }
    return hullPoints;
}

std::vector<glm::dvec2> myGeometry::divideConquerStep(std::vector<glm::dvec2> points)
{
    unsigned int numPoints = points.size();
    if(numPoints == 2)
    {
        return points;
    }
    else if(numPoints == 3)
    {
        std::list<glm::dvec2> hullPoints;
        if(myGeometry::areCounterClockwise(points[0], points[1], points[2]))
        {
            return points;
        }
        else
        {
            std::vector<glm::dvec2> hullPoints(3);
            hullPoints[0] = points[0];
            hullPoints[1] = points[2];
            hullPoints[2] = points[1];
            return hullPoints;
        }
    }
    else
    {
        std::vector<glm::dvec2>::iterator middle = std::next(points.begin(), (int) std::floor(numPoints / 2));
        std::vector<glm::dvec2> firstHalf;
        firstHalf.insert(firstHalf.begin(), points.begin(), middle);
        std::vector<glm::dvec2> secondHalf;
        secondHalf.insert(secondHalf.begin(), middle, points.end());
        return myGeometry::mergeHulls(divideConquerStep(firstHalf), divideConquerStep(secondHalf));
    }
}

std::vector<glm::dvec2> myGeometry::divideConquerConvexHull(std::vector<glm::dvec2> points)
{
    myGeometry::sortByX(points);
    std::vector<glm::dvec2> hullPoints;
    auto timeStart = std::chrono::high_resolution_clock::now();
    hullPoints =  myGeometry::divideConquerStep(points);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - timeStart);
    long long time = duration.count();
    std::cout << "Time elapsed " << time << " microseconds " << std::endl;
    std::cout << "Number of Hull points = " << hullPoints.size() << std::endl;
    return hullPoints;
}









