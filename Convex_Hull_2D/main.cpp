#include <iostream>
#include <list>
#include "../include/geometry.hpp"
#include "../lib/glfw/include/GLFW/glfw3.h"
#include "../include/readShader.h"
#include <string>
#include "../include/myWindow.h"
#include "../lib/glad/include/glad/glad.h"
#include <functional>
#include <algorithm>
#include <iterator>

int main()
{
    unsigned int n = 100;
    std::vector<glm::dvec2> points = myGeometry::randomPoints(n); // generates points in the [-1,1] range
    //std::vector<glm::dvec2> points = myGeometry::randomParabola(n); // generates a parabola with interior points as well
    //std::vector<glm::dvec2> points = myGeometry::createParabola(n); // generates a parabola with only interior points

    std::cout << "\n \n ---------- Using Incremental ---------- \n \n" << std::endl;
    std::list<glm::dvec2> incHull = myGeometry::incrementalConvexHull(points);

    std::cout << "\n \n ---------- Using Gift Wrapping ---------- \n \n" << std::endl;
    std::vector<glm::dvec2> gwHull = myGeometry::giftWrappingConvexHull(points);

    std::cout << "\n \n ---------- Using Graham Scan ---------- \n \n" << std::endl;
    std::vector<glm::dvec2> gsHull = myGeometry::grahamScanConvexHull(points);

    std::cout << "\n \n ---------- Using Divide and Conquer ---------- \n \n" << std::endl;
    std::vector<glm::dvec2> dcHull = myGeometry::divideConquerConvexHull(points);
   
    std::vector<glm::dvec2> Hull = gsHull;

    // --- OpenGL Part ---
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    myWindow::myWindow window(800, 800);
    glfwMakeContextCurrent(window.getWindow());
    if(!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to load GLAD" << std::endl;
        glfwTerminate();
        return -1;
    }
    glViewport(0, 0, 800, 800);

    std::string vertexShaderStr = readShader::readShader("../shaders/vertexShader.txt");
    const char* vertexShaderSrc = vertexShaderStr.c_str();
    unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSrc, NULL);
    glCompileShader(vertexShader);
    int success;
    char infoLog[512];
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
    if(!success)
    {
        glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::VERTEX::COMPILATION_FAIL" << infoLog << std::endl;
    }
    std::string fragShaderStr = readShader::readShader("../shaders/fragShader.txt");
    const char* fragShaderSrc = fragShaderStr.c_str();
    unsigned int fragShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragShader, 1, &fragShaderSrc, NULL);
    glCompileShader(fragShader);
    glGetShaderiv(fragShader, GL_COMPILE_STATUS, &success);
    if(!success)
    {
        glGetShaderInfoLog(fragShader, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::FRAGMENT::COMPILATION_FAIL" << infoLog << std::endl;
    }

    unsigned int shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragShader);
    glLinkProgram(shaderProgram);
    glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
    if(!success)
    {
        glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
        std::cout << "PROGRAM::LINKING_FAIL" << infoLog << std::endl;
    }
    glDeleteShader(vertexShader);
    glDeleteShader(fragShader);

    unsigned int VAO, VBO;
    glCreateVertexArrays(1, &VAO);
    glBindVertexArray(VAO);
    glGenBuffers(1, &VBO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, points.size() * sizeof(glm::dvec2), points.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 2, GL_DOUBLE, GL_FALSE, 2 * sizeof(double), (void*)0);
    glEnableVertexAttribArray(0);

    unsigned int VAO2, VBO2;
    glCreateVertexArrays(1, &VAO2);
    glBindVertexArray(VAO2);
    glGenBuffers(1, &VBO2);
    glBindBuffer(GL_ARRAY_BUFFER, VBO2);
    glBufferData(GL_ARRAY_BUFFER, Hull.size() * sizeof(glm::dvec2), Hull.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 2, GL_DOUBLE, GL_FALSE, 2 * sizeof(double), (void*)0);
    glEnableVertexAttribArray(0);

    glUseProgram(shaderProgram);
    int colorUniformLocation = glGetUniformLocation(shaderProgram, "color");

    glEnable(GL_PROGRAM_POINT_SIZE);
    glLineWidth(10.0f);

    while(!glfwWindowShouldClose(window.getWindow()))
    {
        glClearColor(0.9f, 0.9f, 0.9f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        glUseProgram(shaderProgram);
        glUniform3f(colorUniformLocation, 0.0f, 0.0f, 0.0f);
        glBindVertexArray(VAO);
        glDrawArrays(GL_POINTS, 0, points.size());
        glUniform3f(colorUniformLocation, 0.047f, 0.87f, 0.56f);
        glBindVertexArray(VAO2);
        glDrawArrays(GL_LINE_LOOP, 0, Hull.size());
        glUniform3f(colorUniformLocation, 0.984f, 0.454f, 0.745f);
        glDrawArrays(GL_POINTS, 0, Hull.size());

        glfwSwapBuffers(window.getWindow());
        glfwPollEvents();
    }
    glfwTerminate();
    return 0;
}