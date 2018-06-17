#pragma once

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <iterator>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <unsupported/Eigen/MatrixFunctions>

/// NAMESPACES
using namespace std;
using namespace Eigen;

/// RENDERING CONSTANTS
constexpr int numFrames = 1000;
constexpr int dim = 3;
constexpr bool testing = true;

/// TYPES
typedef double T;
typedef Matrix<T, dim, 1> V; // vertical
typedef Matrix<T, 1, dim> V_horizontal;
typedef Matrix<T, Dynamic, Dynamic> DM;

/// SIMULATION VARIABLES
struct {
    // originally defined
    int constraintIterations = 10;
    V gravity = V(T(0), T(-150), T(0));
    T deltaT = 1e-3;
    T damping = 0.05;
    T compressionStiffness = 1;
    T stretchingStiffness = 1;
    T bendingStiffness = 1;
    T restitution = 0.01; // used for friction
    vector<int> staticParticles = vector<int>();

    // calculated during simulation
    int numParticles = 0;
    int numFaces = 0;
    T originalVolume = 0;

    // not currently using
    T pressure = 1;
    T breakage = 1.5; // how much of the orig length before breaking the length
} sim;

// sphere variables
V_horizontal center = V_horizontal(5, 6.75, -3.5);
V_horizontal centerVelocity = V_horizontal(T(-14) / 300, 0, 0);
T radius = T(3);

/// FILES
string begFileName = "/Users/hbollar/Desktop/cis563/cis563-project1/Projects/example/cNEW_0_1stiff_1bend/";
string fileObj = "../../Vids/pbd/cow.obj";
string filePBDName = "pbdSim_";

/// CONSTANTS
constexpr T EPSILON = 1e-6;
constexpr T PI = 3.141592653589793238462643383279;

/// TYPE METHODS
T dot_horizontal(V_horizontal a, V_horizontal b) {
    return (a * b.transpose())(0, 0);
}
T dot_vertical(V a, V b) {
    return (a.transpose() * b)(0, 0);
}
T clamp(T a, T minVal, T maxVal){
    return max(minVal, min(maxVal, a));
}

