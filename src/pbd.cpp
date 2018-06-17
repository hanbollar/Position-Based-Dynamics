#pragma once

#include "data.cpp"
#include "fileIO.cpp"
#include "constraints.cpp"

using namespace simulationData;

void pbdSimulation() {
    /// setup mesh values
    TriangleMesh<T> mesh = TriangleMesh<T>();

    loadFromObj(mesh.simParticles, mesh.faces);
    mesh.simParticles.resizeVelocityAndMassesFromPositions();

    /// simulation variables
    sim.numParticles = mesh.simParticles.numParticles;
    sim.numFaces = mesh.faces.size();
    sim.originalVolume = calculations::calculateMeshVolume(mesh.faces, mesh.simParticles.positions);
    sim.staticParticles = {}; // fill this with the indices of the static vertices

    mesh.calculateWeights();
    mesh.calculateAdjacentFaces();

    DM x = mesh.simParticles.positions;
    DM v = mesh.simParticles.velocities;
    DM w = mesh.simParticles.inverseMasses;

    // edge pairings only stored i,j s.t. i <= j
    // [removes the issue of checking based on duplicate edgings in the map - i,j and j,i]
    map<int, vector<int>> brokenEdges = map<int, vector<int>>();
    bool allowedToBreak = false;

    /// setup overall constraints
    Constraints allConstraints = Constraints();
    allConstraints.createStretchConstraints(mesh.faces, x, w, allowedToBreak);
    allConstraints.createFaceBendingConstraints(mesh.faces, x, w, allowedToBreak);
    cout<<"finished stretch, bending, volume constraints"<<endl;

    /// temp variables
    DM p = DM(sim.numParticles, sim.numParticles);  // used to store updated position locations
    DM ONESnx1 = DM::Ones(sim.numParticles, 1);     // empty variable used for calculations
    x.col(1) += ONESnx1 * T(4);                     // make mesh go above ground for beg of simulation // 7 plane 3.75 cow

    cout<<"beginning simulation looping"<<endl;
    //////////// BEGIN SIMULATION LOOP
    for (int simulationStep = 0; simulationStep < numFrames; ++simulationStep) {
        /// (5) for all vertices - update velocity by doing vi = vi + deltat*wi*fext(xi)
        ONESnx1 = DM::Ones(sim.numParticles, 1);     // make sure it resizes for edge breakage
        v = v + sim.deltaT * w * ONESnx1 * sim.gravity.transpose();

        /// (6) dampVelocites(v1,...vN)
        calculations::dampVelocities(x, v, w);

        /// (7) for all verticies i find projected point assuming no collisions pi = xi + deltat vi
        p = x + sim.deltaT * v;

        /// (8) for all velocities generate collision constraints (xi -> pi)
        allConstraints.createCollisionConstraints(p, w);

        /// (9) loop the solver the number of desired iterations projecting the constraints
        for (int i = 0; i < sim.constraintIterations; ++i) {
            allConstraints.update(p, brokenEdges);
        }

        /// (12) for all vertices update vi and xi for next overall loop simulation step
        v = (p - x) / sim.deltaT;
        x = p;

        /// (16) velocities update (for friction etc - only on particles involved in this iterations collisions)
        // since collision constraints are the last to be updated in allConstraints update - dont need to reconstrain here
        allConstraints.updateVelocitiesOfCollisions(v);
        //allConstraints.updateRigidBodies(p, x);
        //allConstraints.updateEdges(brokenEdges, &x, &v, w, p, mesh.faces);

        center += centerVelocity; // moving sphere for collisions

        /// - here print out files of current simulation step
        mesh.simParticles.updateValues(x, v, w);
        //writeToFileParticles(string(filePBDName + to_string(simulationStep)), mesh.simParticles);
        //writeToFileSegment(string(filePBDName + to_string(simulationStep)), mesh.simParticles, mesh.faces);
        writeToFileTriangles(filePBDName + to_string(simulationStep), mesh);
        //writeToFileOnePosition(string("sphereCenter_" + to_string(simulationStep)), center); // for debugging sphere location
    }
    //////////// END SIMULATION LOOP
    cout << "simulation completed" << endl;
}