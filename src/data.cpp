#pragma once

#include <Partio.h>
#include "baseInfo.h"

namespace simulationData {

    template <typename T>
    class Particles {
    public:
        int numParticles;
        DM positions;
        DM velocities;
        DM inverseMasses;

        Particles() :
                numParticles(0),
                positions(DM(0, 0)),
                velocities(DM(0, 0)),
                inverseMasses(DM()) {
        }

        ~Particles() {
            positions.resize(0, 0);
            velocities.resize(0, 0);
            inverseMasses.resize(0, 0);
        }

        void setUpSizes(const int &numParticles) {
            // check for resizing
            this->numParticles = numParticles;
            inverseMasses.resize(numParticles, numParticles);
            inverseMasses.block<numParticles, numParticles>(0,0).setIdentity();
            positions.setZero(numParticles, dim);
            velocities.setZero(numParticles, dim);
        }

        void resizeVelocityAndMassesFromPositions() {
            // assuming positions and numParticles is already set up properly
            velocities.setZero(numParticles, dim);
            inverseMasses.setZero(numParticles, numParticles);
        }

        void addPositions(vector<V>& newPositions) {
            positions.resize(numParticles + newPositions.size(), dim);

            for (int i = 0; i < (int)newPositions.size(); ++i) {
                positions(i + numParticles, 0) = newPositions[i][0];
                positions(i + numParticles, 1) = newPositions[i][1];
                positions(i + numParticles, 2) = newPositions[i][2];
            }
            numParticles += newPositions.size();
        }

        DM& getPositions() {
            return positions;
        }

        DM& getVelocities() {
            return velocities;
        }

        DM& getInverseMasses() {
            return inverseMasses;
        }

        void updateValues(DM& x, DM& v, DM& w) {
            positions = x;
            velocities = v;
            inverseMasses = w;
        }
    };

    template <typename T>
    class SegmentMesh {
    public:
        Particles<T> simParticles;
        vector<vector<int>> indices;

        SegmentMesh() :
                simParticles(Particles<T>()),
                indices(vector<vector<int>>()) {}

        ~SegmentMesh() {
            indices.clear();
            simParticles.~Particles<T>();
        }

        void setUpSizes(const int numParticles) {
            simParticles.setUpSizes(numParticles);
        }

        void addSegment(int indexA, int indexB) {
            indices[indexA].push_back(indexB);
            indices[indexB].push_back(indexA);
        }

        void addSegments(int startingA, vector<int>& attachedBs) {
            if (testing && startingA > (int)indices.size()) { throw; }

            if (startingA == (int)indices.size()) {
                indices.push_back(vector<int>());
            }

            for (int i = 0; i < (int)attachedBs.size(); ++i) {
                indices[startingA].push_back(attachedBs[i]);
            }
        }
    };

    class Face {
    public:
        vector<int> indices;
        vector<Face*> attachedFaces;
        int id;

        Face(vector<int> indicesInput, int inputId)
                : indices(indicesInput), attachedFaces(vector<Face*>()), id(inputId) {}
        ~Face() { indices.clear(); attachedFaces.clear(); }

        bool containsIndex(int i) {
            return (indices[0] == i || indices[1] == i || indices[2] == i);
        }

        bool containsIndices(int i, int j) {
            return containsIndex(i) && containsIndex(j);
        }

        bool containsOrderedIndices(int i, int j) {
            return ((indices[0] == i && indices[1] == j)
                    || (indices[1] == i && indices[2] == j)
                    || (indices[2] == i && indices[0] == j));
        }

        bool adjacentToFaceById(int faceId) {
            // returns true if my attached faces vector already contains this face
            for (int i = 0; i < (int) attachedFaces.size(); ++i) {
                if (attachedFaces[i]->id == faceId) {
                    return true;
                }
            }
            return false;
        }

        int getAdjacentFaceIndexByIndices(int sharedI, int sharedJ) {
            for (int i = 0; i < (int) attachedFaces.size(); ++i) {
                if (attachedFaces[i]->containsIndices(sharedI, sharedJ)) {
                    return i;
                }
            }
            if (testing) { throw; }
            return -1;
        }

        bool adjacentToFace(Face* adjacentFace) {
            return adjacentToFaceById(adjacentFace->id);
        }

        bool shouldBeAdjacentToFace(Face* face) {
            // if not already adjacent and contains a shared edge with this face
            // adds face to vector of attachedFaces vector and returns true; false otherwise
            if (!adjacentToFace(face)) {
                if (containsIndices(face->indices[0], face->indices[1])
                        || containsIndices(face->indices[1], face->indices[2])
                        || containsIndices(face->indices[2], face->indices[0])) {
                    attachedFaces.push_back(face);
                    face->shouldBeAdjacentToFace(this);

                    return true;
                }
            }
            return false;
        }

        bool removeAdjacentFaceById(int id) {
            if (!adjacentToFaceById(id)) {
                return false;
            }

            bool ret = false;
            vector<Face*> attachedNew = vector<Face*>();
            for (int i = 0; i < (int) attachedFaces.size() && !ret; ++i) {
                if (attachedFaces[i]->id != id) {
                    attachedNew.push_back(attachedFaces[i]);
                    ret = true;
                }
            }
            attachedFaces = attachedNew;
            return ret;
        }

        Face* getAdjacentById(int faceId) {
            if (testing && !(adjacentToFaceById(faceId))) { throw; }

            for (int i = 0; i < (int) attachedFaces.size(); ++i) {
                if (attachedFaces[i]->id == faceId) {
                    return (attachedFaces[i]);
                }
            }
            if (testing) { throw; }
            return nullptr;
        }

        int getNonListedIndex(int pi_index, int pj_index) {
            // returns index of point making face's triangle that is not one of the listed inputs
            for (int i = 0; i < (int) indices.size(); ++i) {
                if (indices[i] != pi_index && indices[i] != pj_index) {
                    return indices[i];
                }
            }
            if (testing) { throw; }
            return -1;
        }

        void resetIndices(int index0, int index1, int index2) {
            indices[0] = index0;
            indices[1] = index1;
            indices[2] = index2;
        }
    };

    template <typename T>
    class TriangleMesh {
    public:
        Particles<T> simParticles;
        vector<Face*> faces;

        TriangleMesh() :
                simParticles(Particles<T>()),
                faces(vector<Face*>()){}

        ~TriangleMesh() {
            faces.clear();
            simParticles.~Particles<T>();
        }

        void setUpSizes(const int& numParticles) {
            int input = (numParticles == -1) ? simParticles.numParticles : numParticles;
            simParticles.setUpSizes(input);
        }

        void calculateWeights() {
            // weight of particles set to 1 unless static position - if static, weight is 0
            // (maintains static position and prevents other constraints from adjusting them)

            simParticles.inverseMasses.setIdentity(simParticles.numParticles, simParticles.numParticles);
            for (int i = 0; i < (int) sim.staticParticles.size(); ++i) {
                simParticles.inverseMasses(sim.staticParticles[i], sim.staticParticles[i]) = T(0);
            }
        }

        void calculateAdjacentFaces() {
            DM check = DM::Zero(sim.numFaces, sim.numFaces);
            for (int i = 0; i < sim.numFaces; ++i) {
                for (int j = 0; j < sim.numFaces; ++j) {
                    if (i != j && check(i, j) != 1) {

                        check(i, j) = 1;
                        check(j, i) = 1;

                        if ( faces[i]->shouldBeAdjacentToFace(faces[j]) ) {
                            faces[j]->shouldBeAdjacentToFace(faces[i]);
                        }
                    }
                }
            }
        }
    };
}

