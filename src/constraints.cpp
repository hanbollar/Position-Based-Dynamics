#pragma once

#include "baseInfo.h"
#include "data.cpp"
#include "calculations.cpp"

class Constraint {
protected:
     bool remove;
     bool allowedToBreak;
public:
    Constraint(bool inputAllowedToBreak) : remove(false), allowedToBreak(inputAllowedToBreak) {}
    ~Constraint() {}

    virtual void update(DM& p, map<int, vector<int>>& broken) {}
    bool broken() {
        return (allowedToBreak && remove);
    }
};

class StretchConstraint : public Constraint {
protected:
    /// parental items
    // bool remove;
    // bool allowedToBreak;
public:
    int pi_index;
    int pj_index;

    T originalLength;
    T wi;
    T wj;

    StretchConstraint(T inputLen, T inputwi, T inputwj, int inputpi_index, int inputpj_index, bool inputAllowedToBreak)
            : Constraint(inputAllowedToBreak),
              pi_index(inputpi_index), pj_index(inputpj_index),
              originalLength(inputLen),
              wi(inputwi), wj(inputwj) {}

    ~StretchConstraint() { Constraint::~Constraint(); }

    virtual void update(DM& p, map<int, vector<int>>& broken) {
        // check if any fixed points - then dont do constraints calculations
        if (wi == 0 && wj == 0) { return; }
        if (allowedToBreak && (remove || mapCalculations::mapContainsKeyValuePair(pi_index, pj_index, broken))) { remove = true; return; }

        // f_ij = - youngsModulus / originalLength * (( len(xi - xj) / originalLength ) - 1) * (direction of xij)
        //      = - k * (||xi - xj||/L - 1) * n_ij
        // basically the spring force calculation except w1/(w1 + w2) or -w2/(w1 + w2) is E
        // here: delta_pi = - w1/(w1 + w2) * (( len(pi - p2) - d ) * (direction of pij)
        V_horizontal piMinpj = p.row(pi_index) - p.row(pj_index);
        T len = piMinpj.norm();
        if (testing && len == T(0)) {
            cout << "prevent dividing by zero: length is 0 error in calculations.cpp" << endl;
            throw;
        }

        T diff = len - originalLength;
        if (allowedToBreak && diff > 1.5 * originalLength) {
            mapCalculations::addToMap(pi_index, pj_index, broken);
        }

        V_horizontal updateValue = diff * piMinpj / ((wi + wj) * len);
        // compression or stretching stiffness
        T stiffness = (len < originalLength) ? sim.compressionStiffness : sim.stretchingStiffness;
        stiffness = 1 - pow((1 - stiffness), sim.constraintIterations);
        p.row(pi_index) += -wi * updateValue * (stiffness);
        p.row(pj_index) += wj * updateValue * (stiffness);

        if (testing) {
            for (int test_i = 0; test_i < 3; ++test_i) {
                if ( isnan(p.row(pi_index)(test_i)) || isnan(p.row(pj_index)(test_i)) ) { throw; }
            }
        }
    }
};

class FaceBendingConstraint : public Constraint {
protected:
    /// parental items
    // bool remove;
    // bool allowedToBreak;
public:
    // notEdge, edge, edge, notEdge
    vector<int> facesIndices;
    // weights for each of these vertices
    vector<T> w;
    // original angle between the faces
    T phi;

    FaceBendingConstraint(vector<int> inputFacesIndices, vector<T> inputW, T inputPhi, bool inputAllowedToBreak)
            : Constraint(inputAllowedToBreak), facesIndices(inputFacesIndices), w(inputW), phi(inputPhi) {}

    ~FaceBendingConstraint() { Constraint::~Constraint(); }

    virtual void update(DM& p, map<int, vector<int>>& broken) {
        if (w[0] == 0 || w[1] == 0 || w[2] == 0 || w[3] == 0) { return; }
        if (allowedToBreak && (remove || mapCalculations::mapContainsKeyValuePair(facesIndices[0], facesIndices[1], broken))) { remove = true; return; }
        if (allowedToBreak && (remove || mapCalculations::mapContainsKeyValuePair(facesIndices[1], facesIndices[2], broken))) { remove = true; return; }
        if (allowedToBreak && (remove || mapCalculations::mapContainsKeyValuePair(facesIndices[2], facesIndices[3], broken))) { remove = true; return; }

        V p0 = p.row(facesIndices[0]).transpose();
        V p1 = p.row(facesIndices[1]).transpose();
        V p2 = p.row(facesIndices[2]).transpose();
        V p3 = p.row(facesIndices[3]).transpose();

        V n1 = (p2 - p0).cross(p3 - p0); n1 /= n1.squaredNorm(); //n1.normalize();
        V n2 = (p3 - p1).cross(p2 - p1); n2 /= n2.squaredNorm(); //n2 *= T(-1); n2.normalize();

        V p32 = p3 - p2;
        T len_p32 = p32.norm();
        T len_p12 = p1.cross(p2).norm();
        T len_p13 = p1.cross(p3).norm();

        if (len_p12 < EPSILON  || len_p13 < EPSILON) { return; }

        DM q = DM(dim, 4);

        q.col(0) = len_p32 * n1;
        q.col(1) = len_p32 * n2;
        q.col(2) = (p0 - p3).dot(p32) / len_p32 * n1 + (p1 - p3).dot(p32) / len_p32 * n2;
        q.col(3) = (p2 - p0).dot(p32) / len_p32 * n1 + (p2 - p1).dot(p32) / len_p32 * n2;

        n1.normalize();
        n2.normalize();

        T d = clamp(n1.dot(n2), T(-1), T(1));
        if (testing && isnan(d)) { throw; }
        T testPhi = acos(d);
        if (testing && abs(testPhi - phi) < EPSILON) { return; }

        T weighting = 0;
        for (int i = 0; i < 4; ++i) {
            weighting += w[i] * q.col(i).squaredNorm();
        }
        if (weighting == T(0)) { return; }

        // note in below line: (1 - sim.stiffness) --> 1 - (1-k)^iterations bc more than one iter inside the
        // simulation loop [reduces error] (M. Muller et al. / Position Based Dynamics Section 3.3)
        T stiffness = sim.bendingStiffness;
        stiffness = 1 - pow((1 - stiffness), sim.constraintIterations);
        weighting = (testPhi - phi) / weighting * stiffness;
        if (n1.cross(n2).dot(p32) > 0) { weighting *= T(-1); }

        for (int i = 0; i < 4; ++i) {
            p.row(facesIndices[i]) +=  - w[i] * weighting * q.col(i).transpose();
            if (testing) {
                for (int test_i = 0; test_i < 3; ++test_i) {
                    if ( isnan(p.row(facesIndices[i])(test_i)) ) { throw; }
                }
            }
        }
    }
};

class BalloonVolumeConstraint : public Constraint {
protected:
    /// parental items
    // bool remove;
    // bool allowedToBreak;
public:
    BalloonVolumeConstraint(bool inputAllowedToBreak) : Constraint(inputAllowedToBreak) {}
    ~BalloonVolumeConstraint() { Constraint::~Constraint(); }

    virtual void update(DM& p, map<int, vector<int>>& broken) {
        // if allowedToBreak and broken -- return
    }
};

class CollisionConstraint {
public:
    // index of the position needing to be updated
    int index;
    bool ground;
    bool sphere;
    V_horizontal qc;
    V_horizontal nc;
    T w;

    CollisionConstraint(int inputIndex, bool inputGround, bool inputSphere, T inputWeight)
            : index(inputIndex), ground(inputGround), sphere(inputSphere), w(inputWeight) {}

    ~CollisionConstraint() { }

    void update(DM& p) {
        // collisions are an inequality constraint

        // zero length string - if p is within an object OR if p is below ground

        // collision is inequality constraint with stiffness = 1;

        // static: compute the closest surface point qc and surface normal nc at this location
        // continuous: we test for each vertix i the ray xi -> p.
        //      if the ray enters an object, we calculate the constraint between
        //      p and the entry point qc and the surface normal at this point, nc.
        // (diff between static and continuous is where the qc will be places)

        T Cp = 0;
        V_horizontal diff = V_horizontal(0, 0, 0);

        if (ground) {
            // GROUND
            nc = V_horizontal(0, 1, 0);
            qc = V_horizontal(p(index, 0), 0, p(index, 2)); // hardcoding qc for ground constraint
            diff = p.row(index) - qc;
            Cp = dot_horizontal(diff, nc);

            if (Cp < 0) {
                p.row(index) += -Cp * T(1) * nc; // - angle * k-collisionrestitution * in direction of normal //qc;
            }
        } else if (sphere) {
            // SPHERE
            nc = (p.row(index) - center); nc.normalize();
            qc = radius*nc + center; // position - center
            diff = p.row(index) - qc;
            Cp = dot_horizontal(diff, nc); // dot product

            // -angle * k-collisionrestitution * in direction of normal
            if (Cp < 0) {
                p.row(index) += (-Cp * T(1) * nc);
            }
        } else { if (testing) { throw; }}
    }
};

class Constraints {
public:
    vector<Constraint*> stretchBendVolume;
    vector<CollisionConstraint*> collisions;

    Constraints() : stretchBendVolume(vector<Constraint*>()),
                    collisions(vector<CollisionConstraint*>()) {};

    ~Constraints() { stretchBendVolume.clear(); collisions.clear(); }

    void update(DM& p, map<int, vector<int>>& broken) {
        for (Constraint* c : stretchBendVolume) {
            c->update(p, broken);
        }
        for (CollisionConstraint* c : collisions) {
            c->update(p);
        }
    }

    void createStretchConstraints(vector<Face*> &faces, DM &positions, DM& w, bool allowedToBreak) {
        DM check = DM(sim.numParticles, sim.numParticles);
        for (int i = 0; i < (int) faces.size(); ++i) {
            int size = faces[i]->indices.size();
            for (int j = 0; j < size; ++j) {
                int p1 = faces[i]->indices[j];
                int p2 = (j == size - 1) ? faces[i]->indices[0] : faces[i]->indices[j+1];

                if (check(p1, p2) != 1 && p1 != p2) {
                    if (testing && check(p2, p1) == 1) { throw; }

                    check(p1, p2) = 1;
                    check(p2, p1) = 1;

                    T len = ( positions.row(p1) - positions.row(p2) ).norm();
                    stretchBendVolume.push_back(new StretchConstraint(len, w(p1, p1), w(p2, p2), p1, p2, allowedToBreak) );
                }
            }
        }
    }

    void createFaceBendingConstraints(vector<Face*>& faces, DM& p, DM& w, bool allowedToBreak) {
        vector<int> edge;
        vector<int> notEdge;
        vector<int> facesIndices;
        vector<V_horizontal> positions;
        vector<T> weights;
        T phi = 0;

        DM check = DM::Zero(sim.numFaces, sim.numFaces);
        for (int i = 0; i < sim.numFaces; ++i) {
            for (int j = 0; j < sim.numFaces; ++j) {
                if (i != j && check(i, j) != 1 && faces[i]->adjacentToFaceById(j)) {
                    check(i, j) = 1;
                    check(j, i) = 1; // ensures i will always be less than j

                    vector<int> edge = vector<int>();
                    vector<int> notEdge = vector<int>();
                    // add in the first two indices for the shared edge using face j indices
                    for (int k = 0; k < dim; ++k) {
                        if (find(faces[i]->indices.begin(), faces[i]->indices.end(), faces[j]->indices[k]) != faces[i]->indices.end()) {
                            edge.push_back(faces[j]->indices[k]);
                        }
                    }
                    // add in the indices not on the shared edge
                    notEdge.push_back(faces[i]->getNonListedIndex(edge[0], edge[1]));
                    notEdge.push_back(faces[j]->getNonListedIndex(edge[0], edge[1]));
                    if (testing && edge.size() != 2) { throw; }
                    if (testing && notEdge.size() != 2) { throw; }

                    // index ordering for bending = nonEdge, edge, edge, nonEdge
                    weights = { w(notEdge[0], notEdge[0]),  w(edge[0], edge[0]), w(edge[1], edge[1]), w(notEdge[1], notEdge[1]) };
                    facesIndices = {notEdge[0], edge[0], edge[1], notEdge[1]};
                    positions = { p.row(facesIndices[0]), p.row(facesIndices[1]), p.row(facesIndices[2]), p.row(facesIndices[3]) };
                    phi = 0;
                    calculations::mathForFaceBendingConstraintCreation(faces[i], faces[j], facesIndices, positions, phi);

                    stretchBendVolume.push_back(new FaceBendingConstraint(facesIndices, weights, phi, allowedToBreak));
                }
            }
        }
    }

    void createBalloonVolumeConstraints(DM &p, vector<vector<int>>& faces, bool allowedToBreak) {
//        T Cp = 0;
//        V_horizontal v1 = V_horizontal();
//        V_horizontal v2 = V_horizontal();
//        V_horizontal v3 = V_horizontal();
//        for (int i = 0; i < (int) faces.size(); ++i) {
//            //have to save out the variables so eigen will allow the cross math
//            v1 = p.row(faces[i][0]);
//            v2 = p.row(faces[i][1]);
//            v3 = p.row(faces[i][2]);
//            Cp += dot_horizontal(v1.cross(v2), v3);
//        }
//        Cp -= sim.pressure * sim.originalVolume;
//
//        // this is an equality constraint
//
//        // since equality constraint always want to do volume constraint creation
    }

    void createCollisionConstraints(DM& p, DM& w) {
        collisions.clear();

        bool sphereCollision = false;
        bool groundCollision = false;
        for (int i = 0; i < sim.numParticles; ++i) {

            // GROUND
            V_horizontal qc = V_horizontal(0, 0, 0);
            V_horizontal nc = V_horizontal(0, 1, 0);
            T Cp = 0;
            V_horizontal diff = V_horizontal(0, 0, 0);

            qc = V_horizontal(p(i, 0), 0, p(i, 2)); // hardcoding qc bc just grounding constraint for now
            diff = p.row(i) - qc;
            Cp = dot_horizontal(diff, nc);

            if (Cp < 0) {
                groundCollision = true;
            }

            // SPHERE
//            nc = (p.row(i) - center); nc.normalize();
//            qc = radius*nc + center; // position - center
//            diff = p.row(i) - qc;
//            Cp = dot_horizontal(diff, nc); // dot product
//
//            if (Cp < 0) {
//                sphereCollision = true;
//            }

            // w != 0 is so dont create a constraint between all restricted points
            if (groundCollision || sphereCollision) {
                collisions.push_back(new CollisionConstraint(i, groundCollision, sphereCollision, w(i, i)));
            }

            groundCollision = false;
            sphereCollision = false;
        }
    }

    void updateVelocitiesOfCollisions(DM &v) {
        for (CollisionConstraint* c : collisions) {
            V_horizontal nc = c->nc;
            V_horizontal velo = v.row(c->index);

            // Reflect velocity in the direction of the collision normal.
            // [2 * k-restitution * lambert's law * normal]
            velo -= T(2) * sim.restitution * dot_horizontal(velo, nc) * nc;
            V_horizontal friction = -(velo - dot_horizontal(velo, nc) * nc);

            // frictional update
            v.row(c->index) += friction;
        }
    }

    void updateRigidBodies(DM& p, DM& x) {
//        for (CollisionConstraint* c : collisions) {
//            V_horizontal deltaP = p.row(c->index) - x.row(c->index);
//            deltaP /= (c->w * sim.deltaT);
//            // accel = (p.row(c->index) - x.row(c->index)) / (c->w * deltaT);
//            // v = v0 + at
//            // centerVelocity = centerVelocity + accel * deltaT;
//            // shortens to: deltaP * deltaT / c->w
//            // all weights are 1 or zero so fine
//            centerVelocity += (p.row(c->index) - x.row(c->index)) * sim.deltaT;
//            //cout<<x.row(c->index)<<" << x, p >>"<<p.row(c->index)<<" deltaP: "<<deltaP<<endl;
//            cout<<centerVelocity<<endl;
//        } // impulse for sphere position
    }

    void updateEdges(map<int, vector<int>>& broken, DM* x, DM* v, DM& w, DM& p, vector<Face*>& faces) {

        if (!(broken.size() > 0)) { return; }

        std::map<int, vector<int>>::iterator it = broken.begin();
        vector<int> edgePoints;
        int key;

        map<int, tuple<V_horizontal, V_horizontal>> newParticles = map<int, tuple<V_horizontal, V_horizontal>>();
        vector<Constraint*> addingConstraints = vector<Constraint*>();

        int numPositions = sim.numParticles;

        /// loop over all the key values
        while( it != broken.end()) {
            key = it->first;
            edgePoints = it->second;

            /// loop over each connected edgePoint to each key value
            for (int i = 0; i < (int) edgePoints.size(); ++i) {
                /// create edges DE, DF, and CE, CF [using indices]
                int d = key; // to match the notation in my notes
                int c = edgePoints[i];
                int e = numPositions;
                int f = numPositions + 1;
                numPositions += 2;

                /// store new particles x and v components
                // find midpoint of removed edge
                V_horizontal midpoint = (x->row(d) + x->row(c)) / T(2);
                V_horizontal aveVelocity = (v->row(d) + v->row(c)) / T(2);
                // create points E and F at midpoint with aveVelocity
                newParticles[e] = make_tuple(midpoint, aveVelocity);
                newParticles[f] = make_tuple(midpoint, aveVelocity);

                /// create and find indices for computations
                int indexFace0 = -1;
                int indexFace1 = -1;
                for (int k = 0; k < (int) faces.size() && indexFace0 == -1; ++k) {
                    if (faces[k]->containsIndices(c, d)) {
                        // on one of the correct faces
                        indexFace0 = k;
                    }
                }
                if (testing && indexFace0 == -1) { throw; }
                Face *face0 = faces[indexFace0];
                indexFace1 = face0->getAdjacentFaceIndexByIndices(c, d);
                if (testing && indexFace1 == -1) { throw; }
                if (testing && indexFace0 == indexFace1) { throw; }
                Face *face1 = faces[indexFace1];

                // must be triangular faces for our implementation
                if (testing && (face0->indices.size() != 3 || face1->indices.size() != 3)) { throw; }

                /// fill in a and b variables
                // a is the non-face index on face 1; b is the non-face index on face 0
                int a = -1;
                int b = -1;
                for (int index = 0; index < 3; ++index) {
                    if (face1->indices[index] != d && face1->indices[index] != c) {
                        a = face1->indices[index];
                    }
                    if (face0->indices[index] != d && face0->indices[index] != c) {
                        b = face0->indices[index];
                    }
                }
                if (testing && (a == -1 || b == -1)) { throw; }

                /// check swap for proper indexing based on my faces framework
                // for proper triangulation for our notation face 0 goes c->d, face 1 goes d->c
                if (faces[indexFace1]->containsOrderedIndices(c, d)) {
                    if (testing && faces[indexFace0]->containsOrderedIndices(c, d)) { throw; }

                    // face0 and face1 were labeled backwards for our indexing scheme so swap them
                    swap(indexFace0, indexFace1);
                }

                /// update faces indices in right order - CEB[update 0], AFC[update 1], ADF[new 2], BED[new 3]
                // dont update face0 and face1 indices yet
                Face *newFace2 = new Face({a, d, f}, faces.size());
                faces.push_back(newFace2);
                Face *newFace3 = new Face({b, e, d}, faces.size());
                faces.push_back(newFace3);

                /// update new faces' attached faces components for the one prev shared by face0 or face1 but not
                /// involved in other calculations
                // 2 takes 1's old face across DA edge connection [face 4]
                // 3 takes 0's old face across BD edge connection [face 5]
                int indexFace4 = face1->getAdjacentFaceIndexByIndices(d, a);
                int indexFace5 = face0->getAdjacentFaceIndexByIndices(b, d);
                if (testing && (indexFace5 == -1 || indexFace4 == -1)) { throw; }
                Face* face4 = faces[indexFace4];
                Face* face5 = faces[indexFace5];

                /// delete face0 and face1's old attached faces across BD and DA edges
                bool adjacentsRemovedCorrectly = true;
                // delete 0's old face across BD edge connection [face 5]
                // delete 5's shared face connection with 0
                // delete 1's old face across DA edge connection [face 4]
                // delete 4's shared face connection with 1
                adjacentsRemovedCorrectly &= face0->removeAdjacentFaceById(indexFace5);
                adjacentsRemovedCorrectly &= face5->removeAdjacentFaceById(indexFace0);
                adjacentsRemovedCorrectly &= face1->removeAdjacentFaceById(indexFace4);
                adjacentsRemovedCorrectly &= face4->removeAdjacentFaceById(indexFace1); // TODO: THIS PART MIGHT CAUSE ISSUES
                if (testing && !adjacentsRemovedCorrectly) { throw; }

                /// update face0 and face1 have proper position index values
                // resetting for simplicity [dont need to search for changed index]
                face0->resetIndices(c, e, b);
                face1->resetIndices(a, f, c);
                // 0 and 1 SHOULD STILL BE adjacents
                if (testing && !(face0->adjacentToFace(face1))) { throw; }

                /// make sure all other faces in this arrangement [including newly created ones]
                /// are properly adjacents
                bool adjacentsAddedCorrectly = true;
                // new adjacents = 2 3; 0 3; 2 1 [currently connecting these adjacent faces]
                adjacentsAddedCorrectly &= newFace2->shouldBeAdjacentToFace(newFace3);
                adjacentsAddedCorrectly &= face0->shouldBeAdjacentToFace(newFace3);
                adjacentsAddedCorrectly &= face1->shouldBeAdjacentToFace(newFace2);
                if (testing && !adjacentsAddedCorrectly) { throw; }
                                                                                        // TODO: ADD IN CASE FOR IF FACES G AND H DONT EXIST
                /// update mesh face bending constraints
                // remove bend between 05; 14; - done in loop outside of this iter loop
                // add bend between 0 3; 1 2; 3 5; 2 4;
                // face indices is in order of notEdge, edge, edge, notEdge
                T phi_03 = 0; T phi_12 = 0; T phi_35 = 0; T phi_24 = 0;
                int g = face4->getNonListedIndex(d, a);
                int h = face5->getNonListedIndex(b, d);

                // fill in needed weight values and point values for constraint creations
                // if point is created before this edge method then grab it, if not, then set to 1
                // num - 1 bc indexing
                T wa = (a > sim.numParticles - 1) ? T(1) : w(a,a);
                T wb = (b > sim.numParticles - 1) ? T(1) : w(b,b);
                T wc = (c > sim.numParticles - 1) ? T(1) : w(c,c);
                T wd = (d > sim.numParticles - 1) ? T(1) : w(d,d);
                T we = T(1); //new edge for this iteration
                T wf = T(1); //new edge for this iteration
                T wg = (g > sim.numParticles - 1) ? T(1) : w(g,g);
                T wh = (h > sim.numParticles - 1) ? T(1) : w(h,h);
                V_horizontal pa = (a > sim.numParticles - 1) ? get<0>(newParticles[a]) : p.row(a);
                V_horizontal pb = (b > sim.numParticles - 1) ? get<0>(newParticles[b]) : p.row(b);
                V_horizontal pc = (c > sim.numParticles - 1) ? get<0>(newParticles[c]) : p.row(c);
                V_horizontal pd = (d > sim.numParticles - 1) ? get<0>(newParticles[d]) : p.row(d);
                V_horizontal pe = get<0>(newParticles[e]); //new edge for this iteration
                V_horizontal pf = get<0>(newParticles[f]); //new edge for this iteration
                V_horizontal pg = (g > sim.numParticles - 1) ? get<0>(newParticles[g]) : p.row(g);
                V_horizontal ph = (h > sim.numParticles - 1) ? get<0>(newParticles[h]) : p.row(h);

                // facei, facej, faceIndices, phi to be filled in
                calculations::mathForFaceBendingConstraintCreation(face0, newFace3, {c, e, b, d}, {pc, pe, pb, pd}, phi_03);
                calculations::mathForFaceBendingConstraintCreation(face1, newFace2, {d, a, f, c}, {pd, pa, pf, pc}, phi_12);
                calculations::mathForFaceBendingConstraintCreation(newFace3, face5, {e, b, d, h}, {pe, pb, pd, ph}, phi_35);
                calculations::mathForFaceBendingConstraintCreation(newFace2, face4, {f, d, a, g}, {pf, pd, pa, pg}, phi_24);

                // inputFacesIndices, inputW, inputPhi
                addingConstraints.push_back(new FaceBendingConstraint({c, e, b, d}, {wc, we, wb, wd}, phi_03, true));
                addingConstraints.push_back(new FaceBendingConstraint({d, a, f, c}, {wd, wa, wf, wc}, phi_12, true));
                addingConstraints.push_back(new FaceBendingConstraint({e, b, d, h}, {we, wb, wd, wh}, phi_35, true));
                addingConstraints.push_back(new FaceBendingConstraint({f, d, a, g}, {wf, wd, wa, wg}, phi_24, true));

                /// update mesh stretch constraints
                // remove stretch from d,c - done in loop outside of this iter loop
                // add stretch to be; ce; cf; af; de; df;
                stretchBendVolume.push_back(new StretchConstraint((pb - pe).norm(), wb, we, b, e, true) );
                stretchBendVolume.push_back(new StretchConstraint((pc - pe).norm(), wc, we, c, e, true) );
                stretchBendVolume.push_back(new StretchConstraint((pc - pf).norm(), wc, wf, c, f, true) );
                stretchBendVolume.push_back(new StretchConstraint((pa - pf).norm(), wa, wf, a, f, true) );
                stretchBendVolume.push_back(new StretchConstraint((pd - pe).norm(), wd, we, d, e, true) );
                stretchBendVolume.push_back(new StretchConstraint((pd - pf).norm(), wd, wf, d, f, true) );

                // TODO FOR LATER : UPDATE VOLUME CONSTRAINTS

            } // end: going through all attached vertices of this key

            /// step the iterator
            ++it;
        }

        /// remove all the stretchConstraints and faceBendingConstraints that should be for this step
        stretchBendVolume.erase(remove_if(stretchBendVolume.begin(), stretchBendVolume.end(),
                               [](Constraint* c) { return c->broken(); }), stretchBendVolume.end()); // TODO: DOES THIS ACTUALLY DELETE THE CONSTRAINT* AND ALL INSTANCES?

        /// add in all the new stretchConstraints and faceBendingConstraints
        stretchBendVolume.insert(end(stretchBendVolume), begin(addingConstraints), end(addingConstraints));

        /// reclear the edges that were broken since they were fixed
        broken.clear();

        /// calculate updates for all inputs x, v, p, w, faces for new added positions [includes resizing]
        // copy over x into newX
        // copy over v into newV
        DM newX = DM::Zero(numPositions, dim);
        DM newV = DM::Zero(numPositions, dim);
        for (int i = 0; i < sim.numParticles; ++i) {
            newX.row(i) = x->row(i);
            newV.row(i) = v->row(i);
        }
        map<int, tuple<V_horizontal, V_horizontal>>::iterator newParticlesIter = newParticles.begin();
        while (newParticlesIter != newParticles.end()) {
            newX.row(newParticlesIter->first) = get<0>(newParticlesIter->second);
            newV.row(newParticlesIter->first) = get<1>(newParticlesIter->second);
        }
        sim.numParticles = numPositions;

        /// actually updating the inputs
        x->resize(sim.numParticles, dim); *x = newX;
        v->resize(sim.numParticles, dim); *v = newV;

        p.resize(sim.numParticles, dim);
        // below line is the same as mesh.calculate weights
        w.resize(sim.numParticles, sim.numParticles);
        for (int i = 0; i < (int) sim.staticParticles.size(); ++i) {
            w(sim.staticParticles[i], sim.staticParticles[i]) = T(0);
        }
    }

};

