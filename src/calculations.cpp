#pragma once

#include "baseInfo.h"

namespace calculations {
    T calculateMeshVolume(vector<Face*> &faces, DM &p) {
        T volume = 0;

        V_horizontal v1 = V_horizontal();
        V_horizontal v2 = V_horizontal();
        V_horizontal v3 = V_horizontal();

        for(int i = 0; i < (int) faces.size(); ++i) {

            v1 = p.row(faces[i]->indices[0]);
            v2 = p.row(faces[i]->indices[1]);
            v3 = p.row(faces[i]->indices[2]);

            v1 = v1.cross(v2);
            volume += dot_horizontal(v1, v3);
        }

        return abs(volume / T(6));
    }

    void dampVelocities_simple(DM& v) {
        for (int i = 0; i < sim.numParticles; ++i) {
            for (int j = 0; j < dim; ++j) {
                v(i, j) *= 0.98;
            }
        }
    }

    void dampVelocities(DM &x, DM &v, DM &w) {
        // calculate center of mass / center of velocity / overall mass
        DM m = DM(dim, sim.numParticles);

        // since all weights are 1 or 0 --> will let all masses be the same
        // [need to update this if have diff weighted particles]
        m.row(0) = w.diagonal();
        m.row(1) = w.diagonal();
        m.row(2) = w.diagonal();

        DM xw = x * m;
        DM vw = v * m;
        V xwcm =  V( xw.col(0).sum(), xw.col(1).sum(), xw.col(2).sum() );
        V_horizontal vwcm = V_horizontal( vw.col(0).sum(), vw.col(1).sum(), vw.col(2).sum() );
        T mSum = m.sum();
        xwcm /= mSum;
        vwcm /= mSum;

        // computing global linear velocity (L = sum ri x (mivi); I = sum ?????)
        V L = V(0, 0, 0);
        DM I = DM(dim, dim); I = I.setZero(dim, dim);
        // temp variables for computations
        DM I_temp = DM(dim, dim); I_temp = I_temp.Zero(dim, dim);
        V r = V();
        V xTranspose = V();
        V vTranspose = V();
        for (int i = 0; i < sim.numParticles; ++i) {
            // since all weights are 1 or 0 --> will let all masses be the same as weights
            T mi = w(i, i);
            xTranspose = x.row(i).transpose();
            vTranspose = v.row(i).transpose();

            r = xTranspose - xwcm;
            L += r.cross(vTranspose) * mi;

            I_temp(0, 0) = T(0);  I_temp(0, 1) = -r(2); I_temp(0, 2) = r(1);
            I_temp(1, 0) = r(2);  I_temp(1, 1) = T(0);  I_temp(1, 2) = -r(0);
            I_temp(2, 0) = -r(1); I_temp(2, 1) = r(0);  I_temp(2, 2) = T(0);

            I += ((I_temp * I_temp.transpose()) * mi);
            I_temp = I_temp.Zero(dim, dim);
        }

        // calculate inverse of I
        T det = I.determinant();
        I.adjointInPlace();
        I /= det;

        V omega = I * L;
        for (int i = 0; i < sim.numParticles; ++i) {
            r = x.row(i).transpose() - xwcm;
            v.row(i) += sim.damping * (vwcm + omega.cross(r).transpose()  - v.row(i));
        }
    }

    void mathForFaceBendingConstraintCreation(Face* f1, Face* f2, vector<int> facesIndices, vector<V_horizontal> p, T& phi) {
        // using same p ordering as in update for this constraint
        V n1 = (p[2] - p[0]).cross(p[3] - p[0]); n1.normalize();
        V n2 = (p[3] - p[1]).cross(p[2] - p[1]); n2.normalize();

        T d = clamp(n1.dot(n2), T(-1), T(1));
        if (testing && isnan(d)) { throw; }

        phi = acos(d);
        if (testing && abs(acos(d) - phi) < EPSILON) { return; }
    }
}

namespace mapCalculations {
    bool mapContainsKeyValuePair(int i, int j, map<int, vector<int>>& map1) {
        // always want i to be smaller than j for map iterations [removes the issue of duplicating elements when checking
        if (j < i) { swap(i, j); }
        if (map1.find(i) == map1.end()) {
            return false;
        }

        vector<int> vectorAt = map1.at(i);
        return (find(vectorAt.begin(), vectorAt.end(), j) == vectorAt.end());
    }

    void addToMap(int i, int j, map<int, vector<int>>& map1) {
        if (j < i) { swap(i, j); }
        if (testing && mapContainsKeyValuePair(i, j, map1)) { throw; }

        vector<int> adding = (map1.find(i) == map1.end()) ? vector<int>() : map1.at(i);
        adding.push_back(j);

        map1[i] = adding;
    }
}