#pragma once

#include "baseInfo.h"
#include "data.cpp"

using namespace simulationData;

void writeToFileParticles(const string& particleFile, Particles<T> &data){
    string usingFile = begFileName + particleFile + ".bgeo";
    cout<<usingFile<<endl;

    Partio::ParticlesDataMutable* parts = Partio::create();
    Partio::ParticleAttribute posH, vH, mH;
    mH = parts->addAttribute("m", Partio::VECTOR, 1);
    posH = parts->addAttribute("position", Partio::VECTOR, 3);
    vH = parts->addAttribute("v", Partio::VECTOR, 3);

    for (int i = 0; i < data.numParticles; ++i){
        int idx = parts->addParticle();
        float* m = parts->dataWrite<float>(mH, idx);
        float* p = parts->dataWrite<float>(posH, idx);
        float* v = parts->dataWrite<float>(vH, idx);
        m[0] = (T)i;
        for (int k = 0; k < 3; ++k) {
            p[k] = data.positions(i, k);
        }
        for (int k = 0; k < 3; ++k) {
            v[k] = data.velocities(i, k);
        }
    }

    Partio::write(usingFile.c_str(), *parts);
    parts->release();
}

void writeToFileOnePosition(const string& positionFile, V_horizontal centerPosition){
    string usingFile = begFileName + positionFile + ".bgeo";
    cout<<usingFile<<endl;

    Partio::ParticlesDataMutable* parts = Partio::create();
    Partio::ParticleAttribute posH;
    posH = parts->addAttribute("position", Partio::VECTOR, 3);

    int idx = parts->addParticle();
    float* p = parts->dataWrite<float>(posH, idx);
    for (int k = 0; k < 3; ++k) {
        p[k] = centerPosition(k);
    }

    Partio::write(usingFile.c_str(), *parts);
    parts->release();
}

void writeToFileSegment(const string& segmentsFile, const Particles<T> &data, const vector<Face*>& faces){
    ofstream myfile;
    myfile.open (begFileName + segmentsFile + ".poly");
    cout<<begFileName<<segmentsFile<<".poly"<<endl;

    myfile << "POINTS\n";
    for (int i = 0; i < data.numParticles; ++i) {
        // bc obj files index from 1 instead of 0
        myfile << i+1 <<" : "<<data.positions(i, 0)<<" "<<data.positions(i, 1)<<" "<<data.positions(i, 2)<<endl;
    }
    myfile << "POLYS\n";

    for (int i = 0; i < (int)faces.size(); ++i) {
        myfile << i+1 <<" : "<<(faces[i]->indices[0] + 1)<<" "<<(faces[i]->indices[1] + 1)<<" "<<(faces[i]->indices[2] + 1)<<" "<<(faces[i]->indices[0] + 1)<<endl;
    }
    myfile << "\nEND";
    myfile.close();
}

void writeToFileTriangles(const std::string& triangleFile, TriangleMesh<T> &data) {
    ofstream myfile;
    myfile.open (begFileName + triangleFile + ".obj");
    cout<<begFileName<<triangleFile<<".obj"<<endl;

    myfile << "g " << triangleFile << "\n";
    for (int i = 0; i < data.simParticles.numParticles; ++i) {
        myfile << "v ";
        for (int k = 0; k < 3; ++k) {
            myfile << data.simParticles.positions(i, k) << " ";
        }
        myfile << "\n";
    }
    myfile << "\n";
    for (int i = 0; i < (int)data.faces.size(); ++i) { //todo update here
        myfile << "f "<<(data.faces[i]->indices[0] + 1)<<" "<<(data.faces[i]->indices[1] + 1)<<" "<<(data.faces[i]->indices[2] + 1)<<" "<<"\n";
    }

    myfile.close();
}

void loadFromObj(Particles<T>& simParticles, vector<Face*>& faces) {
    string triangleFile = fileObj;

    string line;
    ifstream myfile (triangleFile);

    vector<V> newPositions = vector<V>();

    int numPositions = 0;
    if (myfile.is_open()) {
        while ( getline(myfile,line) ) {
            if (!(line.size() > 0)) { continue; }

            bool vertex = false;
            bool face = false;
            bool loadLine = false;
            if (line.at(0) == 'v' && line.at(1) == ' ') {
                vertex = true;
                loadLine = true;
            } else if (line.at(0) == 'f' && line.at(1) == ' ') {
                face = true;
                loadLine = true;
            }

            if (loadLine) {
                line = line.substr(2);
                int first = line.find(' ', 0);
                T firstVal = stod(line.substr(0, first));
                line = line.substr(first + 1);
                int second = line.find(' ', 0);
                T secondVal = stod(line.substr(0, second));
                line = line.substr(second + 1);
                T thirdVal = stod(line);

                if (vertex) {
                    /// loading positions
                    newPositions.push_back(V(firstVal, secondVal, thirdVal));
                    ++numPositions;
                } else if (face) {
                    // reached faces part of obj - if never been initialized yet, initialize
                    // note : (positions must be initialized before faces bc of obj format)

                    /// loading faces
                    // obj files index from 1 but we index from 0.
                    firstVal -= 1;
                    secondVal -= 1;
                    thirdVal -= 1;

                    faces.push_back( new Face({(int)firstVal, (int)secondVal, (int)thirdVal}, faces.size()) );
                }
            }
        }

        simParticles.addPositions(newPositions);
        myfile.close();
    } else {
        cout << "Unable to open file" << endl;
    }
}