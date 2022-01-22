#include <stdio.h>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>

#include "invert.cpp"
#include "ppm.cpp"
#include "glm/glm.hpp"
#include "glm/ext.hpp"
#include "glm/gtx/intersect.hpp"

using namespace std;

void drawLine();
glm::vec3 raytracing(glm::vec3 eye, glm::vec3 ray, int count);

const int MAX_DEPTH = 3;
const glm::vec3 BLACK(0.f, 0.f, 0.f);
const float PASSED_TIME = 0.00001f;
const glm::vec3 PASSED_TIME_VEC(PASSED_TIME, PASSED_TIME, PASSED_TIME);

//Sphere struct
struct Sphere {
    string name;
    glm::vec3 pos;
    glm::vec3 scale;
    glm::vec3 color;
    //kAmbient, kDiffuse, kSpecular, kReflect
    glm::vec4 k;
    int specularExp;
};

//Light struct
struct Light {
    string name;
    glm::vec3 pos;
    glm::vec3 intensity;
};

//Plane struct
struct Plane {
    float near;
    glm::vec4 viewport;
    int resolution[2];
    glm::vec3 backgroundColor;
    glm::vec3 ambient;
};

struct Plane p;
char output[20];
vector<Sphere> sphereVec;
vector<Light> lightVec;
vector<glm::mat4> matrices;
vector<glm::mat4> inverseMatrices;
glm::vec3 eye(0.f, 0.f, 0.f);

int main(int argc, char **argv) {
    ifstream inputFile;
    inputFile.open(argv[1]);

    //Initialize data structs
    if(inputFile.is_open()) {
        string type;
        //Near and Viewport
        inputFile >> type >> p.near;
        for(int i = 0; i < 4; i++) {
            inputFile >> type >> p.viewport[i];
        }

        //Resolution
        inputFile >> type >> p.resolution[0] >> p.resolution[1];

        //Spheres
        inputFile >> type;
        while(type.compare("SPHERE") == 0) {
            struct Sphere s;
            inputFile >> s.name;
            for(int i = 0; i < 3; i++) {
                inputFile >> s.pos[i];
            }
            for(int i = 0; i < 3; i++) {
                inputFile >> s.scale[i];
            }
            for(int i = 0; i < 3; i++) {
                inputFile >> s.color[i];
            }
            for(int i = 0; i < 4; i++) {
                inputFile >> s.k[i];
            }
            inputFile >> s.specularExp;
            sphereVec.push_back(s);
            inputFile >> type;
        }

        //Lights
        while(type.compare("LIGHT") == 0) {
            struct Light l;
            inputFile >> l.name;
            for(int i = 0; i < 3; i++) {
                inputFile >> l.pos[i];
            }
            for(int i = 0; i < 3; i++) {
                inputFile >> l.intensity[i];
            }
            lightVec.push_back(l);
            inputFile >> type;
        }
        
        //Background Color and Ambient
        for(int i = 0; i < 3; i++) {
            inputFile >> p.backgroundColor[i];
        }
        inputFile >> type;
        for(int i = 0; i < 3; i++) {
            inputFile >> p.ambient[i];
        }

        //Output File Name
        inputFile >> type >> output;

        //Check file reading is done
        //If not, stop the raytrace and print error
        if(!inputFile.eof()) {
            cout << "File is invalid. Execute the raytrace" << endl;
            return 1;
        }

        //Printing Input Data
        drawLine();
        cout << "| Input Data \t\t" << argv[1] << "\t\t\t|" << endl;
        cout << "| Viewport \t\tN: " << p.near << "\tL: " << p.viewport[0] << "\tR: " << p.viewport[1] << "\tB: "
             << p.viewport[2] << "\tT: " << p.viewport[3] << "\t|" << endl;
        cout << "| Resolution \t\tX: " << p.resolution[0] << "\ty: " << p.resolution[1] << "\t\t\t\t|" << endl;
        drawLine();
        for(vector<int>::size_type i = 0; i < sphereVec.size(); i++) {
            cout << "| Sphere \t\t" << sphereVec[i].name << "\t\t\t\t\t|" << endl;
            cout << "| Position \t\tx: " << sphereVec[i].pos[0] << "\ty: " << sphereVec[i].pos[1] << "\tz: "
                 << sphereVec[i].pos[2] << "\t\t\t|" << endl;
            cout << "| Scale \t\tx: " << sphereVec[i].scale[0] << "\ty: " << sphereVec[i].scale[1] << "\tz: "
                 << sphereVec[i].scale[2] << "\t\t\t|" << endl;
            cout << "| Color \t\tr: " << sphereVec[i].color[0] << "\tg: " << sphereVec[i].color[1] << "\tb: "
                 << sphereVec[i].color[2] << "\t\t\t|" << endl;
            cout << "| K Coefficient\t\ta: " << sphereVec[i].k[0] << "\td: " << sphereVec[i].k[1] << "\ts: "
                 << sphereVec[i].k[2] << "\tr: " << sphereVec[i].k[3] << "\t\t|" << endl;
            cout << "| Specular Exponent \tn: " << sphereVec[i].specularExp << "\t\t\t\t\t|" << endl;
            drawLine();
        }
        for(vector<int>::size_type i = 0; i < lightVec.size(); i++) {
            cout << "| Light \t\t" << lightVec[i].name << "\t\t\t\t\t|" << endl;
            cout << "| Position \t\tx: " << lightVec[i].pos[0] << "\ty: " << lightVec[i].pos[1] << "\tz: " << lightVec[i].pos[2]
                 << "\t\t\t|" << endl;
            cout << "| Intensity \t\tx: " << lightVec[i].intensity[0] << "\ty: " << lightVec[i].intensity[1] << "\tz: "
                 << lightVec[i].intensity[2] << "\t\t\t|" << endl;
            drawLine();
        }
        cout << "| Background Color \tr: " << p.backgroundColor[0] << "\tg: " << p.backgroundColor[1] << "\tb: " 
             << p.backgroundColor[2] << "\t\t\t|" << endl;
        cout << "| Ambient \t\tr: " << p.ambient[0] << "\tg: " << p.ambient[1] << "\tb: " << p.ambient[2] << "\t\t\t|" << endl;
        drawLine();
    }
    inputFile.close();

    //Compute Inverse Matrices
    for(int sphIdx = 0; sphIdx < sphereVec.size(); sphIdx++) {
        glm::mat4 transMatScaleMat = glm::translate(glm::mat4(1.f), sphereVec[sphIdx].pos) * glm::scale(glm::mat4(1.f), sphereVec[sphIdx].scale);
        matrices.push_back(transMatScaleMat);
        double tempMatrix[4][4];
        for(int r = 0; r < 4; r++) {
            for(int c = 0; c < 4; c++) {
                tempMatrix[r][c] = transMatScaleMat[r][c];
            }
        }

        glm::mat4 inverseMatrix;
        double tempInvertMatrix[4][4];
        invert_matrix(tempMatrix, tempInvertMatrix);
        for(int r = 0; r < 4; r++) {
            for(int c = 0; c < 4; c++) {
                inverseMatrix[r][c] = tempInvertMatrix[r][c];
            }
        }
        inverseMatrices.push_back(inverseMatrix);
    }

    //Compute color for each pixel
    unsigned char *pixels;
    pixels = new unsigned char[3 * p.resolution[0] * p.resolution[1]];
    int pixIdx = 0;
    for(int r = p.resolution[1] - 1; r >= 0; r--) {
        for(int c = 0; c < p.resolution[0]; c++) {
            float uc = -p.viewport[1] + (2.f * p.viewport[1] * c / p.resolution[0]);
            float vr = -p.viewport[3] + (2.f * p.viewport[3] * r / p.resolution[1]);
            glm::vec3 point(uc, vr, -p.near);
            glm::vec3 ray = glm::normalize(point);
            glm::vec3 color = raytracing(eye, ray, 1);

            for(int i = 0; i < 3; i++) {
                pixels[pixIdx++] = glm::clamp(color[i], 0.f, 1.f) * 255;
            }
        }
    }
    save_imageP6(p.resolution[0], p.resolution[1], output, pixels);
    return 0;
}

glm::vec3 raytracing(glm::vec3 eye, glm::vec3 ray, int depth) {
    glm::vec3 inverseEye;
    glm::vec3 inverseRay;
    glm::vec3 inverseIntersectPos;
    glm::vec3 inverseIntersectNormal;
    glm::vec3 intersectPos;
    glm::vec3 intersectNormal;
    glm::vec3 localColor = BLACK;
    glm::vec3 reflectColor = BLACK;

    //Return black if depth is greater than MAX_DEPTH
    if (depth > MAX_DEPTH)
        return BLACK;

    bool hasIntersect = false;
    float closestDistance = INFINITY;
    Sphere* closestSphere = NULL;
    glm::mat4* closestMatrix = NULL;
    glm::mat4* closestInvMatrix = NULL;

    //Find closest hit point if exists
    for(int sphIdx = 0; sphIdx < sphereVec.size(); sphIdx++) {
        inverseEye = glm::vec3(inverseMatrices[sphIdx] * glm::vec4(eye, 1.f));
        inverseRay = glm::normalize(glm::vec3(inverseMatrices[sphIdx] * glm::vec4(ray, 0.f)));
        if(glm::intersectRaySphere(inverseEye, inverseRay, glm::vec3(0.f, 0.f, 0.f), 1, inverseIntersectPos, inverseIntersectNormal)) {
            hasIntersect = true;
            intersectPos = matrices[sphIdx] * glm::vec4(inverseIntersectPos, 1.f);
            float distance = glm::length(intersectPos);
            if(distance < closestDistance) {
                closestDistance = distance;
                closestSphere = &sphereVec[sphIdx];
                closestMatrix = &matrices[sphIdx];
                closestInvMatrix = &inverseMatrices[sphIdx];
            }
        }
    }

    //When no intersection and depth is 1, return background color
    //When no intersection and depth is not 1, return BLACK
    if(!hasIntersect) {
        return depth == 1 ? p.backgroundColor : BLACK;
    }

    inverseEye = glm::vec3(*closestInvMatrix * glm::vec4(eye, 1.f));
    inverseRay = glm::normalize(glm::vec3(*closestInvMatrix * glm::vec4(ray, 0.f)));
    glm::intersectRaySphere(inverseEye, inverseRay, glm::vec3(0.f, 0.f, 0.f), 1, inverseIntersectPos, inverseIntersectNormal);
    intersectPos = *closestMatrix * glm::vec4(inverseIntersectPos, 1.f);
    intersectNormal = glm::normalize(*closestMatrix * glm::vec4(inverseIntersectNormal, 0.f));

    //Find lights pointing intersection without blocked and Compute Diffuse and Specular
    for (int litIdx = 0; litIdx < lightVec.size(); litIdx++) {
        //Compute shadow rays
        bool isBlocked = false;
        for (int sphIdx = 0; sphIdx < sphereVec.size(); sphIdx++) {
            glm::vec3 lightPosOnSphMatrix(inverseMatrices[sphIdx] * glm::vec4(lightVec[litIdx].pos, 1.f));
            glm::vec3 intersectPosOnSphMatrix(inverseMatrices[sphIdx] * glm::vec4(intersectPos, 1.f));
            glm::vec3 lightDirectOnSphMatrix = glm::normalize(lightPosOnSphMatrix - intersectPosOnSphMatrix);
            float distanceToBlock;

            if (&sphereVec[sphIdx] != closestSphere && glm::intersectRaySphere(intersectPosOnSphMatrix, lightDirectOnSphMatrix, glm::vec3(0.f, 0.f, 0.f), 1, distanceToBlock) && distanceToBlock > 0.f) {
                isBlocked = true;
                break;
            }
        }

        if (!isBlocked) {
            glm::vec3 inverseLightPosition(*closestInvMatrix * glm::vec4(lightVec[litIdx].pos, 1.f));
            for (int i = 0; i < 3; i++) {
                //Compute Diffuse Color
                //Id * Kd * (N dot L)
                glm::vec3 L = inverseLightPosition - inverseIntersectPos;                
                float NdotL = glm::clamp(glm::dot(inverseIntersectNormal, glm::normalize(L)), 0.f, 1.f);
                float colorD = lightVec[litIdx].intensity[i] * closestSphere->k[1] * NdotL * closestSphere->color[i];
                
                //Compute Specular Color
                //Is * Ks * (H dot L)^n
                //Halfway vector H between L, V
                glm::vec3 H = (L + inverseRay) / glm::length(L + inverseRay);
                float HdotL = glm::clamp(glm::dot(H, inverseIntersectNormal), 0.f, 1.f);
                float colorS = lightVec[litIdx].intensity[i] * closestSphere->k[2] * glm::pow(HdotL, closestSphere->specularExp);

                localColor[i] += colorD + colorS;
            }
        }
    }

    //Compute Ambient Light
    for (int i = 0; i < 3; i++) {
        localColor[i] += p.ambient[i] * closestSphere->k[0] * closestSphere->color[i];
    }

    //Reflect Raytract -> Ray_rf(t) = P + vT
    //v = -2 * (N dot c) * N + c
    glm::vec3 reflectRay((-2.f * glm::dot(intersectNormal, ray) * intersectNormal) + ray);
    //Compute reflectPos start point when t = PASSED_TIME to avoiding false intersections
    glm::vec3 reflectPos(intersectPos + (reflectRay * PASSED_TIME_VEC));
    reflectColor = raytracing(reflectPos, glm::normalize(reflectRay), depth + 1);

    return localColor + (closestSphere->k[3] * reflectColor);
}

//Draw a line
void drawLine() {
    cout << "-----------------------------------------------------------------" << endl;
}