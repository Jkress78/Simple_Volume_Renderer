/*
* Volume renderer header file
*/
#pragma once

#include <stdio.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <errno.h>

using namespace std;

/* Function Definitions */
vector<float> findnuv(float v1[], float v2[], float p[]);
void findmWCmCW(float vpn[], float vup[], float vrp[]);

void rayConstruction(int i, int j, float(&p0)[3], float(&v0)[3]);
float magnitude(float v[3]);
float dotProd(float v1[3], float v2[3]);
void swap(float (&t)[2]);

void computeShadingVolume();
int rayBoxIntersection(float(&p0)[3], float(&v0)[3], float(&ts)[2]);
int volumeRayTracing(float p0[3], float v[3], float ts[2]);
void triLinearInterpolation(float t, float p0[3], float v[3], float(&res)[2]);

/* Definition of image buffers */

#define LAYS 128
#define ROWS 128
#define COLS 128
unsigned char	CT[LAYS][ROWS][COLS]; /* a 3D array for CT data */
unsigned char	SHADING[LAYS][ROWS][COLS]; /* a 3D array for shading values */

#define IMG_ROWS 512
#define IMG_COLS 512
unsigned char	out_img[IMG_ROWS][IMG_COLS];

/* Camera parameters */
float VRP[3] = { 100.0, 64.0, 250.0 };
float VPN[3] = { -64.0, 0.0, -186.0 };
float VUP[3] = { 0.0, 1.0, 0.0 };

/* Image Plane Sizes */
float focal = 0.05;	/* 50 mm lens */
float xmin = -0.0175;	/* 35 mm "film" */
float ymin = -0.0175;
float xmax = 0.0175;
float ymax = 0.0175;

/* Light direction (unit length vector) */
float Light[3] = { 0.577, -0.577, -0.577 };
/* Light Intensity */
float Ip = 255.0;

float Mcw[4][4];
