/*
* Volume render
*/

#include "volume_renderer.h"

using namespace std;

int main() {
	FILE* infid, * outfid;
	int n;
	int i, j;
	errno_t err;
	if ((err = fopen_s(&infid, "./smallHead.den", "rb")) != 0) {
		cout << "OPEN CT DATA FILE ERROR!" << endl;
		exit(1);
	}
	for (i = 0; i < LAYS; i++) {
		n = fread(&CT[i][0][0], sizeof(char), ROWS * COLS, infid);
		if (n < ROWS * COLS * sizeof(char)) {
			cout << "READ CT DATA LAYER " << i << " ERROR" << endl;
			exit(1);
		}
	}
	
	/* Initialize Global data structures*/
	float P0[3] = { VRP[0], VRP[1], VRP[2] };
	float v[3];
	findmWCmCW(VPN, VUP, VRP);
	computeShadingVolume();
	
	/*=============================================
		Main Ray Tracing Volume Rendering Part
	=============================================*/

	for (i = 0; i < IMG_ROWS; i++) {
		for (j = 0; j < IMG_COLS; j++) {

			//Construct the ray V
			rayConstruction(i, j, P0, v);

			float ts[2]; //Stores the points t0 and t1

			n = rayBoxIntersection(P0, v, ts);  //n is the number of intersections found

			if (n == 2) {
				if (ts[0] > ts[1]) //if the intersections are not in the array smallest->largest switch their order
					swap(ts);
				out_img[i][j] = volumeRayTracing(P0, v, ts); //output the calculated shading value to the image buffer

			}
		}
	}

	/* Save the output image */
	if ((err = fopen_s(&outfid, "outimage28.raw", "wb")) != 0) {
		cout << "OPEN OUTPUT FILE ERROR" << endl;
		exit(1);
	}
	n = fwrite(out_img, sizeof(char), IMG_ROWS * IMG_COLS, outfid);
	if (n < IMG_ROWS * IMG_COLS * sizeof(char)) {
		cout << "WRITE OUTPUT IMAGE ERROR!" << endl;
		exit(1);
	}
	fclose(infid);
	fclose(outfid);
	return 0;
}

/* swap(float array[2]): 
*  Input: an array of size 2
*  Output: none
*  Description: swaps the position of the values in the inputed array
*/
void swap(float(&t)[2]) {
	float temp = t[0];
	t[0] = t[1];
	t[1] = temp;
}

/* computeShadingVolume():
*  Input: none
*  Output: none
*  Descriptiom: calculates the shading value for each point in the CT
*  data file and stores the calculated values in a 3D array
*/
void computeShadingVolume() {
	float dx;
	float dy;
	float dz;
	float N[3];
	float nMag;
	float epsilon = 20;
	int I;

	for (int z = 1; z < (LAYS - 1); z++) {
		for (int y = 1; y < (ROWS - 1); y++) {
			for (int x = 1; x < (COLS - 1); x++) {

				//Compute the Partial Derivative at [z, y, x]
				dx = 0.5 * (CT[z][y][x + 1] - CT[z][y][x - 1]);
				dy = 0.5 * (CT[z][y + 1][x] - CT[z][y - 1][x]);
				dz = 0.5 * (CT[z + 1][y][x] - CT[z - 1][y][x]);

				N[0] = dx;
				N[1] = dy;
				N[2] = dz;

				if ((nMag = magnitude(N)) < epsilon)
					I = 0;
				else
					I = Ip * (dotProd(N, Light));

				SHADING[z][y][x] = I;
			}
		}
	}

}


/* rayConstruction(int i, int j, float(&p0)[3], float(&v0)[3]):
*  Input: current points of the image plane (i,j), the camara origin p0[3], and the var for the ray to be stored in v0[3]
*  Output: none
*  Descriptiom: constructs a ray that passes through the point (i,j) on the image plane starting at p0
*/
void rayConstruction(int i, int j, float(&p0)[3], float(&v0)[3]) {
	//Convert point in image buffer to camera coords
	float x = ((xmax - xmin) * j) / (IMG_COLS - 1) + xmin;
	float y = ((ymax - ymin) * i) / (IMG_ROWS - 1) + ymin;


	//add focal length to make the converted image buffer coords a 3D point on the image plane
	float P1[4] = { x, y, focal, 1 };

	//Convert image plane coord into world coords storing in a temp variable
	//Temp var is in homogenious coords
	float temp[4] = { 0, 0, 0, 0 };
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			temp[i] += Mcw[i][j] * P1[j];
		}
	}

	//take converted point from temp and place into P1
	//also convert the coords in temp back into regular 3D coords from homogeneous

	for (int i = 0; i < 3; i++) {
		P1[i] = temp[i] / temp[3];
	}


	//Calculate vector from P0 (VPN) to P1 (point on image plane)
	for (int i = 0; i < 3; i++) {
		v0[i] = P1[i] - p0[i];
	}

	//get magnitude of new vector
	float vMag = magnitude(v0);

	//Normalize the new vector

	for (int i = 0; i < 3; i++) {
		v0[i] = v0[i] / vMag;
	}


}
/* rayBoxIntersection(float(&p0)[3], float(&v0)[3], float(&ts)[2]):
*  Input: camera origina p0[3], current ray v0[3], var to hold the two intersections with the box ts[2]
*  Output: n -> # of intersections
*  Description: extend the ray to check which sides of the 'box' the ray intersects with, if any.
*  Iterate through each of the six sides until all checks have been made.
*/
int rayBoxIntersection(float(&p0)[3], float(&v0)[3], float(&ts)[2]) {
	int n = 0;
	float t;
	float x, y, z;

	/* Side x = 0 */
	t = -p0[0] / v0[0];
	y = p0[1] + v0[1] * t;
	z = p0[2] + v0[2] * t;

	if ((y > 0 && y < ROWS) && (z > 0 && z < LAYS)) {
		ts[n] = t;
		n++;
	}

	/* Side x = 127 */
	t = (127 - p0[0]) / v0[0];
	y = p0[1] + v0[1] * t;
	z = p0[2] + v0[2] * t;

	if ((y > 0 && y < ROWS) && (z > 0 && z < LAYS)) {
		ts[n] = t;
		n++;
	}

	/* Side y = 0 */
	t = -p0[1] / v0[1];
	x = p0[0] + v0[0] * t;
	z = p0[2] + v0[2] * t;

	if ((x > 0 && x < COLS) && (z > 0 && z < LAYS)) {
		if (n >= 2) {

			n = 0;
		}

		ts[n] = t;
		n++;
	}

	/* Side y = 127 */
	t = (127 - p0[1]) / v0[1];
	x = p0[0] + v0[0] * t;
	z = p0[2] + v0[2] * t;

	if ((x > 0 && x < COLS) && (z > 0 && z < LAYS)) {
		if (n >= 2) {

			n = 0;
		}
		ts[n] = t;
		n++;
	}

	/* Side z = 0 */
	t = -p0[2] / v0[2];
	x = p0[0] + v0[0] * t;
	y = p0[1] + v0[1] * t;

	if ((x > 0 && x < COLS) && (y > 0 && y < ROWS)) {
		if (n >= 2) {

			n = 0;
		}
		ts[n] = t;
		n++;
	}

	/* Side z = 127 */
	t = (127 - p0[2]) / v0[2];
	x = p0[0] + v0[0] * t;
	y = p0[1] + v0[1] * t;

	if ((x > 0 && x < COLS) && (y > 0 && y < ROWS)) {
		if (n >= 2) {

			n = 0;
		}
		ts[n] = t;
		n++;
	}

	return n;
}

/* volumeRayTracing(float p0[3], float v[3], float ts[2]):
*  Input: camera origin p0[3], current ray v[3], two intersection points ts[2]
*  Output: c -> the colour value for the current point on the image plane
*  Description: step along the ray by a predetermined value and 
*  perform trilinear interpolation with the CT data and the Shading 
*  data to calculate their respective values at the current point along the ray.
*  then calculate the shading value for the current pixel and return.
*/
int volumeRayTracing(float p0[3], float v[3], float ts[2]) {
	float dt = 20.0; //size of steps taken along the ray
	int c = 0;	//accumulated color 
	float trans = 1.0; //accumulated transparency
	float ai, ci; //density and colour values for each point along the current ray
	float res[2]; //used to store the ai & ci values calculated within the function triLinearInterpolation()



	for (float t = ts[0]; t <= ts[1]; t += dt) {
		triLinearInterpolation(t, p0, v, res); //do the interpolation to get values for ai & ci
		ai = res[0]/250; //store the results
		ci = res[1]; 
		c += ci * ai * trans; //update the accumulated values
		
		trans *= (1.0 - ai);
		if (trans <= 0)
			return c;
	}

	return c;
}
/* triLinearInterpolation(float t, float p0[3], float v[3], float(&res)[2]): 
*  Input: current point on ray t, camera origin p0[3], current ray v[3], var to hold the results res[2]
*  Output: none
*  Description: performs trilinear interpolation
*/
void triLinearInterpolation(float t, float p0[3], float v[3], float(&res)[2]) {
	float x, y, z; //point choords
	float dx, dy, dz; //distance from current point to point on grid (dy -> yaxis, dz -> zaxis, dx -> xaxis)
	int i, j, k;	//location of closest point on grid to point on ray (i -> z, j -> y, k -> x)

	//calculate the current point on the ray using the parametric form where;
	//p0 -> vector origin, 
	//v -> unit length vector leaving p0, 
	//t -> the distance along the vector where the point is.
	x = p0[0] + (t * v[0]);
	y = p0[1] + (t * v[1]);
	z = p0[2] + (t * v[2]);

	/* finding the largest integer value for each x,y,z that is smaller than 
	   each of the values of the point (x,y,z) to locate the closest point within the CT array */
	i = (int)z;
	j = (int)y;
	k = (int)x;

	/* Find the distance between the point on the ray to the point found above by subtracting the x,y,z values
	   of the CT point from the point on the ray */
	dx = x - k;
	dy = y - j;
	dz = z - i;


	/* \/ Perform trilinear interpolation to get the color and transparency values for the current point along the ray \/ */

	//transparency value for ai
	res[0] = ((1 - dx) * (1 - dy) * (1 - dz) * CT[i][j][k]) + ((1 - dx) * (1 - dy) * (dz)*CT[i + 1][j][k]) +
		((dx) * (1 - dy) * (dz)*CT[i + 1][j][k + 1]) + ((dx) * (1 - dy) * (1 - dz) * CT[i][j][k + 1]) +
		((1 - dx) * (dy) * (1 - dz) * CT[i][j + 1][k]) + ((1 - dx) * (dy) * (dz)*CT[i + 1][j + 1][k]) +
		((dx) * (dy) * (dz)*CT[i + 1][j + 1][k + 1]) + ((dx) * (dy) * (1 - dz) * CT[i][j + 1][k + 1]);

	//colour value for ci
	res[1] = ((1 - dx) * (1 - dy) * (1 - dz) * SHADING[i][j][k]) + ((1 - dx) * (1 - dy) * (dz)*SHADING[i + 1][j][k]) +
		((dx) * (1 - dy) * (dz)*SHADING[i + 1][j][k + 1]) + ((dx) * (1 - dy) * (1 - dz) * SHADING[i][j][k + 1]) +
		((1 - dx) * (dy) * (1 - dz) * SHADING[i][j + 1][k]) + ((1 - dx) * (dy) * (dz)*SHADING[i + 1][j + 1][k]) +
		((dx) * (dy) * (dz)*SHADING[i + 1][j + 1][k + 1]) + ((dx) * (dy) * (1 - dz) * SHADING[i][j + 1][k + 1]);
	
}

/* magnitude(float v[3]): 
*  Input: vector v[3]
*  Output: vector's magnitude
*  Description: calculates the magnitude of the given vector
*/
float magnitude(float v[3]) {
	float res; //holds result of magnitude computation

	//compute the magnitude of the Vector
	res = sqrt(pow(v[0], 2) + pow(v[1], 2) + pow(v[2], 2));

	//retrun the computed magnitude
	return res;
}

/* dotProd(float v1[3], float v2[3]):
*  Input: vector v1[3], vector v2[3]
*  Output: dot product of the two vectors
*  Description: calculates the dot product of the two given vectors
*/
float dotProd(float v1[3], float v2[3]) {
	float res = 0; //holds the dot product result

	//compute the dot product of the two vectors
	for (int i = 0; i < 3; i++) {
		res += v1[i] * v2[i];
	}

	//return the computed dot product
	return res;
}

/* findnuv(float v1[], float v2[], float p[]):
*  Input: vector v1[], vector v2[], point p[]
*  Output: an array of the three orthogonal unit length vectors to v1 and v2
*  Description: calculate the 3 unit length vectors n,u,v used for creating the 
*  transformation matrices.
*/
vector<float> findnuv(float v1[], float v2[], float p[]) {
	float n[3]; //unit vector n
	float u[3]; //unit vector u
	float v[3]; //unit vector v
	vector<float> result;

	//------------- Calculate n ---------------


	float tM = magnitude(v1);

	for (int i = 0; i < 3; i++)
		n[i] = v1[i] / tM;
	//-------------------------------------------

	//------------- calculate u -----------------


	//cross product between vectors (P3-P1) and (P2-P1)
	float temp[3];
	temp[0] = v2[1] * v1[2] - v2[2] * v1[1];
	temp[1] = v2[2] * v1[0] - v2[0] * v1[2];
	temp[2] = v2[1] * v1[0] - v2[0] * v1[1];

	tM = magnitude(temp);

	for (int i = 0; i < 3; i++)
		u[i] = temp[i] / tM;
	//-----------------------------------------------

	//----------- Calculate v------------------------
	v[0] = n[1] * u[2] - n[2] * u[1];
	v[1] = n[2] * u[0] - n[0] * u[2];
	v[2] = n[0] * u[1] - n[1] * u[0];
	//-----------------------------------------------

	for (int i = 0; i < 3; i++)
		result.push_back(u[i]);

	for (int i = 0; i < 3; i++)
		result.push_back(v[i]);

	for (int i = 0; i < 3; i++)
		result.push_back(n[i]);


	return result;
}

/* findmWCmCW(float vpn[], float vup[], float vrp[]):
*  Input: camera view plane normal vpn[], camera up vector vup[], camera reference point vrp[]
*  Output: none
*  Description: calculate the two transformation world to camera and camera to world to be used 
*  for ray construction
*/
void findmWCmCW(float vpn[], float vup[], float vrp[]) {
	vector<float> nuv = findnuv(vpn, vup, vrp);// holds the u, v, n, values for the given vectors

	//-------------------- Camera to World -----------------------------------
	float invT[4][4] = {
						{1, 0, 0, vrp[0]},
						{0, 1, 0, vrp[1]},
						{0, 0, 1, vrp[2]},
						{0, 0, 0, 1}
	};

	float invR[4][4] = {
					  {nuv[0], nuv[3], nuv[6], 0},
					  {nuv[1], nuv[4], nuv[7], 0},
					  {nuv[2], nuv[5], nuv[8], 0},
					  {     0,      0,      0, 1}
	};



	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {

				Mcw[i][j] += invT[i][k] * invR[k][j];
			}
		}
	}

}
