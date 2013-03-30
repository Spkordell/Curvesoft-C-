/** 
 * @file curves1.c
 * @author Steven Kordell <spkordell@wpi.edu>
 * @version 0.1
 * @section DESCRIPTION
 * Program to calculate the curvatures of a profile of points at various scales as well as perform several additional functions such as creating curvature distributions and conversion of data to lograrithmic scale for more convenient plotting.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>

//The available calculation methods
typedef enum {heronAndCalc, heron, calc, doubleDeriv, parabola} METHOD;

//A structure to hold a 2D data point
struct _point {
  double x;
  double y;
};

//type define the structure for convenience
typedef struct _point Point;


/** Converts a string to all upper case characters
 * @param str Pointer to the c-style string to convert
 * @return Pointer to the returned string
 * @author Steven Kordell
 */
char* convertToUpper(char *str){
    char *newstr, *p;
    p = newstr = strdup(str); //amke a cpy of the string
    while(*p++=toupper(*p)); //iterate over it, uppercasing the characters
    return newstr; //return the resulting string
}

/** Searches an array of strings for the existance of one particular string, not case sensititve
 * @param arr The array of strings to search
 * @param count Number of elements in the array
 * @param key The string to search the array for
 * @return The insex of the sting in the array or -1 if the string wasn't found
 * @author Steven Kordell
 */
int searchArrayForString(char* arr[], int count, char* key) {
  int n;
  char* elementUp;
  char* keyUp;
  for(n = 0; n < count; n++) {
    elementUp = convertToUpper(arr[n]);
    keyUp = convertToUpper(key);
    if (!strcmp(elementUp,keyUp)) {
      free(elementUp);
      free(keyUp);
      return n;
    }
    free(elementUp);
    free(keyUp);
   }
  return -1;
}

/** Reads an ascii text file containing 2D data points (seperated by spaces) to an array of Points
 * @param filename The path of the file to open and read
 * @param size Will be filled with the number of resulting array elements
 * @return A pointer to a dynamically allocated array of points containing the read data
 * @author Steven Kordell
*/
Point* readFileToArray(char* filename, int* size) {
  char temp[30];
  Point* points = NULL;
  *size = 0;

  FILE *infile = fopen(filename, "r");
  if (!infile) {
    printf("Could not open input file\n");
    return NULL;
  } else {
    while(fscanf(infile,"%s",temp) != EOF) {
      points = (Point*) realloc(points,sizeof(Point)*(*size+1));
      (points+(*size))->x=atof(temp);
      fscanf(infile,"%s",temp);
      (points+(*size))->y=atof(temp);
      (*size)++;
    }
    fclose(infile);
    return points;
  }
}

/** Prints an array of points to the terminal
 * @param points The array of points to print
 * @param size The number of array elements
 * @author Steven Kordell
*/
void printPoints(Point*  points, int size) {
  Point* p = points;
  while (size--) {
    if (isnan(p->x)) break;
    printf("(%.15f, %.15f)\n",p->x,p->y);
    p++;
  }
}

/** Writes an array of points to a file
 * @param filename The path of the file to write to
 * @param points The array of points to write to the file
 * @param size The number of elements in the array (or the number to write to the file)
 * @authro Steven Kordell
*/
int writePointsToFile(char* filename, Point* points, int size) {
  Point* p = points;
  FILE *outfile = fopen(filename, "w");
  if (!outfile) {
    printf("Could not open output file\n");
    return -1;
  } else {
    while (size--) {
      if (!isnan(p->x)) {
	fprintf(outfile,"%.15f,%.15f\r\n",p->x,p->y);
      }
      p++;
    }
  }
  fclose(outfile);
  return 0;
}

/** Sets all points values (x,y) in an array of points to NaN
 * @param points The array of points to set to NaN (will be modified in place)
 * @param The number of elements in the array 
 * @author Steven Kordell
 */
void initializePointsToNaN(Point* points, int count) {
  Point* p = points;
  while(count--) {
    p->x = 0.0 / 0.0;
    p->y = 0.0 / 0.0;
    p++;
  }
}


/** Squares a number
 * @param b The number to square (multiply by itself)
 * @return The result of the squaring operation
 */
double square(double b) {
    return  b*b;
}

/** Cubes a number (multiplies it by iteself three times)
 * @param b The number to cube
 * @return the Result of the cubing operation
 */
double cube(double b) {
  return b*b*b;
}

/* PositiveOrZero - If the input is negative, makes it zero, otherwise returns the input
 * @param b The incoming value to poz
 * @return The input value or, if a value is less than 0, returns 0.
 * @author Steven Kordell
*/
double poz(double b) {
  if (b < 0) {
    return 0;
  } else {
    return b;
  }
}

/** Calculates the curvature using herons method given three points
 * @param one The first data point
 * @param two The middle data point
 * @param three The last data point
 * @returns The curvature of the three points as determined using heron's method
 * @author Steven Kordell
*/ 
double heronCurvature(Point* one, Point* two, Point* three) {
  if (((one->y+three->y)/2) == two->y) {
    return 0; //If the three points lie on the same line, the curvature is 0;
  }

  double a = sqrt(square(two->x - one->x) + square(two->y - one->y));
  double b = sqrt(square(three->x - two->x) + square(three->y - two->y));
  double c = sqrt(square(three->x - one->x) + square(three->y - one->y));

  double s = .5 * (a + b + c);
  double k = sqrt(s*poz(s-a)*poz(s-b)*poz(s-c));
  double r = (a*b*c)/(4*k);
  
  if (((one->y+three->y)/2) < two->y) {
    r*=-1;
  }

  return 1/r;
}

/** Calculates the curvature using calculus  method given three points
 * @param one The first data point
 * @param two The middle data point
 * @param three The last data point
 * @returns The curvature of the three points as determined using the calculus method
 * @author Steven Kordell
*/ 
double calculusCurvature(Point* one, Point* two, Point* three) {
  double yprime = (three->y - one->y)/(2* (two->x-one->x)); 
  double ydoubleprime = (three->y-2*two->y+one->y)/square(two->x-one->x);
  return ydoubleprime/cube(sqrt(1+square(yprime)));
}

/** Calculates the curvature using herons method and the calculus method, then determines which is better, given three points
 * @param one The first data point
 * @param two The middle data point
 * @param three The last data point
 * @returns The curvature of the three points as determined using either heron's or the calculus method, whichever is better
 * @author Steven Kordell
*/ 
double chCurvature(Point* one, Point* two, Point* three, double scale) {
  double a = sqrt(square(two->x - one->x) + square(two->y - one->y));
  double b = sqrt(square(three->x - two->x) + square(three->y - two->y));

  double comp = (a*sqrt(2) + b*sqrt(2))/2;

  if (comp > three->x - one->x) {
    return calculusCurvature(one,two,three);
  } else {
    return heronCurvature(one,two,three);
  }
}

/** Calculates the curvature using the double derivative approximation method given three points
 * @param one The first data point
 * @param two The middle data point
 * @param three The last data point
 * @returns The curvature of the three points as determined using the double derivative approximation method
 * @author Steven Kordell
*/ 
double doubleDerivAproxCurvature(Point* one, Point* two,Point* three){
  return (three->y-2*two->y+one->y)/square(two->x-one->x);
}

/** Calculates the curvature using the parabola method given three points
 * @param one The first data point
 * @param two The middle data point
 * @param three The last data point
 * @returns The curvature of the three points as determined using the parabola method
 * @author Steven Kordell
*/ 
double parabolaCurvature(Point* one, Point* two,Point* three){

  double x1 = one->x;
  double x2 = two->x;
  double x3 = three->x;
  double y1 = one->y;
  double y2 = two->y;
  double y3 = three->y;

  //calculate parabola coefficients A,B, and C of form Ax^2+Bx+C;
  double denom = (x1 - x2)*(x1 - x3)*(x2 - x3);
  double A = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
  double B = (square(x3) * (y1 - y2) + square(x2) * (y3 - y1) + square(x1) * (y2 - y3)) / denom;
  //  double C = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;

  double yprime = 2*A*(two->x)+B;
  double ydoubleprime = 2*A;

  return ydoubleprime/sqrt(cube(1+square(yprime)));
}

/** Determines the sampling interval of a profile
 * @param profile A point array representing the profile (must contain at least two points)
 * @return The sampling interval of the profile
 * @author Steven Kordell
 */
double samplingInterval(Point* profile) {
  return ((profile+1)->x) - (profile->x);
}

/** Determines the curvature of a profile at a specific scale using the indicated method
 * @param profile The profile to calculate the curvatures of, as an array of Points
 * @param size The number of elements in the array
 * @param scale The scale to calculate the curvature at
 * @param method The method to use in calculating the curvature
 * @return An array of points containing the positions curvatures were calculated at and the corrosponding curvature 
 * @author Steven Kordell
 */
Point* curvatureAtScale(Point* profile, int size, double scale, METHOD method) {
  Point* one = profile;
  Point* two = profile + (int)(scale/samplingInterval(profile))/2;
  Point* three = profile +  (int)(scale/samplingInterval(profile));

  Point* curves = (Point*) malloc(sizeof(Point) * size);
  initializePointsToNaN(curves,size); //initialize the array
  Point* p = curves;

  while (three < (profile+size)) { //was -1
    p->x = two->x;
    switch(method) {
    case heronAndCalc:
      (p++)->y = chCurvature(one++,two++,three++,scale);
      break;
    case heron:
      (p++)->y = heronCurvature(one++,two++,three++);
      break;
    case calc:
      (p++)->y = calculusCurvature(one++,two++,three++);
      break;
    case doubleDeriv:
      (p++)->y = doubleDerivAproxCurvature(one++,two++,three++);
      break;      
    case parabola:
      (p++)->y = parabolaCurvature(one++,two++,three++);
      break;    
    }  
  }
  return curves;
}

/** Determines the length of a profile 
 * @param profile The array of points describing the profile
 * @param size The number of points in the profile
 * @return The length of the profile
 * @author Steven Kordell
 */
double profileLength(Point* profile, int size) {
  return profile[size-1].x;
}

/** Calculates all the possible scales curvature can be calculated at and returns an array containing them
 * @param profile The profile to determing the available scales of, as an array of points
 * @param size The number of elements in the profile
 * @param scaleCount Will be filled with the number of available scales/the size of the resulting array
 * @return A dynamically allocated array of doubles containing the possible scales
 * @author Steven Kordell
*/
double* allPossibleScales(Point* profile, int size, int* scaleCount) {
  double* scalesCalculated = NULL;
  *scaleCount = 0;
  double scale = samplingInterval(profile) * 2;

  while (scale <= profile[size-1].x) {
    (*scaleCount)++;
    scalesCalculated = (double*) realloc(scalesCalculated,sizeof(double) * (*scaleCount));
    scalesCalculated[(*scaleCount)-1] = scale;
    scale += 2*samplingInterval(profile);
  }

  return scalesCalculated;
}

/** writes all the points referenced to in an array of pointers to points to a file with their corrosponding scale (for use with an array of curvatures calculated at various scales)
 * @param filename The path of the file to write to
 * @param ptrToPoints The Multideimmensional array of pointers to points (an array of curvatures at multiple scales)
 * @param ptrArraySize The size of the array of pointers to arrays of points
 * @param targetArraySize The size of an array pointed to by the array of pointers to arrays of points
 * @param scalesCalculated An array of doubles containing the scales the curvatures were calculated at
 * @return -1 if the operation failed, 0 if the operation was successful
 * @author Steven Kordell
*/
int writeArrayOfPointersToPointsToFile(char* filename, Point** ptrToPoints, int ptrArraySize, int targetArraySize, double* scalesCalculated) {
  int i;
  int j;
  Point* p;
  FILE *outfile = fopen(filename, "w");
  if (!outfile) {
    printf("Could not open output  file\n");
    return -1; //operation failed
  }
  fprintf(outfile,"Position, Scale, Curvature\r\n");
  for(i = 0; i < ptrArraySize; i++) {
    p = ptrToPoints[i];
    j = targetArraySize;
    while (j--) {
      if (!isnan(p->x)) {
	fprintf(outfile,"%.15f,%.15f,%.15f\r\n",p->x,scalesCalculated[i],p->y);
      }
      p++;
    }
  }
  //printf("File Size: %d\n",ftell(outfile));
  if (ftell(outfile) < 2) { //If the file is empty
    printf("Disk is full. Could not write to file");
    return -1;
  }

  fclose(outfile);
  return 0; //operation succeeded
}

/** Frees the memory used by a multidimensional array of pointers to points
 * @param ptrArray The array to free
 * @param the number of elements in the array
 */
void freeArrayOfPointersToPoints(Point** ptrArray, int count) {
  int i;
  for (i = 0; i < count; i++) {
    free(ptrArray[i]);
  }
  free(ptrArray);
}

/** Determines if there is a hyphan found before a comma in an incoming string
 * @param s The string to search
 * @return 1 if a hyphen is found before a comma a comma is found, otherwise  0
 * @author Steven Kordell
*/
char isHyphenBeforeComma(char* s) {
  char* p = s;
  while (*p && *p != ',') {
    if (*p == '-') {
      return 1;
    }
    p++;
  }
  return 0;
}

/** Determines if there is colon found before a comma in an incoming string
 * @param s The string to search
 * @return 1 if a colon  is found before a comma a comma is found, otherwise  0
 * @author Steven Kordell
 */
char isColonBeforeComma(char* s) {
  char* p = s;
  while (*p && *p != ',') {
    if (*p == ':') {
      return 1;
    }
    p++;
  }
  return 0;
}

/** Takes a scale and the sampling interval and determines if the scale is valid
 * @param scale The scale to determine wether or not it is a valid scale
 * @param samplingInterval The sampling interval of the profile
 * @return 1 if the scale is valid, 0 otherwise
 * @author Steven Kordell
*/
char isScaleValid(double scale, double samplingInterval, double maxScale) {
  if (!fmod(scale, 2*samplingInterval) && scale <= maxScale && scale != 0) {
    return 1;
  } else {
    return 0;
  }
}

/** Takes a user defined string of scales to calculate and extracts an array of scales to calculate
 * @param param The incoming parameter string
 * @param scaleCount Will be filled with the number of scales to calculate
 * @param samplingInterval The sampling interval of the profile
 * @return An array of doubles containg the scales the user would like curvature calculated at
 * @author Steven Kordell
 */
double* determineScalesToCalculate(char* param, int* scaleCount, double samplingInterval, double maxScale) {
  double buffer;
  double min;
  double max;
  double increment;
  char reverse; //0 if max is greater than min, 1 otherwise
  char wasRanOnce = 0;
  double* scalesToCalculate = (double*) malloc(sizeof(double));
  *scaleCount = 0;

  char* p;
  p = param;
  while (*(p++)) {
    if (*p == ',' || wasRanOnce == 0) {
      wasRanOnce = 1;
      if (*scaleCount == 0) {
	p-=2;
      }
      p++;
      if(isHyphenBeforeComma(p)){
	if (isColonBeforeComma(p)) {
	  sscanf(p,"%lf-%lf:%lf",&min,&max,&increment);
	} else {
	  increment = 2*samplingInterval;
	  sscanf(p,"%lf-%lf",&min,&max);
	}
	if (max > min) {
	  while (min <= max+samplingInterval) {
	    if (isScaleValid(min,samplingInterval,maxScale)) {
	      (*scaleCount)++;
	      scalesToCalculate = (double*) realloc(scalesToCalculate,sizeof(double)*(*scaleCount));  
	      scalesToCalculate[(*scaleCount)-1] = min;
	    }
	    min += increment;
	  }
	} else {
	  while (min+samplingInterval >=  max) {
	    if (isScaleValid(min,samplingInterval,maxScale)) {
	      (*scaleCount)++;
	      scalesToCalculate = (double*) realloc(scalesToCalculate,sizeof(double)*(*scaleCount));  
	      scalesToCalculate[(*scaleCount)-1] = min;
	    }
	    min -= increment;
	  }	  
	}
      } else {
	sscanf(p,"%lf,",&buffer);
	if (isScaleValid(buffer,samplingInterval,maxScale)) {
	  (*scaleCount)++;
	  scalesToCalculate = (double*) realloc(scalesToCalculate,sizeof(double)*(*scaleCount));
	  scalesToCalculate[(*scaleCount)-1] = buffer;
	}
      }
    }
  } 
  return scalesToCalculate;
}

/** Calculates the curvature of a profile at multiple scales
 * @param profile The profile (an array of Points) to calcuate the curvatures of
 * @param pointCount The number of points in the profile
 * @param scalesToCalculate An array containing the scales to calculate the curvature at
 * @param scaleCount The number of scales to calculate the curvature at
 * @param method The method to use in calculating the curvatures
 * @return An array of pointers to points (a multidimensional array) containing the calculated curvatures at each scale
 * @author Steven Kordell
 */
Point** curvatureAtScales(Point* profile, int pointCount, double* scalesToCalculate, int scaleCount, METHOD method) {
  Point** arrayOfCurvatures = (Point**) malloc(sizeof(Point*) * scaleCount);
  int i;  
  for(i = 0; i < scaleCount; i++) {
    arrayOfCurvatures[i] = curvatureAtScale(profile, pointCount, scalesToCalculate[i], method);  
  }
  return arrayOfCurvatures;
}

/** Converts the curvature to a log scale axis
 * @param curvatureAtMultipleScales An multidimmensional array of curvatures and scales to convert (in place)
 * @param pointCount The number of points in each scale
 * @param scaleCount The number of scales used to build the array
 * @author Steven Kordell
*/
void makeCurvatureLog(Point** curvatureAtMultipleScales,int pointCount, int scaleCount) {
  int i;
  int j;
  for(i = 0; i < scaleCount; i++) {
    for(j = 0; j < pointCount; j++) {
      if (curvatureAtMultipleScales[i][j].y > 0){
	curvatureAtMultipleScales[i][j].y = log(curvatureAtMultipleScales[i][j].y+1);
      } else {
	curvatureAtMultipleScales[i][j].y = -1 * log(fabs(curvatureAtMultipleScales[i][j].y)+1);
      }
    }
  }
}

/** Converts the postion to a log scale axis
 * @param curvatureAtMultipleScales An multidimmensional array of curvatures and scales to convert (in place)
 * @param pointCount The number of points in each scale
 * @param scaleCount The number of scales used to build the array
 * @author Steven Kordell
*/
void makePositionLog(Point** curvatureAtMultipleScales,int pointCount, int scaleCount) {
  int i;
  int j;
  for(i = 0; i < scaleCount; i++) {
    for(j = 0; j < pointCount; j++) {
      if (curvatureAtMultipleScales[i][j].x > 0){
	curvatureAtMultipleScales[i][j].x = log(curvatureAtMultipleScales[i][j].x+1);
      } else {
	curvatureAtMultipleScales[i][j].x = -1 * log(fabs(curvatureAtMultipleScales[i][j].x)+1);
      }
    }
  }
}

/** Converts the scale to a log scale axis
 * @param curvatureAtMultipleScales An multidimmensional array of curvatures and scales to convert (in place)
 * @param pointCount The number of points in each scale
 * @param scaleCount The number of scales used to build the array
 * @author Steven Kordell
*/
void makeScaleLog(double* calculatedScales, int scaleCount) {
  int i;
  for(i = 0; i < scaleCount; i++) {
      if (calculatedScales > 0){
	calculatedScales[i] = log(calculatedScales[i]+1);
      } else {
	calculatedScales[i] = -1 * log(fabs(calculatedScales[i])+1);
      }
    
  }
}

/**Determines the maximum curvature in an array of curvatures at various scales
 * @param curvaturAtMultipleScales The multidimensional array of curvatures at various scales
 * @param pointCount The number of points in each scale
 * @param scaleCount The number of scales used to build the array
 * @return The maximum curvature
 * @author Steven Kordell
*/
double getMaxCurvature(Point** curvatureAtMultipleScales,int pointCount, int scaleCount) {
  int i;
  int j;
  double maxCurvature = -DBL_MAX;
  for(i = 0; i < scaleCount; i++) {
    for(j = 0; j < pointCount; j++) {
      if (curvatureAtMultipleScales[i][j].y > maxCurvature && curvatureAtMultipleScales[i][j].y < DBL_MAX){
	maxCurvature = curvatureAtMultipleScales[i][j].y;
      }
    }
  }
  return maxCurvature;
}

/**Determines the minimum curvature in an array of curvatures at various scales
 * @param curvaturAtMultipleScales The multidimensional array of curvatures at various scales
 * @param pointCount The number of points in each scale
 * @param scaleCount The number of scales used to build the array
 * @return The minimum curvature
 * @author Steven Kordell
*/
double getMinCurvature(Point** curvatureAtMultipleScales,int pointCount, int scaleCount) {
  int i;
  int j;
  double minCurvature = DBL_MAX;
  for(i = 0; i < scaleCount; i++) {
    for(j = 0; j < pointCount; j++) {
      if (curvatureAtMultipleScales[i][j].y < minCurvature && curvatureAtMultipleScales[i][j].y > -DBL_MAX){
	minCurvature = curvatureAtMultipleScales[i][j].y;
      }
    }
  }
  return minCurvature;
}

/** creates a distribution of curvatues
 * @param curvatures The multidimensional array of curvatures at various scales
 * @param pointCount The number of points in each scale
 * @param scaleCount The number of scales used to build the array
 * @param param The sting parameter passed to the function by the user (used to determing how to set up the bins)
 * @param binCount Will be filled with the number of bins needed to store the distribution
 * @retrun An array of points, where x is the low value of the bin (inclusive) and  y is the frequency of that curvature
 * @author Steven Kordell
 */
Point* makeCurvatureDistribution(Point** curvatures, int pointCount,int scaleCount, char* param, int* binCount) {
  double binSize;
  double bin;
  double maxBin = getMaxCurvature(curvatures,pointCount,scaleCount);
  double minBin = getMinCurvature(curvatures,pointCount,scaleCount);
  Point* distribution;
  char keepZeros; //if 0, will turn all 0's to NaN's

  if (strchr(param,'z') != NULL) {
    sscanf(param,"z%lf",&binSize);
    keepZeros = 1;
  } else {
    binSize = atof(param);
    keepZeros = 0;
  }

  maxBin = ceil(maxBin/ binSize) * binSize;
  minBin = floor(minBin/ binSize) * binSize;

  *binCount = 0;
  for (bin = minBin; bin < maxBin; bin+=binSize) {
    (*binCount)++;
  } 
  distribution = (Point*) malloc(sizeof(Point) * (*binCount));
  
  //set up the bins
  Point* k = distribution; //temporaty pointer to hold the distribution data
  bin = minBin;
  while (k != distribution+(*binCount)) {
    k->x = bin;
    k->y = 0;
    bin += binSize;
    k++;
  }

  Point** i = curvatures; //temporaty pointer to the array of pointers to arrays of points
  Point* j; //temporary pointer to the array of points
  while (i != curvatures+scaleCount) {
    j = *i;
    while (j != (*i)+pointCount) {
      if (!isnan(j->y)) {
	(distribution[(int)(((floor((j->y) / binSize) * binSize)-minBin)/binSize)].y)++;
      }
      j++;
    }
    i++;
  }

  //The print function being used won't print points with and x-value of NaN to a file, so remove bins with 0 items to make the file more compact by setting the corrosponging x point to NaN;
  if (!keepZeros) {
    int z;
    for (z = 0; z < *binCount; z++){
      if (distribution[z].y == 0) {
	distribution[z].x = 0.0/0.0;
      }
    }
  }

  return distribution;
}

/** Prints information about how to use the program to the terminal to aid the user
 *@author Steve Kordell
 */
void printHelp() {
    printf("--------------------INFORMATION--------------------\n");
    printf("Created by Steven Kordell at the WPI Surface Metrology Lab: March 2013\n\n");

    printf("Load a profile:\n");
    printf("Enter \"-p\" followed by the path of the profile to load.\n\n");

    printf("Select a calculation method:\n");
    printf("Enter \"-m\" followed by a method idtentifier. c = calculus method, h = herons method, ch = hybrid calculus and herons, d = double derivitive approximation, p = parabola\n\n");    

    printf("Selecting scales to calculate:\n");
    printf("Enter \"-s\" followed by any of the folloing.\n");
    printf("For a single scale, enter the scale. e.g. \"100\".\n");
    printf("For all possible scales, enter \"all\".\n");
    printf("For all possibles scales between two scales (inclusive), enter scale range with hypens, e.g. \"100-150\".\n");
    printf("Use a colon to denote an increment. e.g. \"100-1	exit(22);50:25\" will produce scales of 100,125,and 150.\n");
    printf("Enter ranges backwards to produce reuslt in decrementing order. e.g. \"150-100:25\".\n");
    printf("Enter multiple scales or scale ranges with commas. e.g. \"100,105-110,115,118\" \n\n");

    printf("Optional: Create distribution of curvatures for use in a histogram\n");
    printf("Enter \"-d\" followed by the width of each bin. e.g. \"-d 20\" to create a distribution of size 20. Will output the data to the output file in the form of bin_start,frequency. \n");
    printf("By defualt, empty bins will be removed from the output file. To disable this, precede the bin size with a lowercase z. e.g. \"-d z20\". \n\n");

    printf("Optional: Output results on logarithmic scale:\n");
    printf("Note: Won't work on for curvature of a single scale and won't effect output if a curvature distribution is used\n");
    printf("Enter \"-l\" follwed by a list of items to make logarithmic. c = curvature, p = position, s = scale. e.g. \"-l cs\" will output the curvature and scale  scaled logarithmically.\n\n");

    printf("Saving the results:\n");
    printf("Enter \"-o\" followed by the path of the file to save the results in. File extensions are ignored, but using a *.csv extension could be more convenient for later use. \n\n");
}



int plotProfile(char* filename) {
  FILE *outfile = fopen("temp_profile_plot.yodawg", "w");
  if (!outfile) {
    printf("Could not open temporay plot file\n");
    return -1;
  } else {
    fprintf(outfile,"plot \"%s\"\npause -1\n",filename);
  }
  fclose(outfile);
  system("gnuplot temp_profile_plot.yodawg");
  system("rm temp_profile_plot.yodawg");
  return 1;
}

int plot3D(Point** ptrToPoints, int ptrArraySize, int targetArraySize, double* scalesCalculated) {
  char* filename = "temp_3d_plot_data.yodawg";
  int i;
  int j;
  Point* p;
  FILE* outfile = fopen(filename, "w");
  if (!outfile) {
    printf("Could not open output  file\n");
    return -1; //operation failed
  }
  for(i = 0; i < ptrArraySize; i++) {
    p = ptrToPoints[i];
    j = targetArraySize;
    while (j--) {
      if (!isnan(p->x)) {
	fprintf(outfile,"%.15f %.15f %.15f\n",p->x,scalesCalculated[i],p->y);
      }
      p++;
    }
  }
  if (ftell(outfile) < 2) { //If the file is empty
    printf("Disk is full. Could not write to file");
    return -1;
  }
  fclose(outfile);

  outfile = fopen("temp_3d_plot.yodawg", "w");
  if (!outfile) {
    printf("Could not open temporay plot file\n");
    return -1;
  } else {
    fprintf(outfile,"splot \"%s\"\npause -1\n",filename);
  }
  fclose(outfile);
  system("gnuplot temp_3d_plot.yodawg");
  system("rm temp_3d_plot_data.yodawg");
  system("rm temp_3d_plot.yodawg");

  return 0; //operation succeeded
}

int plotDistribution(Point* points, int size) {
  char* filename = "temp_dist_plot_data.yodawg";

  Point* p = points;
  FILE *outfile = fopen(filename, "w");
  if (!outfile) {
    printf("Could not open output file\n");
    return -1;
  } else {
    while (size--) {
      if (!isnan(p->x)) {
	fprintf(outfile,"%.15f %.15f\n",p->x,p->y);
      }
      p++;
    }
  }
  fclose(outfile);


  outfile = fopen("temp_dist_plot.yodawg", "w");
  if (!outfile) {
    printf("Could not open temporay plot file\n");
    return -1;
  } else {
    fprintf(outfile,"plot \"%s\" \npause -1\n",filename);
  }
  fclose(outfile);
  system("gnuplot temp_dist_plot.yodawg");
  system("rm temp_dist_plot_data.yodawg");
  system("rm temp_dist_plot.yodawg");

  return 0; //operation succeeded
}


int main (int argc, char *argv[]) {
  //profile information
  int pointCount; //The number of points in the input profile
  Point* profile; //An array of Points to store the profile

  //single scale
  // double scale = 0; //A single scale to calculate the curvature at
  // Point* curvatureAtSingleScale; //The resulting curvatures

  //multiple scales
  int scaleCount = 0;//The number of scales to calculate
  Point** curvatureAtMultipleScales; //An array of an array of resulting curvatur
  double* scalesToCalculate = NULL; //An array of the scales that were calculated
  Point* curvatureDistribution = NULL; //An array of curvatures and the frequency that curvature occurs at
  int curvatureDistributionBinCount; //Number of bins in the curvature distribution

  //miscellaneous variables
  int paramIndex; //For processing cammand line arguments
  METHOD method = heronAndCalc; //Stores the calculation method. Defaulted to heron and calc.

  printf("\n");

  //Load the profile
  if ((paramIndex = searchArrayForString(argv,argc,"-p")) != -1) {
    printf("Loading profile...\t\t\t\t\t\t");
    fflush(stdout);
    if((profile = readFileToArray(argv[paramIndex+1],&pointCount)) != NULL) {
      printf("DONE! - Sampling Interval: %f\n",samplingInterval(profile));
      // plotProfile(argv[paramIndex+1]);
    } else {
      printf("FAILED!\n");
      return 1;
    }
  } else {
    printf("usage: %s -p filename [options]\n\n", argv[0]);
    printHelp();
    return 1;
  }

  //Determine method to use
  if ((paramIndex = searchArrayForString(argv,argc,"-m")) != -1) {
    if (!strcmp(argv[paramIndex+1],"h")) {
      method = heron;
    }
    if (!strcmp(argv[paramIndex+1],"c")) {
      method = calc;
    } 
    if (!strcmp(argv[paramIndex+1],"ch") || !strcmp(argv[paramIndex+1],"hc")) {
      method = heronAndCalc;
    }  
    if (!strcmp(argv[paramIndex+1],"d")) {
      method = doubleDeriv;
    }  
    if (!strcmp(argv[paramIndex+1],"p")) {
      method = parabola;
    }  
  }

  //Read in the scale and calculate the curvatures
  if ((paramIndex = searchArrayForString(argv,argc,"-s")) != -1) {
    printf("Calculating curvatures...\t\t\t\t\t");
    fflush(stdout);
    if (!strcmp(argv[paramIndex+1],"all")) {
      scalesToCalculate = allPossibleScales(profile,pointCount, &scaleCount);
      curvatureAtMultipleScales = curvatureAtScales(profile, pointCount, scalesToCalculate, scaleCount, method);
      printf("DONE! - %d scales calculated.\n",scaleCount);
    } else {
      scalesToCalculate = determineScalesToCalculate(argv[paramIndex+1], &scaleCount,samplingInterval(profile),profileLength(profile,pointCount));
      if (scaleCount > 0) {
	curvatureAtMultipleScales = curvatureAtScales(profile, pointCount, scalesToCalculate, scaleCount, method);
	printf("DONE! - %d scales calculated.\n",scaleCount);
      } else {
	printf("NO VALID SCALES ENTERED\n");
      }
    }
  }

  //print the mininum and maximum curvature
  if (scaleCount > 0) {
    printf("Minimum Curvature: %f, Maximum Curvature: %f\n",getMinCurvature(curvatureAtMultipleScales,pointCount,scaleCount),getMaxCurvature(curvatureAtMultipleScales,pointCount,scaleCount));
  }

  //Make a curvature distribution (histogram)
  if (scaleCount > 0) {
    if ((paramIndex = searchArrayForString(argv,argc,"-d")) != -1) {
      printf("Calculating curvature vs. frequency distribution...\t\t");
      fflush(stdout);
      curvatureDistribution = makeCurvatureDistribution(curvatureAtMultipleScales,pointCount,scaleCount,argv[paramIndex+1],&curvatureDistributionBinCount);
      plotDistribution(curvatureDistribution,curvatureDistributionBinCount);
      printf("DONE! - %d bins\n",curvatureDistributionBinCount);
    }
  }

  //convert output to log scale if desired  
  if (scaleCount > 0) {
    if ((paramIndex = searchArrayForString(argv,argc,"-l")) != -1) {
      printf("Converting to log scale...\t\t\t\t\t"); 
      fflush(stdout);
      if(strchr(argv[paramIndex+1],'c') != NULL) {
	makeCurvatureLog(curvatureAtMultipleScales,pointCount,scaleCount);
      }
      if(strchr(argv[paramIndex+1],'p') != NULL) {
	makePositionLog(curvatureAtMultipleScales,pointCount,scaleCount);
      }
      if(strchr(argv[paramIndex+1],'s') != NULL) {
	makeScaleLog(scalesToCalculate,scaleCount);
      }
      printf("DONE!\n");
    }
  }

  //save the results
  if ((paramIndex = searchArrayForString(argv,argc,"-o")) != -1) {
    if (scaleCount > 0) {
      int result;
      printf("Saving results  to file:\n\t %-55s",argv[paramIndex+1]);
      fflush(stdout);
      if (curvatureDistribution == NULL) {
	if (scalesToCalculate != NULL) {
	  result = writeArrayOfPointersToPointsToFile(argv[paramIndex+1],curvatureAtMultipleScales,scaleCount,pointCount,scalesToCalculate);
	  //plot3D(curvatureAtMultipleScales,scaleCount,pointCount,scalesToCalculate);
	} else {
	  if (scaleCount == 1) {//if (scale) {
	    result = writePointsToFile(argv[paramIndex+1],*curvatureAtMultipleScales,pointCount);
	  } else {
	    result = -1;
	  }
	}
      } else {
	result = writePointsToFile(argv[paramIndex+1],curvatureDistribution,curvatureDistributionBinCount);
      }
      if (result != -1) {
	printf("DONE!\n");
      } else {
	printf("FAILED!\n");
      }
    } else {
      printf("No Results computed. Please specify a scale with -s or specify all to compute curvature at all scales. Scales must be evenly divisable by the sampling interval\n");
    }
  }

  printf("\n");

  //free any memory used
  free(profile);
  if (scaleCount > 0) {
    freeArrayOfPointersToPoints(curvatureAtMultipleScales,scaleCount);
    free(scalesToCalculate);
  }
  return 0;
}
