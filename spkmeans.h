#ifndef SPKMEANS_H
#define SPKMEANS_H
#include <math.h>
#include <stdio.h>

typedef struct eigens{
    int index;
    double value;
    double* vector;
} eigens;

void wam(double** vectors_matrix, int n, int vec_length);
double** wamCalc(double** vectors_matrix, int n, int vec_length);
void ddg(double** vectors_matrix, int n, int vec_length);
double** ddgCalc(double** w_matrix, int n);
void lnorm(double** vectors_matrix, int n, int vec_length);
double** lnormCalc(double** w_matrix, double** d_matrix, int n);
void squareMatrixD(double** d_matrix, int n);
int eigengapHeuristic(struct eigens* eigensArray, int n);
void jacobi(double** vectors_matrix, int n);
struct eigens* jacobiCalc(double** A_matrix, int n, int print, int sort);
void copyMatrices(double** matrix_A, double** matrix_B, int n);
double** createMatrixI(int n);
void matrixTranspose(double** mat, int n);
int comparator (const void* first, const void* second);
double** createMatrixT(struct eigens* eigensArray, int n, int k);
double* copyToEigenVectors(double* vec_matrix, int n);
void matrixMult(int rows_num, int columns_num, double** mat1, double** mat2, double** result);
double offFunc(double** A_matrix, double** Atag_matrix, int n);
double** jacobiMatForPrint(struct eigens* eigensArray, int n);
void printJacobi(struct eigens* eigensArray, int n);
void printMatrix(double** matrix, int n, int vector_length);
double** allocateMem(int n, int vector_length);
void errorOccured();
int* retrieveLargestIndexes(double** A_matrix, int n);
double retrieveTheta(double** A_matrix, int lar_i, int lar_j);
double retrieveT(double theta);
double retrieveC(double t);
double retrieveS(double t, double c);
int sign(double theta);
double** createMatrixP(double** A_matrix,int n, int lar_i, int lar_j);
double** createMatrixAtag(double** A_matrix,int n,int lar_i,int lar_j);
void getFinalCentroids(double **centroids, double **elements, int k, int d, int n, int max_iter, double epsilone);
void resetCentroids(double** centroids,int k, int d);
void initClusters(int* elements_loc, int* items_number_clusters, int n, int k);
void saveCentroids(double** old_centeroids, double** centeroids, int k, int d);
void assignCentroids(double **ele,double **cntrds,int* in_clstrs,int* ele_loc,int k,int d,int n);
void updateCentroids(double** cntrds,double** ele,int* in_clstrs,int *ele_loc,int d,int n,int k);
int convergence(double** old_centeroids,double** centeroids,int k,int d,int eps);
#endif 