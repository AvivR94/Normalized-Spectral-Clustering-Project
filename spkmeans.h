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
double** wam_calc(double** vectors_matrix, int n, int vec_length);
void ddg(double** vectors_matrix, int n, int vec_length);
double** ddg_calc(double** w_matrix, int n);
void lnorm(double** vectors_matrix, int n, int vec_length);
double** lnorm_calc(double** w_matrix, double** d_matrix, int n, int vec_length);
int eigengap_heuristic(struct eigens* eigensArray, int n);
void jacobi(double** vectors_matrix, int n, int vec_length);
struct eigens* jacobi_calc(double** vectors_matrix, int n, int vec_length);
double** matrix_Transpose(double** mat, int n, int vec_length);
void spk(double** vectors_matrix, int n, int vec_length, int k);
int comparator (const void* first, const void* second);
double** build_matrix_t_eigen(struct eigens* eigensArray, int n, int k);
double* copy_to_eigen_Values(double* vec_matrix, int n);
void matrix_multiplication(int rows_num, int columns_num, double** mat1, double** mat2, double** result);
double off_func(double** vectors_matrix, double** to_copy, int n);
double** jacobi_mat_for_print(struct eigens* eigensArray, int n, int vector_length);
void print_jacobi(struct eigens* eigensArray, int n, int vector_length);
void print_matrix(double** matrix, int n, int vector_length);
double** allocateMem(int n, int vector_length);
void error_occurred();
int* largest_indexes(double** vectors_matrix, int n);
double retrieveTheta(double** vectors_matrix, int lar_i, int lar_j);
double retrieveT(double theta);
double retrieveC(double t);
double retrieveS(double t, double c);
int sign(double theta);
void p_mat_maker(double** matrixA, double** p_matrix, int n, int lar_i, int lar_j);
void set_tocopy(double** vectors_matrix, double** to_copy, int n);
int v_calculation(double** p_matrix, double** v_matrix, int first,int vec_length, int n);
void A_tag_calc(double** to_copy,double** matrixA,int n, int lar_i, int lar_j);
void A_to_A_tag(double** vectors_matrix,double** to_copy, int n);
double** getFinalCentroids(double **centroids, double **elements, int k, int d, int n, int max_iter, double epsilone);
void resetCentroids(double** centroids,int k, int d);
void initClusters(int* elements_loc, int* items_number_clusters, int n, int k);
void saveCentroids(double** old_centeroids, double** centeroids, int k, int d);
void assignCentroids(double **ele,double **cntrds,int* in_clstrs,int* ele_loc,int k,int d,int n);
void updateCentroids(double** cntrds,double** ele,int* in_clstrs,int *ele_loc,int d,int n,int k);
int Convergence(double** old_centeroids,double** centeroids,int k,int d,int eps);
#endif 