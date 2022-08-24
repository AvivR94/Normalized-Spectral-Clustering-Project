#ifndef SPKMEANS_H
#define SPKMEANS_H
#include <math.h>
#include <stdio.h>

typedef struct eigens{
    int index;
    double value;
    double* vector;
}eigens;



FILE  *f;
struct eigens *eigensArray;

void wam(double** vectors_matrix, int n, int vec_length);
double** wam_calc(double** vectors_matrix, int n, int vec_length);
void ddg(double** vectors_matrix, int n, int vec_length);
double** ddg_calc(double** vectors_matrix, int n, int vec_length);
void lnorm(double** vectors_matrix, int n, int vec_length);
double** lnorm_calc(double** vectors_matrix, int n, int vec_length);
void eigengap_heuristic(double** vectors_matrix, int n, int vec_length);
void jacobi(double** vectors_matrix, int n, int vec_length);
struct eigens* jacobi_calc(double** vectors_matrix, int n);
double** matrix_Transpose(double** mat, int n);
void spk(double** vectors_matrix, int n, int vec_length, int k);
int comparator (const void* first, const void* second)
double** build_matrix_t_eigen(struct eigens* eigensArray, int n, int k);
double* copy_to_eigen_Values(double* vec_matrix, int n);
int get_vec_length(FILE *f);
int get_vec_num(FILE *f);
void matrix_multiplication(int rows_num, int columns_num, double** mat1, double** mat2, double** result);
double off_func(double** vectors_matrix[i][j], double** to_copy[i][j]);
double** jacobi_mat_for_print(struct eigens* eigensArray, int n, int vector_length);
void print_jacobi(struct eigens* eigensArray, int n, int vector_length);
void print_matrix(double** matrix, int n, int vector_length);
double **allocateMem(int n, int d);
void error_occurred();
int** largest_indexes(double** vectors_matrix, int n);
double* get_c_and_s(double** vectors_matrix, int lar_i, int lar_j);
void p_mat_maker(double** vectors_matrix, double** p_matrix, int lar_i, int lar_j);
void set_tocopy(double** vectors_matrix, double** to_copy);
int first_p_to_v(double** vectors_matrix, double** p_matrix, double** v_matrix,double** temp, int first);
void A_tag_calc(double** to_copy,double** vectors_matrix, int lar_i, int lar_j);
void A_to_A_tag(double** vectors_matrix,double** to_copy);
double *getFinalCentroids(double *centroids, double *elements, int k, int d, int n, int max_iter, double epsilone);
void resetCentroids(double* centroids,int k, int d);
void initClusters(int* elements_loc, int* items_number_clusters, int n, int k);
void saveCentroids(double* old_centroids, double* centroids, int k, int d);
void assignCentroids(double *ele,double *cntrds,int* in_clstrs,int* ele_loc,int k,int d,int n);
void updateCentroids(double* cntrds,double* ele,int* in_clstrs,int *ele_loc,int d,int n,int k);
int Convergence(double* old_centroids,double* centroids,int k,int d,int eps);
#endif 