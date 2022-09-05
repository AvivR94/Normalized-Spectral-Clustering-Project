#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "spkmeans.h"
#define EPSILON 0.00001

int main(int argc, char* argv[]){
    int n;
    int i;
    int j;
    int r;
    int vec_length;
    char* goal;
    char* input_filename;
    double** vectors_matrix;
    FILE* input_file;
    double vector_item;
    char comma;
    char tav;
    if (argc != 3){
        printf ("Invalid Input!");
        exit(1);
    }

    goal = argv[1]; /* the enum type */
    input_filename = argv[2]; /* the input */

    input_file = fopen(input_filename, "r");
    if (input_file == NULL)
        errorOccured();

    /* counts the dimension of the vectors(vec_length), and the number of vectors (n) */
    vec_length = 1;
    n = 0;
    tav = fgetc(input_file);
    while (tav != '\n'){
        if (tav == ',')
            vec_length++;
        tav = fgetc(input_file);
    }
    while (tav != EOF){
        if (tav == '\n')
            n++;
        tav = fgetc(input_file);
    }
    fseek(input_file,0,0); /* rewind file */

    /* initiallization of vectors matrix */
    vectors_matrix = allocateMem(n, vec_length);
    if (vectors_matrix == NULL)
        errorOccured();

    /* insert vectors to vectors_matrix */
    for (i = 0; i < n; i++){
        for (j = 0; j < vec_length; j++){
            fscanf(input_file, "%lf%c", &vector_item, &comma);
            vectors_matrix[i][j] = vector_item;
        }
    }
    r = fclose(input_file);
    if (r!=0){
        errorOccured();
    }
    /* check which goal to choose */
        if(strcmp(goal, "wam") == 0){
            wam(vectors_matrix, n, vec_length);
        }
        else if (strcmp(goal, "ddg") == 0){
            ddg(vectors_matrix, n, vec_length);
        }
        else if (strcmp(goal, "lnorm") == 0){
            lnorm(vectors_matrix, n, vec_length);
        }
        else if (strcmp(goal, "jacobi") == 0){
            jacobi(vectors_matrix, n);
        }
        
    for(i=0; i<n; i++)
        free(vectors_matrix[i]);
    free(vectors_matrix);

    return 0;
}

void wam(double** vectors_matrix, int n, int vec_length){
    int i;
    double** w_matrix;
    printf("Start wam\n");
    w_matrix = wamCalc(vectors_matrix, n, vec_length);
    printf("gonna print\n");
    printMatrix(w_matrix, n, n);
    printf("after print\n");
    for(i = 0; i < n; i++)
        free(w_matrix[i]);
    free(w_matrix);
}

double** wamCalc(double** vectors_matrix, int n, int vec_length){
    int i;
    int j;
    int s;
    double sum;
    double** w_matrix;

    w_matrix = allocateMem(n, n);
    if (w_matrix == NULL)
        errorOccured();  
    printf("Start wamCalc\n");
    /* calculating the values in each cell */
    for (i = 0; i < n; i++)
        w_matrix[i][i] = 0; /* the diagonal line in matrix's values are 0 */
    printf("print vectors matrix:\n");
    printMatrix(vectors_matrix, n, vec_length);
    for (i = 0; i < n; i++){
        for (j= i + 1; j < n; j++){
            sum = 0;
            for (s = 0; s < vec_length; s++)
                sum += pow((vectors_matrix[i][s] - vectors_matrix[j][s]),2);
            w_matrix[i][j] = exp((-1)*(sqrt(sum)/2));
            w_matrix[j][i] = exp((-1)*(sqrt(sum)/2));
        }
    }
    return w_matrix;
}

void ddg(double** vectors_matrix, int n, int vec_length){
    double** w_matrix;
    double** d_matrix;
    int i;
    w_matrix = wamCalc(vectors_matrix, n, vec_length);
    d_matrix = ddgCalc(w_matrix, n);

    printMatrix(d_matrix, n, n);

    for(i = 0; i < n ; i++){
        free(w_matrix[i]);
        free(d_matrix[i]);
    }
    free(w_matrix);
    free(d_matrix);
}

double** ddgCalc(double** w_matrix, int n){
    double** d_matrix;
    int i;
    int j;

    d_matrix = allocateMem(n, n);
    if (d_matrix == NULL)
        errorOccured();
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            d_matrix[i][j] = 0;

    /* calculate the sum of each row and place in the diagonal line of d_matrix */
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++)
            d_matrix[i][i] += w_matrix[i][j];
    }

    return d_matrix;
}
        
void lnorm(double** vectors_matrix, int n, int vec_length){
    int i;
    double** d_matrix;
    double** w_matrix;
    double** lnorm_matrix;
    w_matrix = wamCalc(vectors_matrix, n, vec_length);
    d_matrix = ddgCalc(w_matrix, n);
    lnorm_matrix = lnormCalc(w_matrix, d_matrix, n, vec_length);

    printMatrix(lnorm_matrix, n, n);

    for(i = 0; i < n; i++){
        free(w_matrix[i]);
        free(d_matrix[i]);
        free(lnorm_matrix[i]);
    }
    free(w_matrix);
    free(d_matrix);
    free(lnorm_matrix);
}

double** lnormCalc(double** w_matrix, double** d_matrix, int n, int vec_length){
    int i;
    int j;
    double** laplacian_matrix;
    double** union_matrix;
    double** opp_d_matrix;
    double** result_matrix;
    double** temp_matrix;

    /* Laplacian_matrix initiallization */
    laplacian_matrix = allocateMem(n, n);
    if (laplacian_matrix == NULL)
        errorOccured();
    union_matrix = allocateMem(n, n);
    if (union_matrix == NULL)
        errorOccured();
    opp_d_matrix = allocateMem(n, n);
    if (opp_d_matrix == NULL)
        errorOccured();
    result_matrix = allocateMem(n, n);
    if (result_matrix == NULL)
        errorOccured();
    temp_matrix = allocateMem(n, n);
    if (temp_matrix == NULL)
        errorOccured();

    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            /*laplacian_matrix[i][j] = 0;*/
            opp_d_matrix[i][j] = 0;
            if (i == j)
                union_matrix[i][j] = 1; 
            else
                union_matrix[i][j] = 0;
        }
    }


    /* calculate D^(-0.5) */
    for (i = 0; i < n; i++)
        opp_d_matrix[i][i] = 1/(pow((d_matrix[i][i]),0.5));

    /*calculate D^(-1/2) * W * D^(-1/2):*/

    /* D^(-1/2) * W */
    matrixMultipication(n,vec_length,opp_d_matrix,w_matrix,result_matrix);
    /*(D^(-1/2) * W) * D ^ (-1/2) */
    matrixMultipication(n,vec_length,result_matrix,opp_d_matrix,temp_matrix);

    /*calculate final L_norm */
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            laplacian_matrix[i][j] = union_matrix[i][j] - temp_matrix[i][j];

    for(i = 0; i < n ; i++){
        free(union_matrix[i]);
        free(opp_d_matrix[i]);
        free(result_matrix[i]);
        free(temp_matrix[i]);
    }
    free(union_matrix);
    free(opp_d_matrix);
    free(result_matrix);
    free(temp_matrix);

    return laplacian_matrix;
}

int eigengapHeuristic(struct eigens* eigensArray, int n){
    double* deltas;
    int i, argmax_i;
    /*case we need to use 1.3 The Eigengap Heuristic because the input k = 0 */
    deltas = (double*)calloc(n-1, sizeof(double));
    /*Calculating the deltas */
    for (i = 0; i < n - 1; i++)
        deltas[i] = fabs(eigensArray[i].value - eigensArray[i+1].value);
    
    /*k = argmax_i(delta_i), i = 1,...,n/2 */
    argmax_i = deltas[0];
    for (i = 0; i < (int)(n/2); i++){
        if (deltas[i] > argmax_i)
            argmax_i = deltas[i];
    }
    
    free(deltas);

    return (argmax_i+1);
}

void jacobi(double** vectors_matrix, int n){
    int i;
    struct eigens* eigens_arr;

    eigens_arr = jacobiCalc(vectors_matrix, n, 1, 1);
    for(i=0; i<n; i++)
        free(eigens_arr[i].vector);
    free(eigens_arr);
}

/* print = 1 means do print, print = 0 means do not */
struct eigens* jacobiCalc(double** A_matrix, int n, int print, int sort){
    double off;
    int i, lar_i, lar_j, rotations_number;
    double** P_matrix;
    double** Atag_matrix;
    double** V_matrix;
    double** temp_matrix;
    double** V_trans_matrix;
    struct eigens* eigensArray;
    int* lar_arr;

    rotations_number = 0;
    temp_matrix = allocateMem(n,n);
    V_matrix = createMatrixI(n);
    eigensArray = (eigens*) calloc(n, sizeof(struct eigens));
    if (eigensArray == NULL)
        errorOccured();
    
    while(rotations_number < 100){
    /* find latgest item and it's i,j */
    lar_arr = RetrieveLargestIndexes(A_matrix, n);
    lar_i = lar_arr[0];
    lar_j = lar_arr[1];

    /* construct P */
    P_matrix = createMatrixP(A_matrix, n, lar_i, lar_j);  
    
    /* calculate A' (according to step 6)*/
    Atag_matrix = createMatrixAtag(A_matrix, n, lar_i, lar_j);
    
    /*result saved in temp matrix */
    matrixMultipication(n, n, V_matrix, P_matrix, temp_matrix);

    /* V =  temp */
    copyMatrices(V_matrix, temp_matrix, n);

    /* calc off(A) - off(A') */
    off = offFunc(A_matrix, Atag_matrix, n); 
    if (off <= EPSILON)
        break;
    
    /* A = A' */
    copyMatrices(A_matrix, Atag_matrix, n); 

    rotations_number++;
    } 

    V_trans_matrix = matrixTranspose(V_matrix, n, n);

    for(i = 0; i < n; i++){
        eigensArray[i].index = i;
        eigensArray[i].value = A_matrix[i][i];
        eigensArray[i].vector = copyToEigenVectors(V_trans_matrix[i], n);
    }

    if (print == 1)
        printJacobi(eigensArray, n);

    for(i = 0; i < n ; i++){
        free(P_matrix[i]);
        free(Atag_matrix[i]);
        free(temp_matrix[i]);
        free(V_matrix[i]);
        free(V_trans_matrix[i]);
    }
    free(V_trans_matrix);
    free(P_matrix);
    free(temp_matrix);
    free(Atag_matrix);
    free(lar_arr);
    if (sort == 1)
        qsort(eigensArray, n, sizeof(struct eigens), comparator);
    return eigensArray;
}

void copyMatrices(double** matrix_A, double** matrix_B, int n){
    int i,j;
    for (i=0; i<n; i++){
        for (j=0; j<n; j++){
            matrix_A[i][j] = matrix_B[i][j];
        }
    }
}

double** createMatrixI(int n){
    double** I_matrix;
    int i;

    I_matrix = allocateMem(n,n);
    for (i=0; i<n; i++)
        I_matrix[i][i]=1;

    return I_matrix;
}

/* given a matrix this function return the transpose matrix */
double** matrixTranspose(double** mat, int n, int vec_length){
    int i, j;
    double** trans_matrix;

    trans_matrix = allocateMem(n, vec_length);
    if (trans_matrix == NULL)
        errorOccured();

    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            trans_matrix[i][j] = mat[j][i];
        }
    }
    return trans_matrix;
}

/* comparator func in order to sort a list of eigen structs */
int comparator(const void* first, const void* second)
{
    struct eigens* e1;
    struct eigens* e2; 

    e1 =  (struct eigens*) first;
    e2 =  (struct eigens*) second;

    if (e1->value > e2->value)
        return -1;
    else if (e1->value < e2->value) 
        return 1;
    else if (e1->index < e2->index) 
        return -1;
    return 0;
}

/* for spk in cmodule */
double** createMatrixT(struct eigens* eigensArray, int n, int k){
    int i;
    int j;
    double sum;
    double** t_matrix_eigen;
    double** matrix_U;

    t_matrix_eigen = allocateMem(n, k);
    if (t_matrix_eigen == NULL)
        errorOccured();
    matrix_U = allocateMem(n, k);
    if (matrix_U == NULL)
        errorOccured();

    for(i = 0; i < k; i++){
        for(j = 0; j < n; j++){
             matrix_U[j][i] = eigensArray[i].vector[j];
        }
    }

    for(i = 0; i < n; i++){
        sum = 0;
        for(j = 0; j < k; j++){
            sum += matrix_U[i][j] * matrix_U[i][j];
        }
        sum = sqrt(sum);

        for(j = 0; j < k; j++){
            if(sum != 0){
                t_matrix_eigen[i][j] = ((matrix_U[i][j]) / (sum));
            }
            else{
                t_matrix_eigen[i][j] = 0.0;
            }
        }
    }

    for(i=0; i<n; i++){
        free(matrix_U[i]);
    }
    free(matrix_U);

    return t_matrix_eigen;
}

/* make a vector out of a row in a matrix */
double* copyToEigenVectors(double* vec_matrix, int n){
    int i;
    double* vector;

    vector = (double*) calloc(n, sizeof(double));
    if (vector == NULL)
        errorOccured();
    
    for(i = 0; i < n; i++){
        vector[i] = vec_matrix[i];
    }
    return vector;
}

/* calculate the multiplication of two matrices */
void matrixMultipication(int rows_num, int columns_num, double** mat1, double** mat2, double** result)
{
    int i, j, k;

    for (i = 0; i < rows_num; i++){
        for (j = 0; j < columns_num; j++){
            result[i][j] = 0;
            for (k = 0; k < rows_num; k++){
                result[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
}

/* the function calculate off(A)^2 and off(A')^2*/
double offFunc(double** A_matrix, double** Atag_matrix, int n){
    double off_a, off_a_tag;
    int i, j;
    
    off_a = 0;
    off_a_tag = 0;
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            if (i==j)
                continue;
            off_a += pow(A_matrix[i][j], 2); 
            off_a_tag += pow(Atag_matrix[i][j], 2);
        }
    }
    return (off_a - off_a_tag);
}

double** jacobiMatForPrint(struct eigens* eigensArray, int n){
    double** jacobi_matrix;
    int i, j;

    jacobi_matrix = allocateMem(n, n);
    if (jacobi_matrix == NULL)
        errorOccured();

    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++)
            jacobi_matrix[j][i] = eigensArray[i].vector[j];
    }
    return jacobi_matrix;
}

void printJacobi(struct eigens* eigensArray, int n){
    int i;
    double** jacobi_matrix;

    for (i = 0; i < n; i++){
         if(eigensArray[i].value < 0 && eigensArray[i].value > -0.00005)
            printf("0.0000");
        else
            printf("%.4f", eigensArray[i].value);
        if(i != n - 1)
            printf(",");
    }
    printf("\n");
    jacobi_matrix = jacobiMatForPrint(eigensArray, n);
    printMatrix(jacobi_matrix, n, n);
    for(i=0; i < n; i++)
        free(jacobi_matrix[i]);
    free(jacobi_matrix);
}

/* function that print the given matrix. */
void printMatrix(double** matrix, int n, int vector_length){
    int i, j;
    
    for(i = 0; i < n; i++){
        for(j = 0; j < vector_length; j++) {
            if (matrix[i][j] < 0 && matrix[i][j] > -0.00005) {
                printf("0.0000");
            } else {
                printf("%.4f", matrix[i][j]);
            }
            if(j != vector_length - 1){
                printf(",");
            }
        }
        printf("\n");
    }
}

/* this function allocates memory if succeeded, and if not returns */
double** allocateMem(int n, int vector_length)
{ 
    int i;
    int j;
    double** matrix;

    matrix = (double **)calloc(n, sizeof(double *));
    if (matrix == NULL)
    {
        return NULL;
    }

    for (i = 0; i < n; i++)
    {
        matrix[i] = (double *)calloc(vector_length, sizeof(double));
        if (matrix[i] == NULL)
        {
            for (j = 0; j < vector_length; j++)
            {
            free(matrix[j]);
            }
            free(matrix); 
            return NULL;
        }
    }
    return matrix;
}

void errorOccured(){
    printf("An Error Has Occurred");
    exit(1);
}

int* RetrieveLargestIndexes(double** A_matrix, int n){
    double largest;
    int i, j, lar_i, lar_j;
    int* indexes;

    lar_i = 0;
    lar_j= 1;
    indexes = (int*)calloc(2, sizeof(int));
    if (indexes == NULL)
        return NULL;   
    largest = fabs(A_matrix[0][1]);

    for (i = 0; i < n; i++){
        for (j = i + 1; j < n; j++){
            if (fabs(A_matrix[i][j]) > fabs(largest)){
                /* 3. find the Pivot A_ij */
                largest = fabs(A_matrix[i][j]);
                lar_i = i;
                lar_j = j;
            }
        }
    }
    indexes[0] = lar_i;
    indexes[1] = lar_j;
    return indexes;
}

double retrieveTheta(double** A_matrix, int lar_i, int lar_j){
    return (A_matrix[lar_j][lar_j] - A_matrix[lar_i][lar_i])/(2*A_matrix[lar_i][lar_j]);
}

double retrieveT(double theta){
    return sign(theta)/(fabs(theta) + sqrt((pow(theta, 2)) + 1));
}

double retrieveC(double t){
    return 1/sqrt((pow(t, 2)) + 1);
}

double retrieveS(double t, double c){
    return t*c;
}

int sign(double theta){
    if(theta >= 0) /*CHECK if ?: possible*/
        return 1;
    else 
        return -1;
}

/* this function construct P and initialize */
double** createMatrixP(double** A_matrix,int n, int lar_i, int lar_j){
    int i, j;
    double theta, t, c, s;
    double** P_matrix;
    theta = retrieveTheta(A_matrix, lar_i, lar_j);
    t = retrieveT(theta);
    c = retrieveC(t);
    s = retrieveS(t,c);
    P_matrix = allocateMem(n, n);
    if (P_matrix == NULL)
        errorOccured();
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            if (i == j){
                P_matrix[i][j] = 1;
            }
        }
    }
    P_matrix[lar_i][lar_i] = c;
    P_matrix[lar_j][lar_j] = c;
    P_matrix[lar_i][lar_j] = s;
    P_matrix[lar_j][lar_i] = -s;

    return P_matrix;
}

/* this function calculate A' (according to step 6) */
double** createMatrixAtag(double** A_matrix,int n,int lar_i,int lar_j){
    int i, j;
    double theta, t, c, s;
    double** Atag_matrix;
    theta = retrieveTheta(A_matrix,lar_i,lar_j);
    t = retrieveT(theta);
    c = retrieveC(t);
    s = retrieveS (t,c);    
    Atag_matrix = allocateMem(n,n);
    for (i=0; i<n; i++){
        for (j=0; j<n; j++){
        Atag_matrix[i][j] = A_matrix[i][j];
        }
    }

    for(i=0; i<n; i++){
        if((i != lar_i) || (i != lar_j)){/*CHECK if needed*/
                Atag_matrix[i][lar_i] = c*A_matrix[i][lar_i] - s*A_matrix[i][lar_j];
                Atag_matrix[lar_i][i] = c*A_matrix[i][lar_i] - s*A_matrix[i][lar_j];
                Atag_matrix[i][lar_j] = c*A_matrix[i][lar_j] + s*A_matrix[i][lar_i];
                Atag_matrix[lar_j][i] = c*A_matrix[i][lar_j] + s*A_matrix[i][lar_i];
        }
    }
    Atag_matrix[lar_i][lar_i] = c*c*A_matrix[lar_i][lar_i] + s*s*A_matrix[lar_j][lar_j] - 2*s*c*A_matrix[lar_i][lar_j];
    Atag_matrix[lar_j][lar_j] = s*s*A_matrix[lar_i][lar_i] + c*c*A_matrix[lar_j][lar_j] + 2*s*c*A_matrix[lar_i][lar_j];
    Atag_matrix[lar_i][lar_j] = 0;
    Atag_matrix[lar_j][lar_i] = 0;

    return Atag_matrix;
}

/* K-means Functions */
void getFinalCentroids(double **centroids, double **elements, int k, int d, int n, int max_iter, double epsilone)
{
    int bit;
    int i;
    int iteration_number;
    int* elements_location;
    int* items_number_clusters;
    double** old_centroids;

    bit = 1;
    iteration_number = 0;

    elements_location = (int*)calloc(n, sizeof(int));
    if (!elements_location)
        errorOccured();

    items_number_clusters = (int*)calloc(k, sizeof(int));
    if (!items_number_clusters)
        errorOccured();

    old_centroids = (double**)calloc(k, sizeof(double*));
    if (!old_centroids)
        errorOccured();

    for (i=0;i<k;i++){
        old_centroids[i]=(double*) calloc(d,sizeof(double));
            if (!old_centroids[i])
                errorOccured();
    }

    while (bit==1 && max_iter>iteration_number){
        iteration_number++;
        initClusters(elements_location, items_number_clusters, n, k);
        assignCentroids(elements, centroids, items_number_clusters,elements_location,k,d,n);
        saveCentroids(old_centroids, centroids, k, d);
        resetCentroids(centroids, k, d);
        updateCentroids(centroids,elements,items_number_clusters,elements_location,d,n,k);
        bit = convergence(old_centroids, centroids, k, d, epsilone);
    }

    free(elements_location);
    free(items_number_clusters);
    for(i=0; i<k; i++)
        free(old_centroids[i]);
    free(old_centroids);
}

void resetCentroids(double** centroids,int k, int d){
    int i;
    int j;
    for (i=0;i<k;i++)
        for(j=0;j<d;j++)
            centroids[i][j]=0.0;                
}

void initClusters(int* elements_loc, int* items_number_clusters, int n, int k){
    int i;
    for (i=0;i<n;i++)
        elements_loc[i]=0;
    for (i=0;i<k;i++)
        items_number_clusters[i]=0;
}

void saveCentroids(double** old_centroids, double** centroids, int k, int d){
    int i;
    int j;
    for (i=0;i<k;i++){
            for (j=0;j<d;j++)
                old_centroids[i][j] = centroids[i][j];
        }
}

void assignCentroids(double **ele,double **cntrds,int* in_clstrs,int* ele_loc,int k,int d,int n)
{
    int i, j, l;
    double sum=0.0;
    double min;
    int flag;
    int min_index;

    for (i=0; i < n; i++){
        min_index = 0;
        flag = 0;
        for (j=0; j < k; j++){
            sum = 0.0;
                for (l = 0; l < d; l++)
                    sum += pow((ele[i][l] - cntrds[j][l]),2);
            sum = pow(sum,0.5);
            if (flag == 0){
                min = sum;
                flag = 1;
            }
            else if (sum < min){
                min = sum;
                min_index = j;
            }
        }
        in_clstrs[min_index]+=1;
        ele_loc[i]=min_index;
    }
}

void updateCentroids(double** cntrds,double** ele,int* in_clstrs,int *ele_loc,int d,int n,int k)
{
    int m, i, q;

    for (i=0;i<k;i++){
            for (m=0;m<n;m++){
                if (ele_loc[m]==i){
                    for (q=0;q<d;q++)
                        cntrds[i][q]+=ele[m][q];
                }
            }
            for (q=0;q<d;q++){
                if (cntrds[i][q]==0){
                    continue;
                }
                cntrds[i][q]=cntrds[i][q]/(double)in_clstrs[i];
            }
        }
}

int convergence(double** old_centroids,double** centroids,int k,int d,int eps){
    int bit;
    int i, j;
    double sum;

    bit = 0;
    for(i=0; i<k; i++){
        sum = 0.0;
        for (j = 0; j<d; j++)
            sum += pow((old_centroids[i][j] + centroids[i][j]), 2);
        sum = pow(sum, 0.5);
        if (sum >= eps)
            bit = 1;
    }
    return bit;
}
