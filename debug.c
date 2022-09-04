#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define EPSILON 0.00001

void wam(double** vectors_matrix, int n, int vec_length);
double** wam_calc(double** vectors_matrix, int n, int vec_length);
void print_matrix(double** matrix, int n, int vector_length);
double** allocateMem(int n, int vector_length);
void error_occurred();

int main(int argc, char* argv[]){
    int n;
    int i;
    int j;
    int r;
    int vec_length;
    char* goal;
    char* input_filename;
    double** vectors_matrix;
    FILE* temp_file;
    double vector_item;
    char comma;
    if (argc != 3){
        printf ("Invalid Input!");
        exit(1);
    }

    goal = argv[1]; /* the enum type */
    input_filename = argv[2]; /* the input */

    temp_file = fopen(input_filename, "r");
    if (temp_file == NULL)
        error_occurred();

    /* counts the num of items in vectors, vec_length and number of vector */
    vec_length=0;
    n=0;
    while(fscanf(temp_file, "%lf%c",&vector_item,&comma) == 2){
        if (n == 0){
            vec_length++;}
        if (comma == '\n'){
            n++;}
    }
    fseek(temp_file,0,0);

    /* initiallization of vectors matrix */
    vectors_matrix = allocateMem(n, vec_length);
    if (vectors_matrix == NULL)
        error_occurred();

    /* insert vector's items to vectors_matrix */
    i=0;
    j=0;
    while(fscanf(temp_file, "%lf%c",&vector_item,&comma) == 2){
        vectors_matrix[i][j] = vector_item;
        j++;
        if (comma == '\n')
        {
            j=0;
            i++;
        }
    }
    r = fclose(temp_file);
    if (r!=0){
        error_occurred();
    }
    printf("file was read\n");
    printf("n is: %d\n",n);
    printf("vec_length is: %d\n",vec_length);
    /* check which goal to choose */
        if(strcmp(goal, "wam") == 0){
             printf("identified wam\n");
            wam(vectors_matrix, n, vec_length);
        }
    for(i=0; i<n; i++)
        free(vectors_matrix[i]);
    free(vectors_matrix);

    return 0;
}


void wam(double** vectors_matrix, int n, int vec_length){
    int i;
    double** w_matrix;
    w_matrix = wam_calc(vectors_matrix, n, vec_length);
     printf("wam_calc finished\n");
    print_matrix(w_matrix, n, n);
    printf("matrix print finished\n");
    for(i = 0; i < n; i++)
        free(w_matrix[i]);
    free(w_matrix);
    printf("memory freed\n");
}

double** wam_calc(double** vectors_matrix, int n, int vec_length){
    int i;
    int j;
    int s;
    double sum;
    double** w_matrix;

    w_matrix = allocateMem(n, n);
    if (w_matrix == NULL)
        error_occurred();  

    /* calculating the values in each cell */
    for (i = 0; i < n; i++)
        w_matrix[i][i] = 0; /* the diagonal line in matrix's values are 0 */
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

/* function that print the given matrix. */
void print_matrix(double** matrix, int n, int vector_length){
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

void error_occurred(){
    printf("An Error Has Occurred");
    exit(1);
}