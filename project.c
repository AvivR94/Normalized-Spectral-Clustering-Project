#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define EPSILON = 1.0*pow(10,-5)

typedef struct eigens{
    int index;
    double value;
    double* vector;
}eigens;

int main(int argc, char* argv[]){
if (agrc != 4){
    print ("Invalid Input");
    return 0;
}
k = atoi(argv[1]);
char* goal = sys.argv[2] //the enum type
input_filename = sys.argv[3]; //the input

temp_file = fopen(input_filename, "r");
if (temp_file == NULL){
    printf("An Error Has Occurred");
    return 0;
}

//count what is the number of vectors n and their length vec_length
first = 0;
while(fscanf(temp_file, "%lf%c",&vector_item,&comma) == 2)
{
    if (comma == '\n'){
        n++;
        if (first = 0){
                vec_length ++;
                first = 1
            }
    }
    else if (first == 0){
        vec_length ++;
    }
}

//initiallization of vectors matrix
double** vectors_matrix = (double**)calloc(n, sizeof(double*));

//insert vector's items to vectors_matrix
i = 0;
j = 0;
while(fscanf(temp_file, "%lf%c",&vector_item,&comma) == 2){
    vectors_matrix[i][j] = vector_item;
    j++;
    if (comma == '\n')
    {
        j=0;
        i++;
    }
}

//check which goal to choose
    if (strcmp(goal, "spk") == 0){
        spk(vectors_matrix, n, vec_length, sys.argv[3]);    
    else if(strcmp(goal, "wam") == 0){
        wam(vectors_matrix, n, vec_length);
    }
    else if (strcmp(goal, "ddg") == 0){
        ddg(vectors_matrix, n, vec_length);
    }
    else if (strcmp(goal, "lnorm") == 0){
        lnorm(vectors_matrix, n, vec_length);
    }
    else if (strcmp(goal, "jacobi") == 0){
        funcJacobi(vectors_matrix, n, vec_length);
    }


for(i=0; i<n; i++){
       free(vectors_matrix[i]);
   }
free(vectors_matrix);

fclose(temp_file);
return 0;
}

//###########################################################################
//calculate the multiplication of two matrices
def matrix_multiplication(int rows_num, int columns_num, double** mat1, double** mat2, double** result):
    for (i = 0; i < rows_num; i++)
        for (j = 0; j < columns_num; j++)
            for (k = 0; k < rows_num; k++)
                result[i][j] += mat1[i][k] * mat2[k][j];

//###########################################################################
//free the allocations
void double_free(double** p, int size_of_block){
    int i;
    for (i = 0; i < size_of_block; i++){
       free(p[i]);
       }
    free(p);
}


//##################################################################
//case wam: 1.1.1 The Weighted Adjacency Matrix
void wam(double** vectors_matrix, int n, int vec_length){

    //initiallization of matrix W
    double** w_matrix = (double**)calloc(n, sizeof(double*));

    //calculating the values in each cell
    int sum = 0;
    for (i = 0; i < n; i++)
        w_matrix[i][i] = 0 //the diagonal line in matrix's values are 0
    for (i = 0; i < n; i++){
        for (j= i + 1; j < n; j++){
            sum = 0;
            for (s = 0; s < d; s++)
                sum += pow((vectors_matrix[i][s] - vectors_matrix[j][s]),2)
            w_matrix[i][j] = math.exp(pow(sum,(0.5))/-2);
            w_matrix[j][i] = math.exp(pow(sum,(0.5))/-2);
        }
    }

    for(i = 0; i < n; i++){
        free(w_matrix[i]);
    }

    free(w_matrix);

}
//###########################################################################
//case ddg: 1.1.2 The Diagonal Degree Matrix
void ddg(double** vectors_matrix, int n, int vec_length){
    //initiallization
    double** d_matrix = (double**)calloc(n, sizeof(double*));
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            d_matrix[i][j] = 0;

    //calculate the sum of each row and place in the diagonal line of d_matrix
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++)
            d_matrix[i][i] += w_matrix[i][j];
    }

    for(i = 0; i < n ; i++){
        free(d_matrix[i]);
    }
    free(d_matrix);
}
        

//###########################################################################
//case lnorm: 1.1.3  The Normalized Graph Laplacian

void lnorm(double** vectors_matrix, int n, int vec_length){

    //Laplacian_matrix initiallization
    double** laplacian_matrix = (double**)calloc(n, sizeof(double*));
    double** union_matrix = (double**)calloc(n, sizeof(double*));
    double** opp_d_matrix = (double**)calloc(n, sizeof(double*));
    double** result_matrix = (double**)calloc(n, sizeof(double*));
    double** temp_matrix = (double**)calloc(n, sizeof(double*));
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            laplacian_matrix[i][j] = 0;
            opp_d_matrix[i][j] = 0;
            if (i == j)
                union_matrix[i][j] = 1; 
            else
                union_matrix[i][j] = 0;
        }
    }


    //calculate D^(-0.5)
    for (i = 0; i < n; i++)
        opp_d_matrix[i][i] = 1/(pow((d_matrix[i][i]),0.5))

    //calculate D^(-1/2) * W * D^(-1/2)
    matrix_multiplication(n, vec_length, opp_d_matrix, w_matrix, result_matrix) ; //D^(-1/2) * W
    matrix_multiplication(n, vec_length, result_matrix, opp_d_matrix, temp_matrix) ;//(D^(-1/2) * W) * D ^ (-1/2)

    //calculate final L_norm
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            laplacian_matrix[i][j] = union_matrix[i][j] - temp_matrix[i][j];

    for(i = 0; i < n ; i++){
        free(laplacian_matrix[i]);
        free(union_matrix[i]);
        free(opp_d_matrix[i]);
        free(result_matrix[i]);
        free(temp_matrix[i]);
    }
    free(laplacian_matrix);
    free(union_matrix);
    free(opp_d_matrix);
    free(result_matrix);
    free(temp_matrix);
}
###########################################################################
void eigengap_heuristic(double** vectors_matrix, int n, int vec_length){
    //case we need to use 1.3 The Eigengap Heuristic because the input k = 0
    double* eigenvalues = (int*)calloc(n, sizeof(double));
    double* deltas = (int*)calloc(n, sizeof(double));
    int k=0;
    int cnt = 0
    //insert all eigenvalues to eigenvalues array
    double max = -1;
    int argmax_i=0;
    for (i = 0; i < n; i++){
        eigenvalues[i] = laplacian_matrix[i][i]
    }

    qsort(eigensArray, n, sizeof(struct eigens), comparator);

  
    //Calculating the deltas
    for (int i = 0; i < n - 1; i++)
        deltas[i] = abs(eigenvalues[i] - eigenvalues[i+1]);
    
    //k = argmax_i(delta_i), i = 1,...,n/2
    argmax_i = deltas[0];
    for (int i = 0; i < int(n/2); i++)
        if (deltas[i] > argmax_i)
            argmax_i = deltas[i];

    k = argmax_i;

    free(eigenvalues);
    free(deltas);

    return k;
}
//###########################################################################
//jacobi algorithm
void jacobi(double** vectors_matrix, int n, int vec_length){

    double off_a = 0;
    double off_a_tag= 0;
    double largest = 0;
    int largest_i = -1;
    int largest_j = -1;
    double theta = 0;
    double t = 0;
    double c = 0;
    double s = 0;
    int first = 0;
    int rotations_number = 0;
    double** p_matrix = (double**)calloc(n, sizeof(double*));
    double** to_copy = (double**)calloc(n, sizeof(double*));
    double** temp = (double**)calloc(n, sizeof(double*));
    double** v_matrix = (double**)calloc(n, sizeof(double*));

    do
    {
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
            {
                if (abs(vectors_matrix[i][j]) > abs(largest))
                {
                    //3. find the Pivot A_ij
                    largest = vectors_matrix[i][j];
                    largest_i = i;
                    largest_j = j;
                }
            }

    //4. Obtain c,t
    theta = (vectors_matrix[j][j] - vectors_matrix[i][i])/(2*vectors_matrix[i][j]);
    t = 1/(abs(theta) + pow((pow(theta, 2) + 1) , 0.5));
    if (theta < 0)
        t = t * -1;
    c = 1/pow((pow(t, 2) + 1) , 0.5);
    s = t * c;

    //P construction and initiallization of to_copy array which will become A'

    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            to_copy[i][j] = vectors_matrix[i][j];
            if (i == j){
                p_matrix[i][j] = 1;
            }
        }
    }
    p_matrix[largest_i][largest_i] = c;
    p_matrix[largest_j][largest_j] = c;
    if (largest_i<largest_j){
        p_matrix[largest_i][largest_j] = s;
        p_matrix[largest_j][largest_i] = -s;
    }
    else{
        p_matrix[largest_i][largest_j] = -s;
        p_matrix[largest_j][largest_i] = s;
    }

    //save P matrix
    if (first == 0)
    {
        first = 1;
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                v_matrix[i][j] = p_matrix[i][j];
            }
        }
    }

    else //multiple multipication of privious P's with current * p_matrix 
    {
        matrix_multiplication(n, vec_length, v_matrix, p_matrix, temp)
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                v_matrix[i][j] = temp_matrix[i][j];
            }
        }
    }

    //calculate A' (according to step 6)

    to_copy[largest_i][largest_i] = pow(c, 2) * vectors_matrix[largest_i][largest_i] + pow(s, 2) * vectors_matrix[largest_j][largest_j] - 2 * s * c * vectors_matrix[largest_i][largest_j] + pow(c, 2);
    to_copy[largest_j][largest_j] = pow(s, 2) * vectors_matrix[largest_i][largest_i] + pow(c, 2) * vectors_matrix[largest_j][largest_j] + 2 * s * c * vectors_matrix[largest_i][largest_j] + pow(c, 2);
    to_copy[largest_i][largest_j] = 0;
    to_copy[largest_j][largest_i] = 0;
    for (int r = 0; r < n; r++)
        for (int m = r + 1; m < n; m++)
        {
            if ( r != largest_i && r != largest_j)
            {
                if (m == largest_i)
                {
                    to_copy[r][m] = c * vectors_matrix[r][largest_i] - s * vectors_matrix[r][largest_j];
                    to_copy[m][r] = c * vectors_matrix[r][largest_i] - s * vectors_matrix[r][largest_j];
                }
                if (m == largest_j)
                {
                    to_copy[r][m] = c * vectors_matrix[r][largest_j] - s * vectors_matrix[r][largest_i];
                    to_copy[m][r] = c * vectors_matrix[r][largest_j] - s * vectors_matrix[r][largest_i];
                }
            }
            to_copy[i][j] = vectors_matrix[i][j]
        }


    // calculate off(A)^2 and off(A')^2

    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            off_a += pow(vectors_matrix[i][j], 2); 
            off_a_tag += pow(to_copy[i][j], 2);
        }
    }

    //copy the values of to_copy to original array( A = A')
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            vectors_matrix[i][j] = to_copy[i][j];
        }
    }

        rotations_number++;
    
    } while((off_a - off_a_tag) > EPSILON && rotations_number <= 100);


    v_trans_matrix = matrix_Transpose(v_matrix, n);
    
    eigensArray = (eigens*) calloc(n, sizeof(struct eigens));
    assert(eigensArray != NULL);

    for(i = 0; i < n; i++){
        eigensArray[i].index = i;
        eigensArray[i].value = vectors_matrix[i][i];
        eigensArray[i].vector = copy_to_eigen_Values(v_trans_matrix[i], n);
    }


    for(i = 0; i < n ; i++){
        free(p_matrix[i]);
        free(to_copy[i]);
        free(temp[i]);
        free(v_matrix[i]);
        free(v_trans_matrix[i]);
    }
    free(v_trans_matrix)
    free(p_matrix);
    free(union_matrix);
    free(temp);
    free(v_matrix);

    //to add prints

    qsort(eigensArray, n, sizeof(struct eigens), comparator);
    return eigensArray;
}

//###############################################################
//given a matrix this function return the transpose matrix
double** matrix_Transpose(double** mat, int n){
    trans_matrix = (double**) calloc(n, sizeof(double*));
    assert(trans_matrix != NULL);

    for(i = 0; i < n; i++){
        trans_matrix[i] = (double*) calloc(n, sizeof(double));
        assert(trans_matrix[i] != NULL);
    }

    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            trans_matrix[i][j] = mat[j][i];
        }
    }

    return trans_matrix;
}

//###############################################################

void spk(double** vectors_matrix, int n, int vec_length, int k){
    //to build cover functions to each goal
    w_matrix = wam(vectors_matrix, n, vec_length);
    d_matrix = ddg(w_matrix, n, vec_length);
    l_matrix = lnorm(d_matrix, n, vec_length);
    eigens_arr = jacobi(l_matrix, n, vec_length);

    if(k == 0){
        k = eigengap_heuristic(l_matrix, n, vec_length);
    } 

    //k means issues to be completed
    funcJacobi(vectors_matrix, n, vec_length);
        for(i = 0; i < n ; i++){
        free(w_matrix[i]);
        free(d_matrix[i]);
        free(l_matrix[i]);
    }
    free(w_matrix);
    free(d_matrix);
    free(l_matrix);
}


 //comparator func in order to sort a list of eigen structs

int comparator (const void* first, const void* second)
{
    struct eigens* e1 =  (struct eigens*) first;
    struct eigens* e2 =  (struct eigens*) second;

    if (e1->value > e2->value)
        return 1;
    else if (e1->value < e2->value) 
        return -1;
    else if (e1->index < e2->index) 
        return -1;
    else return 0;
}


//
double** build_matrix_t_eigen(struct eigens* eigensArray, int n, int k){
    t_matrix_eigen = (double**) calloc(n, sizeof(double*));
    assert(matTeigen != NULL);
    matU = (double**) calloc(n, sizeof(double*));
    assert(matU != NULL);

    for(i = 0; i < n; i++){
        t_matrix_eigen[i] = calloc(k, sizeof(double));
        assert(t_matrix_eigen[i] != NULL);
        matU[i] = calloc(k, sizeof(double));
        assert(matU[i] != NULL);

    }

    for(i=0; i<k; i++){
        for(j=0; j<n; j++){
             matU[j][i] = eigensArray[i].vector[j];
        }
    }

    for(i=0; i<n; i++){
        dist = 0;
        for(j=0; j<k; j++){
            dist += matU[i][j] * matU[i][j];
        }
        dist = sqrt(dist);
        for(j=0; j<k; j++){
            if(dist != 0){
                t_matrix_eigen[i][j] = ((matU[i][j]) / (dist));
            } else {
                t_matrix_eigen[i][j] = 0.0;
            }
        }
    }

    for(i=0; i<n; i++){
        free(matU[i]);
    }
    free(matU);

    return t_matrix_eigen;
}

 //copy vector into eigen struct.
 
double* copy_to_eigen_Values(double* vec_matrix, int n){
    vector = (double*) calloc(n, sizeof(double));
    assert(vector != NULL);
    for(j = 0; j < n; j++){
        vector[j] = vec_matrix[j];
    }
    return vector;
}

}