#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define EPSILON = 0.00001

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

//this function counts the num of items in vectors, vec_length
n = get_vec_num(temp_file);
vec_length = get_vec_length(temp_file);

//initiallization of vectors matrix
double** vectors_matrix = allocateMem(n, vec_length);
if (vectors_matrix == NULL)
    error_occurred();


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
    }  
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
        jacobi(vectors_matrix, n, vec_length);
    }


for(i=0; i<n; i++){
       free(vectors_matrix[i]);
   }
free(vectors_matrix);

fclose(temp_file);
return 0;
}


void wam(double** vectors_matrix, int n, int vec_length){
{
    w_matrix = wam_calc(vectors_matrix, n, vec_length);

    Print_matrix(w_matrix, n, n);

    for(i = 0; i < n; i++){
        free(w_matrix[i]);
    }
    free(w_matrix);
}


double** wam_calc(double** vectors_matrix, int n, int vec_length){

    //initiallization of matrix W
    double** w_matrix = allocateMem(n, vec_length);
    if (w_matrix == NULL)
        error_occurred();  

    //calculating the values in each cell
    double sum = 0;
    for (int i = 0; i < n; i++)
        w_matrix[i][i] = 0 //the diagonal line in matrix's values are 0
    for (int i = 0; i < n; i++){
        for (int j= i + 1; j < n; j++){
            sum = 0;
            for (int s = 0; s < vec_length; s++)
                sum += pow((vectors_matrix[i][s] - vectors_matrix[j][s]),2);
            w_matrix[i][j] = math.exp(pow(sum,(0.5))/-2);
            w_matrix[j][i] = math.exp(pow(sum,(0.5))/-2);
        }
    }
    return w_matrix;
}

void ddg(double** vectors_matrix, int n, int vec_length){
    w_matrix = wam_calc(vectors_matrix, n, vec_length);
    d_matrix = ddg_calc(w_matrix, n);

    Print_matrix(d_matrix, n, n);

    for(i = 0; i < n ; i++){
        free(w_matrix[i]);
        free(d_matrix[i]);
    }
    free(w_matrix);
    free(d_matrix);
}

double** ddg_calc(double** vectors_matrix, int n, int vec_length){
    //initiallization
    double** d_matrix = allocateMem(n, vec_length);
    if (d_matrix == NULL)
        error_occurred();
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            d_matrix[i][j] = 0;

    //calculate the sum of each row and place in the diagonal line of d_matrix
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++)
            d_matrix[i][i] += w_matrix[i][j];
    }

    return d_matrix;
}
        

void lnorm(double** vectors_matrix, int n, int vec_length){
    w_matrix = wam_calc(vectors_matrix, n, vec_length);
    d_matrix = ddg_calc(w_matrix, n);
    lnorm_calc = lnorm_calc(vectors_matrix, n, vec_length);

    Print_matrix(lnorm_calc, n, n);

    for(i = 0; i < n; i++){
        free(w_matrix[i]);
        free(d_matrix[i]);
        free(lnorm_calc[i]);
    }
    free(w_matrix);
    free(d_matrix);
    free(lnorm_calc);}

double** lnorm_calc(double** vectors_matrix, int n, int vec_length){

    //Laplacian_matrix initiallization
    double** laplacian_matrix = allocateMem(n, vec_length);
    if (laplacian_matrix == NULL)
        error_occurred();
    double** union_matrix = allocateMem(n, vec_length);
    if (union_matrix == NULL)
        error_occurred();
    double** opp_d_matrix = allocateMem(n, vec_length);
    if (opp_d_matrix == NULL)
        error_occurred();
    double** result_matrix = allocateMem(n, vec_length);
    if (result_matrix == NULL)
        error_occurred();
    double** temp_matrix = allocateMem(n, vec_length);
    if (temp_matrix == NULL)
        error_occurred();

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
    
    //D^(-1/2) * W
    matrix_multiplication(n,vec_length,opp_d_matrix,w_matrix,result_matrix);
    //(D^(-1/2) * W) * D ^ (-1/2)
    matrix_multiplication(n,vec_length,result_matrix,opp_d_matrix,temp_matrix);

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

    return laplacian_matrix;
}


void eigengap_heuristic(double** vectors_matrix, int n, int vec_length){
    //case we need to use 1.3 The Eigengap Heuristic because the input k = 0
    double* eigenvalues = (double*)calloc(n, sizeof(double));
    double* deltas = (double*)calloc(n, sizeof(double));
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


void jacobi(double** vectors_matrix, int n, int vec_length){

    eigens_arr = jacobi_calc(vectors_matrix, n);

    print_jacobi(eigens_arr, n, vec_length);

    for(i = 0; i < n; i++){
        free(w_matrix[i]);
        free(d_matrix[i]);
        free(lnorm_calc[i]);
        free(eigens_arr[i]);

    }
    free(w_matrix);
    free(d_matrix);
    free(lnorm_calc);
    free(eigens_arr);}
}

struct eigens* jacobi_calc(double** vectors_matrix, int n){

    double theta = 0;
    double t = 0;
    double c = 0;
    double s = 0;
    int lar_arr[2] = {0,0};
    int cs_arr[2] = {0,0};
    int first = 0;
    int lar_i = -1;
    int lar_j = -1;
    double off;
    int rotations_number = 0;

    double** p_matrix = allocateMem(n, vec_length);
    if (p_matrix == NULL)
        error_occurred();
    double** to_copy = allocateMem(n, vec_length);
    if (to_copy == NULL)
        error_occurred();
    double** temp = allocateMem(n, vec_length);
    if (temp == NULL)
        error_occurred();
    double** v_matrix = allocateMem(n, vec_length);
    if (v_matrix == NULL)
        error_occurred();

    do
    {
    //find latgest item and it's i,j
    lar_arr = largest_indexes(vectors_matrix, n);
    lar_i = lar_arr[0];
    lar_j = lar_arr[1];

    //calculate c,s by theta and t
    cs_arr = get_c_and_s(vectors_matrix, lar_i, lar_j);
    c = cs_arr[0];
    s = cs_arr[1];

    p_mat_maker(vectors_matrix, p_matrix, lar_i, lar_j);  //construct P
    
    set_tocopy(vectors_matrix, to_copy);  //copy vectors_matrix to to_copy which will become to A'

    //this function calculate the current v
    first = v_calculation(vectors_matrix, p_matrix, v_matrix, temp, first);
    
    A_tag_calc(to_copy, vectors_matrix, lar_i, lar_j); //calculate A' (according to step 6)

    off = off_func(vectors_matrix, to_copy); //calc off(A) - off(A')

    A_to_A_tag(vectors_matrix, to_copy); //copy the values of to_copy to original array( A = A')

    rotations_number++;
    
    } while((off) > EPSILON && rotations_number <= 100);


    v_trans_matrix = matrix_Transpose(v_matrix, n);
    
    eigensArray = (eigens*) calloc(n, sizeof(struct eigens));
    if (eigensArray == NULL)
        error_occurred();  
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

    qsort(eigensArray, n, sizeof(struct eigens), comparator);
    return eigensArray;
}

//given a matrix this function return the transpose matrix
double** matrix_Transpose(double** mat, int n){
    trans_matrix = allocateMem(n, vec_length);
    if (trans_matrix == NULL)
        error_occurred();

    for(i = 0; i < n; i++){
        trans_matrix[i] = (double*) calloc(n, sizeof(double));
        if (trans_matrix[i] == NULL)
            error_occurred();
      }

    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            trans_matrix[i][j] = mat[j][i];
        }
    }
    return trans_matrix;
}

void spk(double** vectors_matrix, int n, int vec_length, int k){
    //to build cover functions to each goal
    w_matrix = wam_calc(vectors_matrix, n, vec_length);
    d_matrix = ddg_calc(w_matrix, n, vec_length);
    l_matrix = lnorm_calc(d_matrix, n, vec_length);
    eigens_arr = jacobi_calc(l_matrix, n, vec_length);

    if(k == 0){
        k = eigengap_heuristic(l_matrix, n, vec_length);
    } 

    build_matrix_t_eigen(eigensArray, n, k)
    kmeans //to be completed

    for(i = 0; i < n ; i++){
        free(w_matrix[i]);
        free(d_matrix[i]);
        free(l_matrix[i]);
        free(eigens_arr[i]);
    }
    free(w_matrix);
    free(d_matrix);
    free(l_matrix);
    free(eigens_arr);
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
    t_matrix_eigen = allocateMem(n, vec_length);
    if (t_matrix_eigen == NULL)
        error_occurred();
    u_matrix = allocateMem(n, vec_length);
    if (u_matrix == NULL)
        error_occurred();

    for(i = 0; i < n; i++){
        t_matrix_eigen[i] = calloc(k, sizeof(double));
        if (t_matrix_eigen[i] == NULL)
            error_occurred();  
        u_matrix[i] = calloc(k, sizeof(double));
        if (u_matrix[i] == NULL)
            error_occurred();    
        }

    for(i = 0; i < k; i++){
        for(j = 0; j < n; j++){
             u_matrix[j][i] = eigensArray[i].vector[j];
        }
    }

    int sum = 0;
    for(i = 0; i < n; i++){
        sum = 0;
        for(j = 0; j < k; j++)
            sum += u_matrix[i][j] * u_matrix[i][j];

        sum = sqrt(sum);

        for(j = 0; j < k; j++){
            if(dist != 0)
                t_matrix_eigen[i][j] = ((u_matrix[i][j]) / (sum));
            else
                t_matrix_eigen[i][j] = 0.0;
        }
    }

    for(i=0; i<n; i++){
        free(u_matrix[i]);
    }

    free(u_matrix);

    return t_matrix_eigen;
}

 //copy vector into eigen struct.
 
double* copy_to_eigen_Values(double* vec_matrix, int n){
    vector = (double*) calloc(n, sizeof(double));
    if (vector == NULL)
        error_occurred();
    for(j = 0; j < n; j++){
        vector[j] = vec_matrix[j];
    }
    return vector;
}

//this function counts what is the length of vectors
int get_vec_length(FILE *f){
    int vec_length = 0;
    int first = 0;
    while(fscanf(temp_file, "%lf%c",&vector_item,&comma) == 2)
    {
        if (comma == '\n'){
                    vec_length ++;
                    return vec_length;
        }
    }
}
//this function counts what is the number of vectors n
int get_vec_num(FILE *f){
    int n = 0;
    while(fscanf(temp_file, "%lf%c",&vector_item,&comma) == 2){
        if (comma == '\n'){
            n++;
        }
    return n
}

//calculate the multiplication of two matrices
void matrix_multiplication(int rows_num, int columns_num, double** mat1, double** mat2, double** result):
    for (i = 0; i < rows_num; i++)
        for (j = 0; j < columns_num; j++)
            for (k = 0; k < rows_num; k++)
                result[i][j] += mat1[i][k] * mat2[k][j];


//the function calculate off(A)^2 and off(A')^2 and return off_a - off_a_tag
double off_func(double** vectors_matrix[i][j], double** to_copy[i][j])
{
    double off_a = 0;
    double off_a_tag = 0;
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            off_a += pow(vectors_matrix[i][j], 2); 
            off_a_tag += pow(to_copy[i][j], 2);
        }
    }
    return off_a - off_a_tag;
}



double** jacobi_mat_for_print(struct eigens* eigensArray, int n, int vector_length){
    double** jacobi_matrix = allocateMem(vec_length, n);
    if (jacobi_matrix == NULL)
        error_occurred();
    for(int i = 0; i < n; i++)
        for(int j = 0; j < vector_length; j++)
        {
            jacobi_matrix[j][i] = eigensArray.vector[i][j];
        }
    return jacobi_matrix;
}


void print_jacobi(struct eigens* eigensArray, int n, int vector_length){
    for(int i = 0; i < n; i++)
    {
        printf("%f ", eigensArray.value[i]);
    }
    
    jacobi_matrix = jacobi_mat_for_print(eigensArray, n, vector_length);

    for(int i = 0; i < n; i++)
        for(int j = 0; j < vector_length; j++)
        {
            printf("%f ", jacobi_matrix[i][j]);
        }
}

// function that print the given matrix.
void print_matrix(double** matrix, int n, int vector_length){
    for(i = 0; i < n; i++){
        for(j = 0; j < d; j++) {
            if (mat[i][j] < 0 && mat[i][j] > -0.00005) {
                printf("0.0000");
            } else {
                printf("%.4f", mat[i][j]);
            }
            if(j != d - 1){
                printf(",");
            }
        }
        printf("\n");
    }
}

//this function allocates memory if succeeded, and if not returns
double **allocateMem(int n, int d)
{ 
    int i;
    double **matrix;
    matrix = (double **)calloc(n, sizeof(double *));
    if (matrix == NULL)
        return NULL;

    for (i = 0; i < n; i++)
    {
        matrix[i] = (double *)calloc(d, sizeof(double));
        if (matrix[i] == NULL)
        {
            freeLoop(matrix, i);
            free(matrix); 
            return NULL;
        }
    }
    return ans;
}

void error_occurred()
{
    printf("An Error Has Occurred");
    exit(1);
}

int** largest_indexes(double** vectors_matrix, int n)
{
int eigenvalues[2] = {0, 0};
double largest = 0;
int lar_i = -1;
int lar_j = -1;
for (int i = 0; i < n; i++)
{
    for (int j = i + 1; j < n; j++)
    {
        if (abs(vectors_matrix[i][j]) > abs(largest))
        {
            //3. find the Pivot A_ij
            largest = vectors_matrix[i][j];
            lar_i = i;
            lar_j = j;
        }
    }
}
eigenvalues[0] = lar_i;
eigenvalues[1] = lar_j;
return eigenvalues;
}

//this function is for obtain c and s
double* get_c_and_s(double** vectors_matrix, int lar_i, int lar_j)
{
double cs_arr[2] = {0,0}; 

theta = (vectors_matrix[lar_j][lar_j] - vectors_matrix[lar_i][lar_i])/(2*vectors_matrix[lar_i][lar_j]);
t = 1/(abs(theta) + pow((pow(theta, 2) + 1) , 0.5));
if (theta < 0)
    t = t * -1;
c = 1/pow((pow(t, 2) + 1) , 0.5);
s = t * c;

cs_arr[0] = c;
cs_arr[1] = s;
return cs_arr;
}

//this function construct P and initialize
void p_mat_maker(double** vectors_matrix, double** p_matrix, int lar_i, int lar_j)
{
for (i = 0; i < n; i++){
    for (j = 0; j < n; j++){
        if (i == j){
            p_matrix[i][j] = 1;
        }
    }
}

p_matrix[lar_i][lar_i] = c;
p_matrix[lar_j][lar_j] = c;

if (lar_i<lar_j){
    p_matrix[lar_i][lar_j] = s;
    p_matrix[lar_j][lar_i] = -s;
    }

else{
    p_matrix[lar_i][lar_j] = -s;
    p_matrix[lar_j][lar_i] = s;
    }
}

//this function copy vectors_matrix to to_copy which will become to A'
void set_tocopy(double** vectors_matrix, double** to_copy)
{
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            to_copy[i][j] = vectors_matrix[i][j];
        }
    }
}

int first_p_to_v(double** vectors_matrix, double** p_matrix, double** v_matrix,double** temp, int first)
{
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
    return first;
}

//this function calculate A' (according to step 6)
void A_tag_calc(double** to_copy,double** vectors_matrix, int lar_i, int lar_j)
{
    //calculate cell i,i
    to_copy[lar_i][lar_i] = pow(c, 2) * vectors_matrix[lar_i][lar_i];
    to_copy[lar_i][lar_i] += pow(s, 2) * vectors_matrix[lar_j][lar_j];
    to_copy[lar_i][lar_i] += - 2*s*c*vectors_matrix[lar_i][lar_j]+pow(c, 2);

    //calculate cell j,j
    to_copy[lar_j][lar_j] = pow(s, 2) * vectors_matrix[lar_i][lar_i]
    to_copy[lar_j][lar_j] += pow(c, 2) * vectors_matrix[lar_j][lar_j]
    to_copy[lar_j][lar_j] +=  2*s*c*vectors_matrix[lar_i][lar_j] + pow(c, 2);

    //set cells i,j and j,i
    to_copy[lar_i][lar_j] = 0;
    to_copy[lar_j][lar_i] = 0;

    for (int r = 0; r < n; r++)
        for (int m = r + 1; m < n; m++)
        {
            if ( r != lar_i && r != lar_j)
            {
                if (m == lar_i)
                {
                    to_copy[r][m] = c * vectors_matrix[r][lar_i] - s * vectors_matrix[r][lar_j];
                    to_copy[m][r] = c * vectors_matrix[r][lar_i] - s * vectors_matrix[r][lar_j];
                }
                if (m == lar_j)
                {
                    to_copy[r][m] = c * vectors_matrix[r][lar_j] - s * vectors_matrix[r][lar_i];
                    to_copy[m][r] = c * vectors_matrix[r][lar_j] - s * vectors_matrix[r][lar_i];
                }
            }
            to_copy[i][j] = vectors_matrix[i][j]
        }

}

//this function copy the values of to_copy to original array( A = A')
void A_to_A_tag(double** vectors_matrix,double** to_copy);
{
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            vectors_matrix[i][j] = to_copy[i][j];
        }
    }
}

//K-means Functions
double *getFinalCentroids(double *centroids, double *elements, int k, int d, int n, int max_iter, double epsilone)
{
    int bit=1;
    int i;
    int iteration_number=0;
    int min_index = 0;
    int flag = 0;
    double sum = 0.0;
    double min = 0.0;
    int* elements_location;
    int* items_number_clusters;
    double** old_centroids;

    elements_location = (int*)calloc(n, sizeof(int));
    if (!elements_location)
        error_occurred();

    items_number_clusters = (int*)calloc(k, sizeof(int));
    if (!items_number_clusters)
        error_occurred();

    old_centeroids = (double**)calloc(k, sizeof(double*));
    if (!old_centeroids)
        error_occurred();

    for (i=0;i<k;i++){
        old_centeroids[i]=(double*) calloc(d,sizeof(double));
            if (!old_centeroids[i])
                error_occurred();
    }

    while (bit==1 && max_iter>iteration_number)
    {
        iteration_number++;
        initClusters(elements_location, items_number_clusters, n, k);
        assignCentroids(elements, centroids, items_number_clusters,elements_location,k,d,n);
        saveCentroids(old_centroids, centroids, k, d);
        resetCentroids(centroids, k, d);
        updateCentroids(centroids,elements,items_number_clusters,elements_location,d,n,k);
        bit = Convergence(old_centroids, centroids, k, d, epsilone)
    }

    free(elements_location);
    free(items_number_clusters);
    for(i=0; i<k; i++)
        free(old_centeroids[i]);
    free(old_centroids);

    return centroids;
}

void resetCentroids(double* centroids,int k, int d){
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

void saveCentroids(double* old_centroids, double* centroids, int k, int d){
    int i;
    int j;
    for (i=0;i<k;i++){
            for (j=0;j<d;j++)
                old_centeroids[i][j] = centeroids[i][j];
        }
}


void assignCentroids(double *ele,double *cntrds,int* in_clstrs,int* ele_loc,int k,int d,int n)
{
    int i;
    int j;
    int l;
    double sum=0.0;
    int flag;
    int min_index;
    for (i=0; i < n; i++){
        min_index = 0;
        flag = 0;
        for (j=0; j < k; j++){
            sum = 0.0;
                for (l = 0; l < d; l++)
                    sum += pow((elements[i][l] - centeroids[j][l]),2);
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

void updateCentroids(double* cntrds,double* ele,int* in_clstrs,int *ele_loc,int d,int n,int k)
{
    int m;
    int i;
    int q;
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

int Convergence(double* old_centroids,double* centroids,int k,int d,int eps){
    int bit = 0;
    int i;
    double sum;
    for(i=0; i<k; i++){
        sum = 0.0;
        for (j = 0; j<d; j++)
            sum += pow((old_centeroids[i][j] + centeroids[i][j]), 2);
        sum = pow(sum, 0.5);
        if (sum >= eps)
            bit = 1;
    }
    return bit;
}