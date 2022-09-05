#define PY_SSIZE_T_CLEAN
#include "spkmeans.h"
#include <Python.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* helper function */
static PyObject* makeListFromMatrix(double** matrix, int n, int k){
    int i, j;
    PyObject *lst;
    PyObject *row;

    lst = PyList_New(0);
    if(!lst)
        return NULL;
    for(i = 0; i < n; i++){
        row = PyList_New(0);
        if(!row){
            Py_DECREF(lst);
            return NULL;
        }
        for(j = 0; j < k; j++){
            PyList_Append(row,PyFloat_FromDouble((double)matrix[i][j]));
        }
        PyList_Append(lst, row);
    }
    return lst;      
}

/**
 * function called from python.
 * performs the spkmeans by C code,
 * and returns the matrix needed for python.
 */
static PyObject* getMatrixByGoal(PyObject *self, PyObject *args){
    PyObject *python_matrix;
    PyObject *final_mat;
    PyObject *eigens_list;
    PyObject *row;
    int n, k, i, j;
    int vector_length;
    char* goal;
    double** original_matrix;
    double** w_matrix;
    double** d_matrix;
    double** l_matrix;
    double** t_matrix;
    struct eigens* eigens_arr;
    if (!PyArg_ParseTuple(args, "iiiOs",&k, &n, &vector_length, &python_matrix, &goal)){
        errorOccured();}

    original_matrix = allocateMem(n, vector_length);
    if (original_matrix == NULL)
        errorOccured();
    for(i = 0; i < n; i++){
        for (j = 0; j < vector_length ; j++) {
            PyObject* num = PyList_GetItem(python_matrix, (vector_length*i) + j);
            original_matrix[i][j] = PyFloat_AsDouble(num);
        }
    }
    final_mat = PyList_New(0);
    
    /* check which goal to choose */

    if(strcmp(goal, "wam") == 0){
        w_matrix = wamCalc(original_matrix, n, vector_length);
        final_mat = makeListFromMatrix(w_matrix, n, vector_length);
        for(i = 0; i < n ; i++){
            free(w_matrix[i]);
        }
        free(w_matrix);
    }

    else if (strcmp(goal, "ddg") == 0){
        w_matrix = wamCalc(original_matrix, n, vector_length);
        d_matrix = ddgCalc(w_matrix, n);
        final_mat = makeListFromMatrix(d_matrix, n, vector_length);
        for(i = 0; i < n ; i++){
            free(w_matrix[i]);
            free(d_matrix[i]);
        }
        free(w_matrix);
        free(d_matrix);
    }

    else if (strcmp(goal, "lnorm") == 0){
        w_matrix = wamCalc(original_matrix, n, vector_length);
        d_matrix = ddgCalc(w_matrix, n);
        l_matrix = lnormCalc(w_matrix, d_matrix, n, vector_length);
        final_mat = makeListFromMatrix(l_matrix, n, vector_length);
        for(i = 0; i < n ; i++){
            free(w_matrix[i]);
            free(d_matrix[i]);
            free(l_matrix[i]);
        }
        free(w_matrix);
        free(d_matrix);
        free(l_matrix);
    }

    else if (strcmp(goal, "jacobi") == 0){
        double** jacobi_matrix;
        eigens_arr = jacobiCalc(original_matrix, n, 0, 0);
        jacobi_matrix = jacobiMatForPrint(eigens_arr, n);
        eigens_list = PyList_New(0);
        if(!eigens_list)
            errorOccured();
        for(i=0; i<n; i++){
            PyList_Append(eigens_list, PyFloat_FromDouble((double)eigens_arr[i].value));
        }
        PyList_Append(final_mat, eigens_list);
        for(i = 0; i < n; i++){
            row = PyList_New(0);
            if(!row){
                Py_DECREF(final_mat);
                errorOccured();
            }
            for(j = 0; j < k; j++){
                PyList_Append(row,PyFloat_FromDouble((double)jacobi_matrix[i][j]));
            }
            PyList_Append(final_mat, row);
        }
        for(i = 0; i < n ; i++){
            free(eigens_arr[i].vector);
            free(jacobi_matrix[i]);
        }
        free(eigens_arr);
        free(jacobi_matrix);
    }

    else if (strcmp(goal, "spk") == 0){
        w_matrix = wamCalc(original_matrix, n, vector_length);
        d_matrix = ddgCalc(w_matrix, n);
        l_matrix = lnormCalc(w_matrix, d_matrix, n, vector_length);
        eigens_arr = jacobiCalc(l_matrix, n, 0, 1);
        if(k == 0){
            k = eigengapHeuristic(eigens_arr, n);
        } 
        t_matrix = createMatrixT(eigens_arr, n, k);
        final_mat = makeListFromMatrix(t_matrix, n, k);
        for(i = 0; i < n ; i++){
            free(w_matrix[i]);
            free(d_matrix[i]);
            free(l_matrix[i]);
            free(eigens_arr[i].vector);
            free(t_matrix[i]);
        }
        free(w_matrix);
        free(d_matrix);
        free(l_matrix);
        free(eigens_arr);
        free(t_matrix);
    }

    for(i=0; i<n; i++) {
        free(original_matrix[i]);
    }
    free(original_matrix);
    return Py_BuildValue("O", final_mat);
}

/**
 * function called from python.
 * performs the kmeans++ algorithm,
 * returns the centroids of the clusters.
 */
static PyObject* fit(PyObject *self, PyObject *args){
    PyObject *items_of_py;
    PyObject *centroids_arr;
    PyObject *centroids_from_py;
    PyObject *vectors_from_py;
    int n, k, d; 
    int i, j;
    double** elements;
    double** centroids;

    if (!PyArg_ParseTuple(args, "iiiOO", &k, &n, &d, &centroids_from_py, &vectors_from_py)){
        errorOccured();
    }

    elements = (double**)calloc(n, sizeof(double*));
    if (!elements)
        errorOccured();

    for (i=0;i<n;i++){
        elements[i]=(double*) calloc(d,sizeof(double));
            if (!elements[i])
                errorOccured();
    }

    for (i = 0; i < n; i++) {
        for (j = 0; j < d; j++) {
            items_of_py = PyList_GetItem(vectors_from_py, i * d + j);
            elements[i][j] = PyFloat_AsDouble(items_of_py);
        }
    }

    /*init centroids*/
    centroids = (double**)calloc(k, sizeof(double*));
    if (!centroids)
        errorOccured();

    for (i=0;i<k;i++){
        centroids[i]=(double*) calloc(d,sizeof(double));
            if (!centroids[i])
                errorOccured();
    }

     for (i = 0; i < k; i++) {
        for (j = 0; j < d; j++) {
            items_of_py = PyList_GetItem(centroids_from_py, i * d + j);
            centroids[i][j] = PyFloat_AsDouble(items_of_py);
        }
    }
    /*now we have first centroids from python in C array*/
    getFinalCentroids(centroids, elements, k, d, n, 300, 0);
    /*creating final centroids array for python in matrix form*/
    centroids_arr = makeListFromMatrix(centroids,k,k);
    for(i=0; i<n; i++)
        free(elements[i]);
    free(elements);
    for(i=0; i<k; i++)
        free(centroids[i]);
    free(centroids);
    return Py_BuildValue("O", centroids_arr);
}


 /* Python C-API functions */
static PyMethodDef spkmeansmoduleMethods[] = {
        {"getMatrixByGoal", (PyCFunction) getMatrixByGoal, METH_VARARGS, PyDoc_STR
        ("c-api for python, retrieving the matrix accoring to the goal")},
        {"fit", (PyCFunction) fit, METH_VARARGS, PyDoc_STR("Kmeans++")},
        
        {NULL,  NULL, 0, NULL}
};

static struct PyModuleDef _moduledef = {
        PyModuleDef_HEAD_INIT,
        "spkmeansmodule",
        NULL,
        -1,
        spkmeansmoduleMethods
};

PyMODINIT_FUNC PyInit_spkmeansmodule(void)
{
    PyObject *m;
    m = PyModule_Create(&_moduledef);
    if(!m){
        return NULL;
    }
    return m;
}
