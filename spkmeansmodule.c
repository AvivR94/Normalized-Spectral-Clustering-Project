#define PY_SSIZE_T_CLEAN
#include "spkmeans.h"
#include <Python.h>


//function built for python usage
static PyObject* python_get_matrix(PyObject *self, PyObject *args){
    PyObject *python_matrix;
    int n;
    int k;
    int vector_length;
    char *func;
    
    if (!PyArg_ParseTuple(args, "Oiiis", &python_matrix, &k, &n, &d, &func))
        return NULL;
    
    matrix = (double**) calloc(n, sizeof(double*));
    assert(matrix != NULL);
    
    for(i = 0; i < n; i++){
        for (j = 0; j < vector_length ; j++) {
            PyObject* num = PyList_GetItem(python_matrix, (vector_length*i) + j);
            matrix[i][j] = PyFloat_AsDouble(num);
        }
    }

    PyObject *lst = PyList_New(0);
    if(!lst)
        return NULL;
    
    //check which goal to choose

    if(strcmp(goal, "wam") == 0){
        wam(matrix, n, vec_length);
    }
    else if (strcmp(goal, "ddg") == 0){
        ddg(matrix, n, vec_length);
    }
    else if (strcmp(goal, "lnorm") == 0){
        lnorm(matrix, n, vec_length);
    }
    else if (strcmp(goal, "jacobi") == 0){
        funcJacobi(matrix, n, vec_length);
    }
    else if (strcmp(goal, "spk") == 0){
        w_matrix = wam_calc(vectors_matrix, n, vec_length);
        d_matrix = ddg_calc(w_matrix, n, vec_length);
        l_matrix = lnorm(calc_d_matrix, n, vec_length);
        eigens_arr = jacobi(l_matrix, n, vec_length);

        if(k == 0){
            k = eigengap_heuristic(l_matrix, n, vec_length);
        } 

        matrix_t_eigen = build_matrix_t_eigen(eigensArray, n, k);

        for(i = 0; i < n; i++){
            for(j = 0; j < k; j++){
                PyObject* num = PyFloat_FromDouble((double)matrix_t_eigen[i][j]);
                if (!num) {
                    Py_DECREF(lst);
                    return NULL;
                }
                PyList_Append(lst, num);
            }
        }

        PyObject *num = PyFloat_FromDouble((double) k);
            if (!num) {
                Py_DECREF(lst);
                return NULL;
            }

        PyList_Append(lst, num);

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
    return lst;
}

 //Python C-API functions

static PyMethodDef myspkmeansmoduleMethods[] = 
{
        {"kmeanspp", (PyCFunction) kmeanspp, METH_VARARGS, PyDoc_STR("Kmeans++")},
        {"pythonGetMatR", (PyCFunction) pythonGetMatR, METH_VARARGS, PyDoc_STR("c usage for getting the matrixes")},
        {NULL,  NULL, 0, NULL}
};

static struct PyModuleDef _moduledef = 
{
        PyModuleDef_HEAD_INIT,
        "spkmeansmodule",
        NULL,
        -1,
        myspkmeansmoduleMethods
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

