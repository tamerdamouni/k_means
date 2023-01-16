#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <Python.h>
#include "spkmeans.h"

void convert_object_to_matrix(PyObject *all_points, int d, int numPoints, double **all_points_matrix);
void convert_matrix_to_object(double **matrix, int rows, int cols, PyObject *py_result);

/**
 * this function is called when the goal is spk and returns the centroids
 * @param k
 * @param max_iter
 * @param d
 * @param all_points
 * @param init_centroids
 * @return
 */
static PyObject *finalspk(int k , int max_iter, int d, PyObject *all_points, PyObject *init_centroids){
    double **vectors_matrix, **initial_centroids;
    PyObject  *py_result, *tmp;
    Py_ssize_t i, j;
    int num_of_vectors;
    num_of_vectors = (int) PyObject_Length(all_points);
    vectors_matrix = (double**)malloc( num_of_vectors*sizeof(double*));
    if(vectors_matrix == NULL){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    initial_centroids = (double**) malloc(k*sizeof (double*));
    if(initial_centroids == NULL){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    convert_object_to_matrix(all_points,d,num_of_vectors,vectors_matrix);
    convert_object_to_matrix(init_centroids,k,d,initial_centroids);
    execute(k, max_iter, d, num_of_vectors, vectors_matrix, initial_centroids);
    py_result = PyList_New(k);
    for(i=0;i<k;i++){
        tmp = PyList_New(d);
        for(j=0;j<d;j++)
        {
            PyList_SetItem(tmp,j,PyFloat_FromDouble(centroids[i][j]));
        }
        PyList_SetItem(py_result,i,tmp);
    }
    for(i=0;i<k;i++){free(centroids[i]);}
    free(centroids);
    return py_result;
}

/**
 * the api for the above function
 * @param self
 * @param args
 * @return
 */
static PyObject *fitfinalspk_api(PyObject *self,PyObject *args){
    int k, max_iter,d;
    PyObject *all_points;
    PyObject *init_centroids;
    if(!PyArg_ParseTuple(args,"iiiOO",&k,&max_iter,&d,&all_points,&init_centroids));
    return Py_BuildValue("O", finalspk(k,max_iter,d,all_points,init_centroids));
}

/**
 * this matrix gets all points and calculates T  then returns it to python.
 * it is called in case goal == spk
 * @param all_points
 * @param k
 * @param numPoints
 * @param d
 * @return
 */
static PyObject *fitgetTmat(PyObject *all_points, int k, int numPoints, int d ){
    double **all_points_matrix;
    PyObject *py_result, *tmp;
    int i,j;
    all_points_matrix = (double**) malloc(numPoints * sizeof (double*));
    if(all_points_matrix == NULL){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for(i=0;i<numPoints;i++){
        all_points_matrix[i] = (double*) malloc(d*sizeof (double));
        if(all_points_matrix[i] == NULL){
            printf("An Error Has Occurred\n");
            exit(1);
        }
    }
    convert_object_to_matrix(all_points,d,numPoints,all_points_matrix);
    all_points_matrix = matT(all_points_matrix,k,numPoints,d);
    /* matT returns a matrix of dimensions nXk  */
    if(k==0) k=K;
    if(K==1){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    py_result = PyList_New(numPoints);
    for(i=0;i<numPoints;i++){
        tmp = PyList_New(k);
        for(j=0;j<k;j++)
        {
            PyList_SetItem(tmp,j,PyFloat_FromDouble(all_points_matrix[i][j]));
        }
        PyList_SetItem(py_result,i,tmp);
    }
    return py_result;
}
/**
 * the api of the above function
 * @param self
 * @param args
 * @return
 */
static PyObject *fitgetTmat_api(PyObject *self,PyObject *args){
    int k, numPoints,d;
    PyObject *all_points;
    if(!PyArg_ParseTuple(args,"Oiii",&all_points,&k,&numPoints,&d));
    return Py_BuildValue("O", fitgetTmat(all_points,k,numPoints,d));
}
/**
 * this function is called when any goal is given as an input.
 * it returns the result to python which prints it.
 * @param all_points
 * @param goal
 * @param numPoints
 * @param d
 * @return
 */
static PyObject *fit(PyObject *all_points, int  goal, int numPoints, int d){
    double **all_points_matrix, **result;
    PyObject *py_result;
    all_points_matrix = (double**) malloc(numPoints * sizeof (double*));
    if (all_points_matrix == NULL){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    convert_object_to_matrix(all_points,d,numPoints,all_points_matrix);
    result = specified_execution(all_points_matrix,goal,numPoints,d);
    if(goal == 3 ){
        py_result = PyList_New(numPoints+1);
        convert_matrix_to_object(result,numPoints+1,numPoints,py_result);
        freeMatrix(result,numPoints+1);
    }
    else {
        py_result = PyList_New(numPoints);
        convert_matrix_to_object(result,numPoints,numPoints,py_result);
        freeMatrix(result,numPoints);
        freeMatrix(all_points_matrix,numPoints);
    }
    return py_result;
}

/**
 * THE API for the above function
 * @param self
 * @param args
 * @return
 */
static PyObject *fit_api(PyObject *self,PyObject *args){
    int numPoints, d, goal;
    PyObject *all_points;
    if(!PyArg_ParseTuple(args,"Oiii",&all_points,&goal,&numPoints,&d));
    return Py_BuildValue("O",fit(all_points,goal,numPoints,d));
}


static PyMethodDef fit_apiMethods[]={
        {"fit",(PyCFunction)fit_api, METH_VARARGS,
                        PyDoc_STR("kmeans spk")},
        {"fitgetTmat",(PyCFunction)fitgetTmat_api, METH_VARARGS,
                        PyDoc_STR("kmeans spk")},
        {"finalspk",(PyCFunction)fitfinalspk_api, METH_VARARGS,
                        PyDoc_STR("kmeans spk")},
        {NULL, NULL, 0, NULL}

};


static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "kmeansspk",
        NULL,
        -1,
        fit_apiMethods
};

PyMODINIT_FUNC
PyInit_kmeansspk(void){
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}

/**
 *
 * @param all_points
 * @param d
 * @param numPoints
 * @param all_points_matrix
 */
void convert_object_to_matrix(PyObject *all_points, int d, int numPoints, double **all_points_matrix){
    int i,j;
    PyObject  *py_curr_vector;
    for (i=0;i<numPoints;i++){
        all_points_matrix[i] = (double*) malloc(d * sizeof(double ));
        if(all_points_matrix[i] == NULL){
            printf("An Error Has Occurred\n");
            exit(1);
        }
    }
    for(i=0;i<numPoints;i++){
        py_curr_vector = PyList_GetItem(all_points, i);
        for (j=0;j<d;j++){
            all_points_matrix[i][j] = PyFloat_AsDouble(PyList_GetItem(py_curr_vector, j));
        }
    }
}
/**
 *
 * @param matrix
 * @param rows
 * @param cols
 * @param py_result
 */
void convert_matrix_to_object(double **matrix, int rows, int cols, PyObject *py_result){
    int i , j;
    PyObject *tmp;
    for(i=0;i<rows;i++){
        tmp = PyList_New(cols);
        for(j=0;j<cols;j++)
        {
            PyList_SetItem(tmp,j,PyFloat_FromDouble(matrix[i][j]));
        }
        PyList_SetItem(py_result,i,tmp);
    }
}