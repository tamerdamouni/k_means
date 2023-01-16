

#ifndef FPROJ_SPKMEANS_H
#define FPROJ_SPKMEANS_H


int K;
double **centroids; /*defined as a global variable in order to be able to free the space it occupied */
/**
 * gets text file and converts it to a matrix
 * @param input_file
 * @param num_Dim_points
 * @return  points matrix
 *
 */
double** get_points(char* input_file,int* num_Dim_points);
/**
 * takes the matrix all_points and runs the corect function acording to the given goal on the matrix
 * @param all_points
 * @param goal
 * @param numPoints
 * @param points_dimension
 * @return the matrix that you get from runing the function of the goal
 *
 */
double **specified_execution(double** all_points, int goal, int numPoints, int points_dimension);
/**
 * calculates Weighted Adjacency Matrix matrix
 * @param mat
 * @param points
 * @param numPoints
 * @param points_dimension
 * @return wam matrix
 *
 */
double** Weighted_Adjacency_Matrix(double** mat, double** points,int numPoints, int points_dimension);
/**
 * calculates Diagonal Degree Matrix matrix
 * @param mat
 * @param numPoints
 * @return DDG
 */
double**  Diagonal_Degree_Matrix(double** mat,int numPoints);
/**
 * calculates Laplacian matrix
 * @param W
 * @param D
 * @param n
 * @return
 */
double** Build_Laplacian(double** W, double **D, int n);

double** Perform_Jacobi(double** lnorm,int n,double* eigenvalues);
/**
 * calculates K acording to the algorithm of Eigengap Heuristic and the orignal order of the eigenvectors is saved in the indexes array
 * @param eigenvalues
 * @param n
 * @param indexes
 * @return K
 */
int The_Eigengap_Heuristic(double *eigenvalues, int n, int *indexes);
/**
 * normal mergesort that sorts a double and an int array simultaneously
 * @param arr_double
 * @param len
 * @param arr_int
 */
void mergeSort(double *arr_double, int len,int* arr_int);
double* slice_double(double *arr, int start, int end);
int* slice_int(int *arr, int start, int end);
void merge(double *result, double *left_double, double *right_double, int leftLen, int rightLen,int* ind,int* left_int,int* right_int);
/**
 * builds a new matrix with row_num rows and cul_num colmns
 * @param rows_num
 * @param cul_num
 * @return new matrix
 */
double** Matrix_init(int rows_num , int cul_num);
double norm2(double* point1,double* point2,int points_dimension);
double** calc_U(int n, double** V, int k, int *indexes);
double** calc_T(int n, int k, double** U);
void print_matrix(double **M, int row_number, int col_number);
void Matrix_multiplication(double **M1, double **M2, int n, double **answer);
double* divide(double* vector,int num,int vector_dim);
double diff_norm_pow2(double* vector1, double* vector2, int vector_dim);
int calc_norm(double** points, double* point, int k, int vector_dim);
void freeMatrix(double** mat, int len);
/**
 * checks if matrix is symetric
 * @param mat
 * @param rows
 * @param cols
 * @return
 */
int is_sym(double **mat,int rows, int cols);
/**
 * returns (n+1)*n matrix where the first row is the eigenvalues and the rest is V.
 * @param V
 * @param eigenvalues
 * @param numPoints
 * @return
 */
double **row_concatenation(double **V, double *eigenvalues, int numPoints);
void execute(int k, int max_iter, int d, int numPoints, double **all_points, double **init_centroids);
double** matT(double **all_points, int k, int numPoints, int points_dimension);
#endif /*FPROJ_SPKMEANS_H*/
