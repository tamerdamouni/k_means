
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include "spkmeans.h"


int main(int argc,char **argv){
    int* num_Dim_points;
    double **all_points, **result;
    int goal, num_points,dim_points, i;
    if (argc != 3){
        printf("Invalid Input!");
        return 0;
    }
    if(strcmp(argv[1],"wam")==0) goal = 0;
    else if(strcmp(argv[1],"ddg")==0) goal = 1;
    else if(strcmp(argv[1],"lnorm")==0) goal = 2;
    else if(strcmp(argv[1],"jacobi")==0) goal = 3;
    else{
        printf("Invalid Input!");
        return 0;
    }
    num_Dim_points=(int*) malloc(sizeof(int)*2);
    all_points = get_points(argv[2],num_Dim_points);
    num_points=num_Dim_points[0];
    dim_points=num_Dim_points[1];
    free(num_Dim_points);
    result = specified_execution(all_points, goal, num_points, dim_points);
    if(goal == 3){
        for(i=0;i<dim_points;i++){
            if(result[0][i]<0 && result[0][i] > - 0.00005){ result[0][i]=0 ; }
        }
        print_matrix(result,num_points+1,dim_points);
        freeMatrix(result,num_points+1);
    }
    else {
        print_matrix(result,num_points,num_points);
        freeMatrix(result,num_points);
        freeMatrix(all_points,num_points);
    }
    return 1;
}



double** get_points(char* input_file,int* num_Dim_points){
    int all_vectors_num, sum_of_vector_dims, vector_dim, i, j;
    char c;
    double charc;
    FILE *file;
    double** all_points;
    file = fopen(input_file, "r");
    if (file == NULL) exit(EXIT_FAILURE);
    all_vectors_num=0;
    sum_of_vector_dims=0;
    while((c= fgetc(file))!=EOF){
        if(c=='\n'){
            all_vectors_num++;
            sum_of_vector_dims++;

        }
        if(c==','){
            sum_of_vector_dims++;
        }
    }
    rewind(file);/* go back to the beginning of the file */
    vector_dim = sum_of_vector_dims/all_vectors_num; /*find each vector dimension*/
    num_Dim_points[1] = vector_dim;
    all_points = (double**) malloc(all_vectors_num*sizeof (double*)); /*build a 2d array with all the vectors*/
    if(all_points == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i=0;i< all_vectors_num;i++ ){
        all_points[i] = (double*) malloc(vector_dim*sizeof (double));
        if(all_points[i] == NULL) {
            printf("An Error Has Occurred\n");
            exit(1);
        }
    }
    /*read file => */
    for (i=0;i<all_vectors_num;i++){
        for(j=0;j<vector_dim;j++){
            fscanf(file,"%lf",&charc);
            all_points[i][j] = (double) charc;
            fgetc(file);
        }
    }
    num_Dim_points[0]=all_vectors_num;
    fclose(file);
    return all_points;
}


double **specified_execution(double **all_points, int goal, int numPoints, int points_dimension) {
    double **mat, **W, **D, **L, **V, **result;
    double *eigenvalues;
    if( goal == 3){
        if (!is_sym(all_points,numPoints,points_dimension)){
            freeMatrix(all_points,numPoints);
            printf("Invalid Input!");
            exit(1);
        }
        eigenvalues = (double *) malloc(numPoints* sizeof (double));
        if(eigenvalues==NULL){
            printf("An Error Has Occurred\n");
            exit(1);
        }
        V = Perform_Jacobi(all_points,numPoints,eigenvalues);
        result = row_concatenation(V,eigenvalues,numPoints);
        freeMatrix(V,numPoints);
        free(eigenvalues);
        return result;
    }
    mat = Matrix_init(numPoints,numPoints);
    W = Weighted_Adjacency_Matrix(mat, all_points, numPoints, points_dimension);
    if(goal == 0) {
        return W;
    }
    D = Diagonal_Degree_Matrix(W,numPoints);
    if(goal == 1 ){
        freeMatrix(W,numPoints);
        return D;
    }
    L = Build_Laplacian(W,D,numPoints);
    if (goal==2) {
        freeMatrix(W,numPoints);
        freeMatrix(D,numPoints);
        return L;
    }
    else {
        printf("An Error Has Occurred\n");
        exit(1);
    }
}



double** Weighted_Adjacency_Matrix(double** mat, double** points,int numPoints,int points_dimension ){
    int i,j;
    if(mat == NULL || points == NULL){
        exit(1);
    }
    for (i=0; i<numPoints; i++){
        for (j= 0; j < numPoints ; j++) {
            if (i==j){
                continue;
            }
            mat[i][j]=exp(-((norm2(points[i],points[j], points_dimension)/2)));
        }
    }
    return mat;
}



double**  Diagonal_Degree_Matrix(double **mat,int numPoints){
    int i,j;
    double sum=0;
    double **DDM = Matrix_init(numPoints, numPoints);
    for (i= 0; i <numPoints ; i++) {
        for (j=0;j<numPoints;j++){
            sum+=mat[i][j];
        }
        DDM[i][i]=sum;
        sum=0;
    }
    return DDM;
}

double **Build_Laplacian(double** W, double **D, int n){
    double **L;
    double **I;
    double **DWD,**temp1,**temp2;
    int i,j ;
    I = Matrix_init(n, n);
    L = Matrix_init(n, n);
    /*building I + calculating D^-0.5*/
    for(i=0;i<n;i++){
        I[i][i] = 1;
        D[i][i] = pow(D[i][i],-0.5);
    }

    temp1=Matrix_init(n,n);
    temp2=Matrix_init(n,n);

    /*by now => D actually represents D^-0.5 and now we calculate DWD*/
    Matrix_multiplication(D,W,n,temp1);
    Matrix_multiplication(temp1,D,n,temp2);

    DWD=temp2;
    for (i=0;i<n;i++){
        for(j=0;j<n;j++){
            L[i][j] = I[i][j]-DWD[i][j];
        }
    }
    freeMatrix(temp1,n);
    freeMatrix(temp2,n);
    freeMatrix(I,n);
    return L;
}

void Build_Rotating_Matrix(double **A, double**P, double **PT, int n){
    int i, j, sign_theta;
    int max_row, max_col;
    double max, c, theta, t, s;
    max = A[0][1];
    max_row=0;
    max_col=1;
    for(i=0;i<n;i++) {
        for (j = 0; j < n; j++) {
            if (fabs(A[i][j]) > fabs(max) && i != j) {
                max = A[i][j];
                max_row = i;
                max_col = j;
            }
            if (i == j) {
                P[i][j] = 1;
                PT[i][j] = 1;
            }
        }
    }
    theta = (A[max_col][max_col]-A[max_row][max_row])/(2*A[max_row][max_col]);
    if(theta>=0) sign_theta=1;
    else sign_theta=-1;

    t = sign_theta/ (fabs(theta)+ sqrt(pow(theta,2)+1));
    c = 1/sqrt(pow(t,2)+1);
    s = t*c;

    P[max_row][max_row] = c;
    P[max_row][max_col]=s;
    P[max_col][max_col]=c;
    P[max_col][max_row]=-s;

    PT[max_row][max_row] = c;
    PT[max_row][max_col]=-s;
    PT[max_col][max_col]=c;
    PT[max_col][max_row]=s;


}



void Rotate(double **A_PRIME, double **PT, double **A, double **P, int n){
    double **temp1;
    temp1 = Matrix_init(n,n);

    Matrix_multiplication(PT,A,n,temp1);
    Matrix_multiplication(temp1,P,n,A_PRIME);
    freeMatrix(temp1,n);
}

int check_convergence (double **A, double **A_PRIME,int n){
    int i,j;
    double sumA,sumA_PRIME;
    sumA_PRIME =0;
    sumA=0;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(i!=j) {
                sumA += pow(A[i][j], 2);
                sumA_PRIME += pow(A_PRIME[i][j],2);
            }
        }
    }
    if(sumA-sumA_PRIME<= pow(10,-5)) return 1;
    return 0;
}
double **Perform_Jacobi(double **lnorm,int n,double* eigenvalues) {
    double **P, **PT, **A, **V, **A_PRIME, **temp;
    int i, j, counter;
    counter = 0;
    A = lnorm;
    V = Matrix_init(n, n);
    for (i = 0; i < n; ++i) {
        V[i][i] = 1;
    }
    P = Matrix_init(n, n);
    PT = Matrix_init(n, n);
    while (counter < 100) {
        for(i=0;i<n;i++){
            for(j=0;j<n;j++){
                P[i][j]=0;
                PT[i][j]=0;
            }
        }
        Build_Rotating_Matrix(A, P, PT, n);
        temp = Matrix_init(n,n);
        Matrix_multiplication(V, P, n,temp);
        freeMatrix(V,n);
        V=temp;
        A_PRIME = Matrix_init(n,n);
        Rotate(A_PRIME,PT, A, P, n);
        if (check_convergence(A, A_PRIME, n)) {
            freeMatrix(A,n);
            A = A_PRIME;
            break;
        }
        freeMatrix(A,n);
        A = A_PRIME;
        counter++;
    }
    for (i = 0; i < n; i++) {
        eigenvalues[i] = A[i][i];
    }
    freeMatrix(A,n);
    freeMatrix(P,n);
    freeMatrix(PT,n);
    return V;
}




/* **************** functions used by python ********************* */


int The_Eigengap_Heuristic(double *eigenvalues, int n, int *indexes) {
    int i, max_ind, k;
    double max_eig, *tmp_eigenvalues;
    tmp_eigenvalues = (double *) malloc(n * sizeof (double));
    for(i=0;i<n;i++){
        tmp_eigenvalues[i] = eigenvalues[i];
    }
    mergeSort(tmp_eigenvalues,n,indexes);
    max_ind=0;
    max_eig = fabs(tmp_eigenvalues[0] - tmp_eigenvalues[1]);
    for (i = 0; i < floor(n / 2); i++) {
        if (fabs(tmp_eigenvalues[i] - tmp_eigenvalues[i + 1]) > max_eig) {
            max_ind = i;
            max_eig = fabs(tmp_eigenvalues[i] - tmp_eigenvalues[i + 1]);
        }
    }
    k = max_ind;
    free(tmp_eigenvalues);
    return k+1;
}

double** matT(double **all_points, int k, int numPoints, int points_dimension){
    double **W, **D, **V, **L, **U, **T;
    double *eigenvalues;
    int *indexes, i;
    indexes = (int *) malloc(numPoints * sizeof(int));
    for (i=0;i<numPoints;i++){
        indexes[i]=i;
    }
    eigenvalues = (double*) malloc(sizeof(double)*numPoints);
    if (eigenvalues == NULL){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    W = Weighted_Adjacency_Matrix(Matrix_init(numPoints,numPoints), all_points,numPoints, points_dimension);
    D = Diagonal_Degree_Matrix(W,numPoints);
    L = Build_Laplacian(W, D, numPoints);
    V = Perform_Jacobi(L, numPoints, eigenvalues);
    if(k==0){
        k = The_Eigengap_Heuristic(eigenvalues, numPoints,indexes);
        K = k;
    }
    else mergeSort(eigenvalues,numPoints,indexes);
    U = calc_U(numPoints, V, k, indexes);
    T = calc_T(numPoints, k, U);
    free(eigenvalues);
    return T;
}

double** calc_U(int n, double** V, int k, int *indexes) {
    int i, j;
    double **U;
    U = Matrix_init(n, k);
    /*creating U */
    for (i = 0; i < n; i++) {
        for (j = 0; j < k; j++) {
            U[i][j] = V[i][indexes[j]];
        }
    }
    free(indexes);
    return U;
}
double** calc_T(int n, int k, double** U){
    int i, j;
    double sum;
    double** T;
    T = Matrix_init(n,k);
    for(i=0;i<n;i++){
        sum = 0;
        for(j=0;j<k;j++){
            sum += pow(U[i][j],2);
        }
        for(j=0;j<k;j++){
            if(sum!=0) {
                T[i][j] = U[i][j] / sqrt(sum);
            }
            if(sum==0){
                T[i][j] = U[i][j];
            }
        }
    }
    return T;
}

void execute(int k , int max_iter, int d, int numPoints, double **all_points, double **init_centroids) {
    int vector_dim, i, j, centroid_index, iteration, all_vectors_num;
    double **sum_mat, **old_centroids;
    int *clusters_size;
    all_vectors_num=numPoints;
    vector_dim = d;
    iteration=0;
    centroids = (double **) malloc(k*sizeof(double*));
    if (centroids == NULL){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i=0;i<k;i++){
        centroids[i] = (double*) malloc(d*sizeof(double));
        if(centroids[i] == NULL){
            printf("An Error Has Occurred\n");
            exit(1);
        }
    }
    for(i=0;i<k;i++){
        for (j=0;j<d;j++) {
            centroids[i][j] = init_centroids[i][j];
        }
    }
    freeMatrix(init_centroids,k);
    sum_mat = (double **) malloc(k * sizeof(double*));
    if (sum_mat == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i = 0; i < k; i++) {
        sum_mat[i] = (double *) malloc(vector_dim * sizeof(double));
        if (sum_mat[i] == NULL) {
            printf("An Error Has Occurred\n");
            exit(1);
        }
    }
    /*initializing clusters sizes to zero*/
    clusters_size = (int *) malloc(k * sizeof(int));
    if (clusters_size == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    /*memorizing the old centroids in order to calc the norm after updating*/
    old_centroids = (double **) malloc(k * sizeof(double *));
    if (old_centroids == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (j = 0; j < k; j++) {
        old_centroids[j] = (double *) malloc(vector_dim * sizeof(double));
        if(old_centroids[j]==NULL){
            printf("An Error Has Occurred\n");
            exit(1);
        }
    }

    /*ALGORITHM LOOP*/
    while (iteration < max_iter) {
        /* reset cluster sizes and sum matrix */
        for (i = 0; i < k; i++) {
            clusters_size[i] = 0;
            for (j = 0; j < vector_dim; j++) {
                sum_mat[i][j] = 0;
            }
        }
        /*calculating new sum matrix with the new clusters' sizes */
        for (j = 0; j < all_vectors_num; j++) {
            centroid_index = calc_norm(centroids, all_points[j], k,vector_dim);/*calc norm returns index of the current closest centroid to the vector*/
            clusters_size[centroid_index] += 1;
            for (i = 0; i < vector_dim; i++) {
                sum_mat[centroid_index][i] += all_points[j][i];
            }
        }
        /*substituting centroids into old centroids before updating centroids*/
        for (i = 0; i < k; i++) {
            for (j = 0; j < vector_dim; j++) {
                old_centroids[i][j] = centroids[i][j];
            }
        }
        /*calculating the new centroids */
        for (i = 0; i < k; i++) {
            free(centroids[i]);
            centroids[i] = divide(sum_mat[i], clusters_size[i], vector_dim);
        }
        iteration += 1;
    }
    /* memory releasing */
    freeMatrix(all_points, all_vectors_num);
    freeMatrix(old_centroids, k);
    freeMatrix(sum_mat,k);
    free(clusters_size);

    /* will free centroids later */
}

/* **************************** AUX ******************************  */

double** Matrix_init(int rows_num , int cul_num){
    double** mat;
    int i,j;
    mat= (double**) malloc(sizeof(double*)*rows_num);
    if(mat == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i=0;i<rows_num;i++){
        mat[i]= (double*) malloc(sizeof(double)*cul_num);
        if(mat[i] == NULL) {
            printf("An Error Has Occurred\n");
            exit(1);
        }
        for (j=0;j<cul_num;j++){
            mat[i][j]=0;
        }
    }
    return mat;
}

double norm2(double* point1,double* point2,int points_dimension){
    double sum;
    int i;
    sum=0;
    for (i=0;i<points_dimension;i++){
        sum+= pow(point1[i]-point2[i],2);
    }
    sum = sqrt(sum);
    return sum;
}

void print_matrix(double **M, int row_number, int col_number){
    int i,j;
    for(i=0;i<row_number;i++){
        for(j=0;j<col_number;j++){
            printf("%.4f",M[i][j]);
            if(j!=col_number-1) printf(",");
        }
        printf("\n");
    }
}

void Matrix_multiplication(double **M1, double **M2, int n, double **answer){
    int i , j, k;
    for (i=0;i<n;i++){
        for(j=0;j<n;j++){
            for(k=0;k<n;k++)
                answer[i][j]+=M1[i][k]*M2[k][j];
        }
    }
}
int calc_norm(double** points, double* point, int k, int vector_dim){
    double closest_norm_value, x;
    int i, closest_norm;
    closest_norm =0;
    closest_norm_value = diff_norm_pow2(points[0], point, vector_dim);
    for(i=1;i<k;i++){
        x = diff_norm_pow2(points[i], point, vector_dim);
        if(x < closest_norm_value){
            closest_norm=i;
            closest_norm_value =x;
        }
    }
    return closest_norm;
}
double diff_norm_pow2(double* vector1, double* vector2, int vector_dim){
    double* arr;
    int i;
    double sum;
    arr = (double*) malloc(vector_dim*sizeof (double));
    if(arr ==NULL){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for(i=0;i<vector_dim;i++){
        arr[i]=vector1[i]-vector2[i];
    }
    sum=0;
    for(i=0;i<vector_dim;i++){
        sum = sum + arr[i]*arr[i];
    }
    free(arr);
    return sum;
}
double* divide(double* vector,int num,int vector_dim){
    double *result;
    int i;
    result = (double*) malloc(vector_dim*sizeof(double ));
    if(result == NULL){
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i = 0; i < vector_dim; i++) {
        result[i]=  vector[i]/num;
    }
    return result;
}

void mergeSort(double *arr_double, int len,int* arr_int)
{
    double *left_double,*right_double;
    int *left_int,*right_int;
    if (len <= 1)
    {
        return;
    }
    left_double = slice_double(arr_double, 0, len / 2 + 1);
    right_double = slice_double(arr_double, len / 2, len);
    left_int=slice_int(arr_int,0,len/2 + 1);
    right_int=slice_int(arr_int,len/2,len);


    mergeSort(left_double, len / 2,left_int);
    mergeSort(right_double, len - (len / 2),right_int);

    merge(arr_double, left_double, right_double, len / 2, len - (len / 2),arr_int,left_int,right_int);
}

void merge(double *result, double *left_double, double *right_double, int leftLen, int rightLen,int* ind,int* left_int,int* right_int)
{
    int i = 0, j = 0;
    while(i < leftLen && j < rightLen)
    {
        if (left_double[i] < right_double[j])
        {
            result[i + j] = left_double[i];
            ind[i + j] = left_int[i];
            i++;
        }
        else
        {
            result[i + j] = right_double[j];
            ind[i + j] = right_int[j];
            j++;
        }
    }

    for(; i < leftLen; i++)
    {
        result[i + j] = left_double[i];
        ind[i + j] = left_int[i];
    }
    for(; j < rightLen; j++)
    {
        result[i + j] = right_double[j];
        ind[i + j] = right_int[j];
    }

    free(left_double);
    free(right_double);
    free(left_int);
    free(right_int);
}


int* slice_int(int *arr, int start, int end)
{
    int *result;
    int i;
    result = (int *) malloc((end - start) * sizeof(int));
    for (i = start; i < end; i++)
    {
        result[i - start] = arr[i];
    }
    return result;
}
double* slice_double(double *arr, int start, int end)
{
    double *result;
    int i;
    result = (double *) malloc((end - start) * sizeof(double));
    for (i = start; i < end; i++)
    {
        result[i - start] = arr[i];
    }
    return result;
}

double **row_concatenation(double **V, double *eigenvalues, int numPoints){
    int i,j;
    double **result;
    result = Matrix_init(numPoints+1, numPoints);
    for(i=0;i<numPoints+1;i++){
        for(j=0;j<numPoints;j++){
            if(i==0) result[i][j] = eigenvalues[j];
            else result[i][j] = V[i-1][j];
        }
    }
    return result;
}

int is_sym(double **mat, int rows, int cols){
    int i,j;
    if(rows != cols) return 0;
    for(i=0;i<rows;i++){
        for(j=0;j<cols;j++){
            if(mat[i][j]!=mat[j][i]) return 0;
        }
    }
    return 1;
}

void freeMatrix(double** mat, int len){
    int i;
    for (i = 0; i < len; ++i) {
        free(mat[i]);
    }
    free(mat);
}