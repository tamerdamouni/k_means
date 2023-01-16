import pandas as pd
import numpy as np
import kmeansspk as kmp
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("k", type=int)
parser.add_argument("goal", type=str)
parser.add_argument("input", type=str)
args = parser.parse_args()

max_iter = 300
epsilon=0

def calc_norm(centroids, point):  # norm calculator functions below
    closest_norm = 0
    closest_norm_value = diff_norm(centroids[0], point)
    for i in range(len(centroids)):
        x = diff_norm(centroids[i], point)
        if x < closest_norm_value:
            closest_norm = i
            closest_norm_value = x
    return closest_norm


def diff_norm(vector1, vector2):
    arr = [0 for i in range(len(vector2))]
    for i in range(len(vector2)):
        arr[i] = vector1[i] - vector2[i]
    summation = 0
    for i in range(len(vector2)):
        summation += pow(arr[i], 2)
    return summation


def calc_sum(list1, list2):  # A function that returns the sum of two vectors
    summation = []
    for i in range(len(list1)):
        summation.append(list1[i] + list2[i])
    return summation


def divide(vector, num):  # dividing all elements in a vector
    for i in range(len(vector)):
        vector[i] = vector[i] / num
    return vector

def init_centroids(k, all_points):
    N = len(all_points)
    i=0
    indexes = [0 for i in range(k)]
    np.random.seed(0)
    indexes[0] = np.random.choice(N)
    mu0 = all_points[indexes[0]]
    mu = np.zeros((k, len(all_points[0])))
    mu[0] = mu0
    while i < k:
        i+=1
        if(i==k) :break
        D = [0 for s in range(N)]
        for l in range(0,N):
            minimum = diff_norm(all_points[l], mu[0])
            for j in range(0,i):
                calculation = diff_norm(all_points[l], mu[j])
                if calculation < minimum:
                    minimum = calculation
            D[l] = minimum
        denominator = 0
        for m in range(0,N):
            denominator = denominator + D[m]
        p = [0 for i in range(N)]
        for l in range(0,N):
            p[l] = D[l] / denominator
        indexes[i] = np.random.choice(N, p=p)
        mu[i] = all_points[indexes[i]]
    print(','.join(map(str, indexes)))
    return mu, indexes


def printMatrix(matrix):
    for i in range(len(matrix)):
        print(','.join([format(matrix[i][j], ".4f") for j in range(len(matrix[i]))]))


def calc_centroids(k, max_iter, d,  all_points, initial_centroids):
    initial_centroids = initial_centroids.tolist()
    all_points = all_points.tolist()
    centroids = kmp.finalspk(k, max_iter, d, all_points, initial_centroids)
    printMatrix(centroids)


def execute(k, goal, max_iter, all_points):
    if goal == 'wam':
        mat = kmp.fit(all_points, 0, len(all_points), len(all_points[0]))
    elif goal == 'ddg':
        mat = kmp.fit(all_points, 1, len(all_points), len(all_points[0]))
    elif goal == 'lnorm':
        mat = kmp.fit(all_points, 2, len(all_points), len(all_points[0]))
    elif goal == 'jacobi':
        mat = kmp.fit(all_points, 3, len(all_points), len(all_points[0]))
        for i in range(len(mat[0])):
            if(mat[0][i]<0 and mat[0][i]> -0.00005): mat[0][i]=0
    elif goal == 'spk':
        Tmat = kmp.fitgetTmat(all_points, k, len(all_points), len(all_points[0]))
        Tmat = np.array(Tmat)
        k = len(Tmat[0])
        centroids = init_centroids(k, Tmat)[0]
        calc_centroids(k, max_iter, k, Tmat, centroids)
    else:
        print("invalid input")
        return
    if goal != 'spk': printMatrix(mat)


def valid(input):
    if not (str(input[0]).isdigit()):
        return 0
    if(input[0]<0): return 0
    return 1


input_args = [args.k, args.goal, args.input]
file = open(args.input,"r")
lines_list = []
curr_line = 0
# file reading:
for line in file:
    str_line = line.split(",")
    int_line = []
    for i in range(len(str_line)):
        int_line.append(float(str_line[i]))
    lines_list.append(int_line)  # lines list[i] will contain the i'th vector in the list
    if line != "\n":
        curr_line += 1
all_points = lines_list

if(len(all_points) < args.k) :
    print("Invalid Input!")
if(not valid(input_args)):
    print("Invalid Input!")
execute(args.k, args.goal, 300, all_points)
