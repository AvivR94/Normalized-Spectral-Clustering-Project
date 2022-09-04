import sys
import numpy as np
import spkmeansmodule as spkmm

def initialize_centroids(vectors, k, n): #points is matrix of points
    np.random.seed(0)
    centroids = []
    random_index = np.random.randint(n)
    centroids.append(vectors[random_index])
    indexes_used = [random_index]
    index_list=[x for x in range(n)]
    i = 1
    while (i<k):
        distances = retrieve_distances(vectors, centroids)
        prob = retrieve_prob(distances)
        new_index = (np.random.choice(index_list, 1, p=prob))[0]
        centroids.append(vectors[new_index])
        indexes_used.append(new_index)
        i += 1
    return centroids, indexes_used

def input_to_float(input,n):
    for i in range(n):
        vector = input[i].split(",")
        for j in range(len(vector)):
            vector[j] = float(vector[j])
        input[i] = vector

def retrieve_distances(vectors, centroids):
    n = len(vectors)
    distances = [0 for i in range(n)]
    for j in range(n):
        min=float('inf')
        for q in range(len(centroids)):
            dist = distance_calc(vectors[j], centroids[q])
            if dist < min:
                min = dist
        distances[j] = min
    return distances
    
def distance_calc(x, y):
    dist=0
    for i in range(len(x)):
        dist+=(float(x[i])-float(y[i]))**2
    return dist

def retrieve_prob(distances):
    n = len(distances)
    Sum = sum(distances)
    prob = [0.0 for i in range(n)]
    for i in range(n):
        prob[i] = distances[i]/Sum
    return prob
    
def retrieve_flatten_mat(mat,r,c):
    lst = []
    for i in range (r):
        for j in range (c):
            lst.append(mat[i][j])
    return lst

#helping functions to print matrix in python
def spk_print(final_centroids,indexes):    
    for i in range(len(indexes)):
        indexes[i] = str(indexes[i])

    indexes = ",".join(indexes)
    print(indexes)
    print_matrix(final_centroids, len(final_centroids), len(final_centroids[0]))

def print_matrix(mat, r, c):
    for i in range(r):
        for j in range(c):
            mat[i][j] = '%.4f'%mat[i][j]
        result = []
    for i in range(r):
        result.append(",".join(mat[i]))
    for i in range(r):
        if (i != r-1):
            print(",".join(result[i]))
        else:
            print(",".join(result[i])+"\n")

def print_eigens(eigens):
    for i in range(len(eigens)):
        eigens[i] = str(eigens[i])
    print(",".join(eigens))

# main 
if __name__ == '__main__':
    
    if len(sys.argv) != 4:
        print("Invalid Input!")
        sys.exit()
    #read arguments
    try:
        k = int(sys.argv[1]) 
        max_iter = 300
        goal = sys.argv[2]
        possible_goals = {"wam","ddg","lnorm","jacobi","spk"}
        if goal not in possible_goals:
            print("Invalid Input!")
            sys.exit()
        filename = sys.argv[3]
    except ValueError:
        print("Invalid Input!")
        sys.exit()
    #read the input file
    try:
        f = open(filename, "r")
        reading_file = f.read()
        input = reading_file.splitlines()
        f.close()
    except:
        print("An Error Has Occurred")
        sys.exit()
    
    n = len(input)
    #check for validity of k
    if (k==1 and goal=="spk") or k>=n or k<0: 
        print("Invalid Input!")
        sys.exit()
        
    input_to_float(input,n)
    d = len(input[0])
    flatten_input = retrieve_flatten_mat(input,n,d)
    final_mat = spkmm.get_matrix(k, n, d, flatten_input, goal)
    r = len(final_mat)
    c = len(final_mat[0])
    #printing goal matrix we got from C's goal router here
    #if goal is not spk we will print the matrix we got from goal router here
    if goal in {"wam", "ddg", "lnorm"}:
        print_matrix(final_mat,r,c)

    # if goal == 'jacobi', print is made on C code
    # An instrucor in the moodle forum approved this
    
    #if goal is spk we will send T (that we got from goal router) to Kmeans algorithm in C
    if goal == "spk":
        k = len(final_mat[0])
        centroids, indexes = initialize_centroids(final_mat, k, n) # these are the first centroids and their indexes
        initial_centroids_fit_input = retrieve_flatten_mat(centroids,k,k)
        T_fit_input = retrieve_flatten_mat(final_mat,n,k)
        # getting final centroids out of  T as input points :
        final_centroids= spkmm.fit(k, n, k,initial_centroids_fit_input, T_fit_input)
        # printing the initial indexes and the final centroids (as did in HW2):
        spk_print(final_centroids,indexes)