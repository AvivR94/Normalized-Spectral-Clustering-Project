import pandas as pd
import numpy as np
import sys
import spkmeansmodule as spkmm

np.random.seed(0)
def initialize_centroids(vectors, k):
    centroids = []
    randomindex = np.random.randint(n)
    centroids.append(vectors[randomindex])
    prob=[0.0 for i in range(n)]
    distances=[0.0 for i in range(n)]
    index_list=[x for x in range(n)]
    indexes_used=[randomindex]
    randomindex = np.random.randint(n)
    i=0
    while (i<k):
        for j in range(n):
            min=float('inf')
            for q in range(len(centroids)):
                distance=0
                for m in range(d):
                    distance+=(float(vectors[j][m])-float(centroids[q][m]))**2
                if distance<min:
                    min = distance
            distances[j]=min
        Sum = sum(distances)
        for j in range(n):
            prob[j]=distances[j]/Sum
        new_index = (np.random.choice(index_list, 1, p=prob))[0]
        indexes_used.append(new_index)
        centroids.append(vectors[new_index])
        i+=1
    return centroids, indexes_used

#main

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Invalid Input!")
        sys.exit()
    try:
        k = int(sys.argv[1]) 
        max_iter = 300
        goal = sys.argv[2]
        possible_goals = {"wam","ddg","lnorm","jacobi","spk"}
        if not goal in possible_goals:
            print("Invalid Input!")
            sys.exit()
        input = np.loadtxt(sys.argv[3], delimiter=',')
    except ValueError:
        print("Invalid Input!")
        sys.exit()
        
    n = input.shape[0]
    d = input.shape[1]
    if (k>=n or k<0):
        print ("Invalid Input!")
        sys.exit()
    if(k==1 and goal=="spk"):
        print("Invalid Input!")
        sys.exit()
    matT = spkmm.pythonGetMatR(input.flatten().tolist(), k, n, d, goal)
    if(goal == 'spk'):
            np.set_printoptions(formatter={'float_kind': '{:.4f}'.format})
            if(k == 0):
                k = int(matT.pop())
            else:
                matT.pop()
            matT = np.array(matT).reshape(-1, k)
            dim = np.shape(matT)
            centroids, indexes_used = initialize_centroids(matT, k)
            importlib.reload(spkmm)
            lst = spkmm.kmeanspp(centroids.tolist(), matT.tolist(), k, 300, dim[0], dim[1])
            result = np.ndarray(shape=(k, dim[1]))
            for i in range(k):
                for j in range(dim[1]):
                    result[i][j] = lst[(i * dim[1]) + j]
            res = ""
            for i in indexes_used:
                res += str(i) + ","
            print(res[:-1])
            res = ""
            for i in result:
                for j in i:
                    res += str(format(np.round(j, 4), '.4f')) + ','
                print(res[:-1])
                res = ""


