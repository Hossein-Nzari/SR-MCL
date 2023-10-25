# Import necessary libraries
from scipy.sparse import isspmatrix, dok_matrix, csc_matrix
import sklearn.preprocessing
import scipy
import math
import numpy as np
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt

# Load data from an Excel file
data = pd.read_excel('pro200600448_1_s.xls')
df = pd.DataFrame(data, columns= ['Protein A','Protein B','Score'])
G = nx.Graph()

# Create a graph and add nodes and edges based on the data
for i in range(len(df['Protein A'])):
    G.add_node(df['Protein A'][i])
    G.add_node(df['Protein B'][i])
    G.add_edge(df['Protein A'][i], df['Protein B'][i], weight=df['Score'][i])

# Convert the graph to an adjacency matrix
m = nx.adjacency_matrix(G)
nodes = list(G.nodes)

# Set the diagonal elements of the adjacency matrix to their maximum values
a = nx.to_numpy_array(G)
for i in range(len(nodes)):
    a[i][i] = max(a[i])

# Normalize the adjacency matrix
mg = a.copy()
sum = np.sum(a, axis=0)
for i in range(len(nodes)):
    for j in range(len(nodes)):
        mg[i][j] = a[i][j] / sum[j]

# Create a new graph from the normalized adjacency matrix
G_mg = nx.from_numpy_array(mg)
mg = nx.adjacency_matrix(G_mg)

# Set some parameters
n = len(G)
b = 0.5
r = 2
B = 1.25
om = 0.8
t = 3
c = []
count = []
for i in range(n):
    count.append(0)

# Function to prune matrix values below a threshold
def prune(matrix, threshold):

    if isspmatrix(matrix):
        pruned = dok_matrix(matrix.shape)
        pruned[matrix >= threshold] = matrix[matrix >= threshold]
        pruned = pruned.tocsc()
    else:
        pruned = matrix.copy()
        pruned[pruned < threshold] = 0

    num_cols = matrix.shape[1]
    row_indices = matrix.argmax(axis=0).reshape((num_cols,))
    col_indices = np.arange(num_cols)
    pruned[row_indices, col_indices] = matrix[row_indices, col_indices]

    return pruned

# Function to normalize a matrix
def normalize(matrix):
    return sklearn.preprocessing.normalize(matrix, norm='l2', axis=0)

# Function to calculate the regularization matrix
def regularizationMatrix(m, mg, b):
    arr = []
    global n
    for i in range(n):
        arr.append(1)
    cv1 = np.array([arr]).T
    mass = m.dot(cv1)

    mt = m.T
    p = mt.dot(mass)
    P = np.zeros((n, n), int)
    np.fill_diagonal(P, p)

    mr = mg.dot(scipy.linalg.fractional_matrix_power(P, -b))
    mr = normalize(mr)
    return mr

# Function to apply inflation to a matrix
def inflate(matrix, power, count, B):

    if isspmatrix(matrix):
        return normalize(matrix.power(power))


    for i in range(len(count)):
        if count[i] > 0:
            for j in range(len(count)):
                matrix[i][j] += r * B ** count[i]
                if matrix[i][j] <= 0:
                    matrix[i][j] = 0
    return normalize(np.power(matrix, power))


def sparse_allclose(a, b, rtol=1e-5, atol=1e-8):
    c = np.abs(a - b) - rtol * np.abs(b)
    return c.max() <= atol

# Function to check if two matrices have converged
def converged(matrix1, matrix2):
    if isspmatrix(matrix1) or isspmatrix(matrix2):
        return sparse_allclose(matrix1, matrix2)

    return np.allclose(matrix1, matrix2)

# Function to get clusters from a matrix
def get_clusters(matrix):

    if not isspmatrix(matrix):
        matrix = csc_matrix(matrix)

    global attractors
    attractors = matrix.diagonal().nonzero()[0]
    for i in attractors:
        count[i] += 1
    clusters = set()
    for attractor in attractors:
        cluster = tuple(matrix.getrow(attractor).nonzero()[1].tolist())
        clusters.add(cluster)
    return sorted(list(clusters))

# Iteratively apply SR-MCL
for i in range(t):

    lastMatrix = m.copy()
    mr = regularizationMatrix(m, mg, b)
    m = m.dot(mr)
    m = inflate(m, r, count, B)
    m = prune(m, om)

    ci = get_clusters(m)
    c = list(set(c) | set(ci))
    if converged(m, lastMatrix):
        break

new_c = []
for i in c:
    new_c.append(list(i))
c = new_c

for i in range(len(c)):
    for j in range(len(c[i])):
        c[i][j] = nodes[c[i][j]]

# Post-processing to refine clusters
def postProcessing(c, om, p):
    num = 0
    while num != len(c):
        if len(c[num]) <= 2:
            del c[num]
            num -= 1
        num += 1

    num = 0
    sortqf = []

    while num != len(c):
        G_clu = nx.Graph()
        for j in c[num]:
            G_clu.add_node(j)
        for j in c[num]:
            for k in c[num]:
                if j != k and (G.has_edge(j, k)):
                    G_clu.add_edge(j, k, weight=G[j][k]["weight"])
        if nx.density(G_clu) * math.sqrt(len(list(G_clu.nodes))) <= om:
            del c[num]
            num -= 1

        else:
            sortqf.append(nx.density(G_clu) * math.sqrt(len(list(G_clu.nodes))))
            #nx.draw(G_clu, with_labels = True)
            #plt.show()

        num += 1
    print('mean qf: ', np.sum(sortqf) / len(sortqf))
    sortedqf = [x for _, x in sorted(zip(sortqf, c))]
    sortedqf.reverse()

    num = 0
    while num != len(c):
        for j in range(num + 1, len(c)):
            intersection = [value for value in c[num] if value in c[j]]
            if (len(intersection)) ** 2 / (len(c[j]) * len(c[j])) >= p:
                del c[num]
                j -= 1
        num += 1

    return c

# Apply post-processing to the clusters
c = postProcessing(c, om, 0.6)

# Display cluster information
print('number of clusters: ', len(c))
print(c)

#nx.draw(G)
#plt.show()

'''
tunable parameters = r, p, ω
Regarding all three parameters, increasing the parameters leads to a decrease in the number of clusters and an increase in the average value of qf.
'''