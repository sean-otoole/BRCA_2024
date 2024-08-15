import numpy as np
from numpy.linalg import norm

def calculate_similarity(A,B):
    cosine = np.dot(A,B)/(norm(A)*norm(B))
    return cosine