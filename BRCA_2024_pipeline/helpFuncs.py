import numpy as np  # Import NumPy for mathematical operations
from numpy.linalg import norm  # Import the norm function from numpy.linalg to compute vector norms

def calculate_similarity(A, B):
    """
    Function to calculate the cosine similarity between two vectors A and B.
    
    Cosine similarity is a measure of similarity between two non-zero vectors
    based on the cosine of the angle between them. It is often used in text mining
    and other fields to measure similarity between vectors.
    
    Args:
        A (numpy array): First input vector (e.g., gene expression profile of a cell).
        B (numpy array): Second input vector (e.g., reference pseudobulk profile).
        
    Returns:
        float: Cosine similarity value between A and B. Ranges from -1 (opposite direction)
               to 1 (same direction), with 0 meaning no similarity.
    """
    # Compute the cosine similarity between A and B using the dot product and vector norms
    cosine = np.dot(A, B) / (norm(A) * norm(B))  # Dot product of A and B divided by the product of their norms
    return cosine  # Return the computed cosine similarity
