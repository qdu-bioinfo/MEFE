import numpy as np
import pandas as pd


def compute_weights(abundance_matrix, adjacency_list, alpha=0.8, missing_value_strategy='mean'):
    '''
    Computes local attention weights for each OTU based on its neighbors' similarities.
    Updates the abundance matrix with new feature values using the attention mechanism.

    Parameters:
    abundance_matrix (pd.DataFrame): The input matrix containing OTU abundances.
    adjacency_list (dict): A dictionary where keys are OTUs and values are dictionaries of neighboring OTUs and their similarity values.
    alpha (float): A weighting factor for the original abundance (default 0.8).
    missing_value_strategy (str): The strategy to handle missing values ('mean', 'zero', 'remove').

    Returns:
    pd.DataFrame: A new matrix with updated feature values.
    '''
    # Handle missing values in the abundance matrix
    if missing_value_strategy == 'mean':
        abundance_matrix.fillna(abundance_matrix.mean(), inplace=True)
    elif missing_value_strategy == 'zero':
        abundance_matrix.fillna(0, inplace=True)
    elif missing_value_strategy == 'remove':
        abundance_matrix.dropna(inplace=True)

    new_features = pd.DataFrame(index=abundance_matrix.index,
                                columns=abundance_matrix.columns)

    for otu, neighbors in adjacency_list.items():
        if otu not in abundance_matrix.columns:
            continue

        similarities = list(neighbors.values())

        # Compute attention weights using softmax
        attention_weights = softmax(np.array(similarities)) if similarities else []

        new_feature_values = np.zeros(len(abundance_matrix))

        for i, neighbor_otu in enumerate(neighbors):
            if neighbor_otu in abundance_matrix.columns:
                new_feature_values += attention_weights[i] * abundance_matrix[neighbor_otu].values

        new_features[otu] = alpha * abundance_matrix[otu] + (1 - alpha) * new_feature_values

    return new_features


def softmax(x):
    '''
    Computes the softmax function for a given input array.

    Parameters:
    x (numpy.ndarray): The input array of values.

    Returns:
    numpy.ndarray: The softmax transformed values.
    '''
    e_x = np.exp(x - np.max(x))  # Subtract max for numerical stability
    return e_x / e_x.sum()
