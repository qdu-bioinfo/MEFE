from data_processing.data_preprocessing import (read_id_file,
                                                parse_adjacency_file,
                                                parse_abundance_file,
                                                save_selected_features_to_file,
                                                read_data,
                                                wilcox_test_r)

from data_processing.data_preprocessing import (preprocess_data,
                                                statistical_tests)
from model.model import compute_weights

import pandas as pd


def main():
    '''
    Main function that orchestrates the reading of data, processing, selection of significant features,
    and saving the results to files.
    '''

    # Base directory for the project
    base_dir = './'  # Current directory (the root of the project)

    # Define relative file paths
    database_dir = 'database'
    data_dir = 'data'
    output_dir = 'output'

    id_file_path = f'{database_dir}/id.txt'
    adjacency_file_path = f'{database_dir}/gg2_index.txt'

    abundance_file_path = f'{data_dir}/Abd.tab'
    meta_file_path = f'{data_dir}/meta.txt'

    elastic_file_path = f'{output_dir}/elastic_matrix.tab'
    elastic_matrix_file_path = elastic_file_path
    biomarker_ate_path = f'{output_dir}/biomarker_matrix.txt'
    biomarker_p_values_path = f'{output_dir}/biomarkers.txt'

    # Read ID file and get seq_to_id and id_to_seq mappings
    seq_to_id, id_to_seq = read_id_file(id_file_path)

    # Parse adjacency list and abundance matrix
    adjacency_list = parse_adjacency_file(adjacency_file_path, seq_to_id)
    abundance_matrix = parse_abundance_file(abundance_file_path, seq_to_id)

    # Compute the new abundance matrix
    new_abundance_matrix = compute_weights(abundance_matrix, adjacency_list)

    # Save the updated abundance matrix
    new_abundance_matrix.to_csv(elastic_file_path, sep='\t')
    print(f"Updated abundance matrix saved to {elastic_file_path}")

    # Read data and preprocess
    meta_df, csv_df, matrix_df, file_features = read_data(meta_file_path,
                                                              elastic_matrix_file_path,
                                                              elastic_matrix_file_path)

    # Preprocess and select significant features
    X_train, Y, label_encoder = preprocess_data(meta_df, csv_df)
    significant_features = file_features if file_features else statistical_tests(X_train, Y, paired=False)

    print(f"Number of significant features from Wilcoxon test: {len(significant_features)}")

    # Select features from the new matrix based on significant features
    X_selected = matrix_df[significant_features] if matrix_df is not None else None

    # Map features back to Seq_id
    significant_features_seqid = [seq_to_id.get(feature, feature) for feature in significant_features]

    # Insert sample IDs
    meta_sample_ids = meta_df.iloc[1:, 0].values
    X_selected.insert(0, 'Sample_ID', meta_sample_ids[:len(X_selected)])

    # Map columns (Seq_id) back to original seq_id for biomarker_ate.txt
    original_seq_id_dict = id_to_seq
    X_selected.columns = [original_seq_id_dict.get(col, col) for col in X_selected.columns]

    # Save selected features to biomarker_ate.txt with restored seq_id
    save_selected_features_to_file(X_selected, biomarker_ate_path)

    # Compute Wilcoxon p-values
    p_values_wilcoxon = []
    for feature in significant_features:
        group1 = X_train[Y == 0][feature]
        group2 = X_train[Y == 1][feature]
        p_value = wilcox_test_r(group1, group2, paired=False)
        p_values_wilcoxon.append(p_value)

    # Create a DataFrame to store Seq_id and p-value
    p_values_df = pd.DataFrame({
        'Seq_id': significant_features_seqid,
        'p_value': p_values_wilcoxon
    })

    # Map Seq_id back to original seq_id using the id_dict
    p_values_df['Seq_id'] = p_values_df['Seq_id'].map(id_to_seq)

    # Save the DataFrame to biomarker_p_values.txt
    p_values_df.to_csv(biomarker_p_values_path, sep='\t', index=False)
    print(f"Biomarkers to {biomarker_p_values_path}")




if __name__ == "__main__":
    main()
