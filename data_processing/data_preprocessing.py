import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from sklearn.preprocessing import LabelEncoder


pandas2ri.activate()


def read_id_file(file_path):

    '''
    Reads the ID file and returns two dictionaries:
    1. seq_to_id: Maps seq_id to id_value.
    2. id_to_seq: Maps id_value to seq_id.
    '''

    seq_to_id = {}  # Mapping from seq_id to id_value
    id_to_seq = {}  # Mapping from id_value to seq_id

    with open(file_path, 'r') as file:
        for line in file:
            seq_id, id_value = line.strip().split('\t')
            seq_to_id[seq_id] = id_value
            id_to_seq[id_value] = seq_id

    return seq_to_id, id_to_seq




def parse_abundance_file(file_path, id_dict):

    abundance_matrix = pd.read_csv(file_path,
                                   sep='\t',
                                   index_col=0)
    abundance_matrix.columns = [id_dict.get(col, col)
                                for col in abundance_matrix.columns]

    return abundance_matrix


def parse_adjacency_file(file_path, id_dict):

    adjacency_list = {}

    with open(file_path, 'r') as file:

        for line in file:

            parts = line.strip().split('\t')
            otu = id_dict.get(parts[0], parts[0])
            neighbors = {id_dict.get(parts[i], parts[i]):
                             float(parts[i + 1]) * 0.01

                         for i in range(1, len(parts), 2)}

            adjacency_list[otu] = neighbors

    return adjacency_list


def read_data(meta_path, csv_path, new_matrix_path=None, features_path=None):

    meta_df = pd.read_csv(meta_path,
                          header=None,
                          sep='\t')
    csv_df = pd.read_csv(csv_path,
                         index_col=False,
                         sep='\t')

    new_matrix_df = pd.read_csv(new_matrix_path,
                                index_col=False,
                                sep='\t') \
        if new_matrix_path else None
    significant_features = pd.read_csv(features_path,
                                       header=None,
                                       sep='\t').iloc[1:, 0].tolist() \
        if features_path else None

    return meta_df, csv_df, new_matrix_df, significant_features



def save_selected_features_to_file(X_selected, output_path):

    X_selected.to_csv(output_path,
                      sep='\t',
                      index=False)

    print(f"Selected features saved to {output_path}")



def preprocess_data(meta_df, csv_df):

    '''
    Preprocesses the meta and abundance data:
    - Maps sample IDs to their groups.
    - Fills missing values with 0.
    - Returns processed X_train and Y.
    '''

    label_map = dict(zip(meta_df.iloc[:, 0], meta_df.iloc[:, 1]))

    csv_df.insert(1, 'Group', csv_df.iloc[:, 0].map(label_map))

    X_train = csv_df.iloc[:, 2:]
    Y = csv_df.iloc[:, 1]

    label_encoder = LabelEncoder()

    Y = label_encoder.fit_transform(Y)
    X_train = X_train.fillna(0)

    return X_train, Y, label_encoder



def preprocess_data_with_filter(X_train, sample_threshold=0.1, abundance_cutoff=0.00001):

    sample_count_threshold = int(sample_threshold * X_train.shape[0])
    filtered_data = X_train.loc[:, (X_train > abundance_cutoff).sum(axis=0) >= sample_count_threshold]

    print(f"Original number of features: {X_train.shape[1]}")
    print(f"Final number of features: {filtered_data.shape[1]}")

    return filtered_data



def wilcox_test_r(group1, group2, paired=False):

    '''
    Performs a Wilcoxon test (either paired or unpaired) on two groups of data.
    Returns the p-value of the test.
    '''

    robjects.globalenv['group1'] = pandas2ri.py2rpy(group1)
    robjects.globalenv['group2'] = pandas2ri.py2rpy(group2)

    test_command = 'wilcox.test(group1, group2, paired=TRUE, exact=FALSE, correct=FALSE)$p.value' if paired else \
                   'wilcox.test(group1, group2, exact=FALSE, correct=FALSE)$p.value'
    p_value = robjects.r(test_command)

    return p_value[0]




def statistical_tests(X_train, Y, paired=False):

    X_train_filtered = preprocess_data_with_filter(X_train)

    p_values_wilcoxon = [wilcox_test_r(X_train_filtered[Y == 0][col],
                                       X_train_filtered[Y == 1][col],
                                       paired)
                         for col in X_train_filtered.columns]

    significant_features_wilcoxon = [col for col, p in zip(X_train_filtered.columns, p_values_wilcoxon) if p < 0.01]
    print(f"Number of significant features (Wilcoxon test): {len(significant_features_wilcoxon)}")

    return significant_features_wilcoxon

