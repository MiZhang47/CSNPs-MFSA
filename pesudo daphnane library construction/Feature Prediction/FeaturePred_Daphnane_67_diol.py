"""
Last updated: 2024/12
Structure-based MS feature prediction for normal daphnane with 6,7-diol
"""

import os
import pandas as pd

# Define shift value of A parts
A_shift_set = {
    'A1': {'C': 0, 'H': -4, 'O': 0},  # 1,2-en-3-one
    'A2': {'C': 0, 'H': -2, 'O': 0},  # 1,2-dihydro-3-one
    'A3': {'C': 0, 'H': -2, 'O': 0},  # 1,2-en-3-ol
    'A4': {'C': 0, 'H': 0, 'O': 0},  # 1,2-dihydro-3-ol, reference A ring structure
    'A5': {'C': 0, 'H': -4, 'O': 0},  # 1,10-en-3-one
}

# Define shift value of B parts
B_shift_set = {
    'B2': {'C': 0, 'H': 0, 'O': 0},  # 5-ol-6,7-diol, reference B ring structure
    'B5': {'C': 0, 'H': 0, 'O': -1},  # 5-dehydro-6,7-diol
}

# Define shift value of C parts
C_shift_set = {
    'C1': {'C': 0, 'H': 0, 'O': 0},  # 9,13,14-orthoester, reference C ring structure
    'C2': {'C': 0, 'H': 0, 'O': 1},  # 12-ol and 9,13,14-orthoester
    'C3': {'C': 0, 'H': 0, 'O': 0},  # 9,13,14-triol
    'C4': {'C': 0, 'H': 0, 'O': 1},  # 12-ol and 9,13,14-triol
}

# Define reference feature molecular formula for normal daphnane with 6,7-epoxy
# Based on repetition count for six daphnane standards without C-12 acylation
ref_feature_setB = [  # reference MS features from DPSft-2
    {'C': 20, 'H': 28, 'O': 6},  # charge = +1, Theo.mass 365.19587
    {'C': 20, 'H': 26, 'O': 5},  # charge = +1, Theo.mass 347.18530
    {'C': 20, 'H': 24, 'O': 4},  # charge = +1, Theo.mass 329.17474
    {'C': 20, 'H': 22, 'O': 3},  # charge = +1, Theo.mass 311.16417
    {'C': 19, 'H': 24, 'O': 3},  # charge = +1, Theo.mass 301.17982
    {'C': 20, 'H': 20, 'O': 2},  # charge = +1, Theo.mass 293.15361
    {'C': 19, 'H': 22, 'O': 2},  # charge = +1, Theo.mass 283.16926
    {'C': 19, 'H': 20, 'O': 2},  # charge = +1, Theo.mass 281.15361
    {'C': 19, 'H': 20, 'O': 1},  # charge = +1, Theo.mass 265.15869
    {'C': 17, 'H': 20, 'O': 2},  # charge = +1, Theo.mass 257.15361
]

path_csv = os.path.join('C:/Users/zhang/Desktop/StrucSimulation5/', 'D_67_diol(40).csv')
df = pd.read_csv(path_csv)  # Read the CSV file using the pandas and store the data in a DataFrame object named df.

for index, row in df.iterrows():  # Iterate over each row in the DataFrame.
    new_feature_set = []  # Create an empty list to store the calculated new feature molecular formulas.
    for feature in ref_feature_setB:  # Traverse the reference feature set to calculate the new molecular formula.
        # Create a dictionary to store the calculated new feature molecular formulas.
        modified_feature = {
            'C': feature['C'] + A_shift_set[row['A parts']]['C'] + B_shift_set[row['B parts']]['C']
                + C_shift_set[row['C parts']]['C'],
            'H': feature['H'] + A_shift_set[row['A parts']]['H'] + B_shift_set[row['B parts']]['H']
                + C_shift_set[row['C parts']]['H'],
            'O': feature['O'] + A_shift_set[row['A parts']]['O'] + B_shift_set[row['B parts']]['O']
                + C_shift_set[row['C parts']]['O'],
        }
        new_feature_set.append(modified_feature)  # Add the calculated new feature formula to the new_feature_set list.

    for i, feature in enumerate(new_feature_set, 1):
        # Iterate over each element in the new_feature_set list.
        # Return an enumeration object that yields a pair of values:
        # the index (starting from 1) and the corresponding element (dictionary of feature molecular formulas)
        if feature['O'] == 0:  # Check whether the number of O atoms in the current feature molecular formula is 0.
            df.loc[index, f'MS Feature {i}'] = f"C{feature['C']}H{feature['H']}"
            # Store the new feature formula (excluding the O atom) as a string in the current row (index)
            # of the DataFrame with column name "MS Feature i".
        else:
            df.loc[index, f'MS Feature {i}'] = f"C{feature['C']}H{feature['H']}O{feature['O']}"
            # Store the new feature formula (including the O atom) as a string in the current row (index)
            # of the DataFrame with column name "MS Feature i".

output_path = os.path.join('C:/Users/zhang/Desktop/FeaturePredict5/', 'D_67_diol(40)_feature.csv')
df.to_csv(output_path, index=False)