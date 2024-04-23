import os
import pandas as pd
from collections import defaultdict

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# FUNCTION BLOCK.
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
def list_files_with_names(path, names):    
    if not isinstance(names, list):
        raise ValueError("The 'names' argument should be a list of strings.")
    #
    file_list = os.listdir(path)
    file_list = sorted(file_list)
    #
    matching_files = []
    for filename in file_list:
        if all(name in filename for name in names):
            matching_files.append(filename)
    return matching_files


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def list_folders_as_multilevel_df_sorted(path):
    folder_paths = []
    folder_levels = []

    for root, dirs, files in os.walk(path, topdown=True):
        dirs.sort()
        
        normalized_root = os.path.normpath(root)
        # Split the path into parts
        parts = normalized_root.split(os.sep)
        folder_paths.append(normalized_root)
        folder_levels.append(parts)

    folder_levels_sorted = sorted(folder_levels, key=lambda x: tuple(x))
    folder_paths_sorted = ['/'.join(level) for level in folder_levels_sorted]

    multi_index = pd.MultiIndex.from_tuples(folder_levels_sorted, names=[f'Level_{i+1}' for i in range(max(len(x) for x in folder_levels))])

    df_folders = pd.DataFrame(index=multi_index)
    df_folders['Folder Path'] = folder_paths_sorted

    return df_folders


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# paths grouped by the first 5 levels
def list_folders_as_multilevel_df_grouped(path):
    grouped_paths = defaultdict(list)
    for root, dirs, files in os.walk(path, topdown=True):
        dirs.sort()
        normalized_root = os.path.normpath(root)
        parts = normalized_root.split(os.sep)
        
        # Group by the first 5 levels and aggregate the rest under Level 6
        key = tuple(parts[:5])
        grouped_paths[key].append('/'.join(parts[5:])) if len(parts) > 5 else grouped_paths[key].append('/')

    data = []
    for key, values in grouped_paths.items():
        row = list(key) + [values]  
        data.append(row)
    
    # Find the maximum depth for naming columns, considering the aggregated lists as the last level
    max_depth = max(len(row) for row in data)
    column_names = [f'Level_{i+1}' for i in range(max_depth-1)] + ['Level_6_Subfolders']
    
    df_folders = pd.DataFrame(data, columns=column_names)
    
    return df_folders



from collections import defaultdict

def list_folders_as_multilevel_df_grouped1(path, aggregation_level=5):
    """
    Lists folders in a multi-level DataFrame, grouped by a specified level with remaining paths aggregated.
    
    Parameters:
    - path (str): The root path to start listing folders from.
    - aggregation_level (int): The depth level at which to start aggregating subfolders into a list. Defaults to 5.
    
    Returns:
    - pd.DataFrame: A DataFrame with each row representing a path, grouped up to the specified level,
                    and additional subfolders aggregated into a list in the last column.
    """
    grouped_paths = defaultdict(list)
    for root, dirs, files in os.walk(path, topdown=True):
        dirs.sort()
        normalized_root = os.path.normpath(root)
        parts = normalized_root.split(os.sep)
        
        # Adjust to group by the specified levels and aggregate the rest
        key = tuple(parts[:aggregation_level])
        grouped_paths[key].append('/'.join(parts[aggregation_level:])) if len(parts) > aggregation_level else grouped_paths[key].append('/')
    
    data = []
    for key, values in grouped_paths.items():
        row = list(key) + [values]  
        data.append(row)
    
    # Adjust the column names based on the specified aggregation level
    column_names = [f'Level_{i+1}' for i in range(aggregation_level)] + [f'Level_{aggregation_level}_Subfolders']
    
    df_folders = pd.DataFrame(data, columns=column_names)
    
    return df_folders
