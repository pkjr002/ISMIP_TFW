import os
import pandas as pd

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
    # Store the folder paths and their split components
    folder_paths = []
    folder_levels = []

    for root, dirs, files in os.walk(path, topdown=True):
        # Sort directories in-place to ensure alphabetical order in traversal
        dirs.sort()
        
        # Normalize root with consistent separators and remove trailing separator
        normalized_root = os.path.normpath(root)
        # Split the path into parts
        parts = normalized_root.split(os.sep)
        folder_paths.append(normalized_root)
        folder_levels.append(parts)

    # Creating a sorted list from the folder_levels to ensure alphabetical order
    folder_levels_sorted = sorted(folder_levels, key=lambda x: tuple(x))
    folder_paths_sorted = ['/'.join(level) for level in folder_levels_sorted]

    # Convert the sorted list of folder levels into a MultiIndex
    multi_index = pd.MultiIndex.from_tuples(folder_levels_sorted, names=[f'Level_{i+1}' for i in range(max(len(x) for x in folder_levels))])

    # Create a DataFrame with the MultiIndex
    df_folders = pd.DataFrame(index=multi_index)
    df_folders['Folder Path'] = folder_paths_sorted

    return df_folders