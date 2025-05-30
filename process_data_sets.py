import pandas as pd
import numpy as np
import imageio
import os

def pre_process_data_set(file_name, nucleus, cytoplasm, resolution_z, resolution_xy, cytokine_analysis):
    data_set = pd.read_csv(file_name)
    if nucleus:
        data_set = data_set[data_set.nuclei_mask > 0.0]
    elif cytoplasm:
        data_set = data_set[data_set.nuclei_mask == 0]

    if cytokine_analysis == 'TNF':
        data_set = data_set[data_set['cytokine'] == 'TNF']
    if cytokine_analysis == 'IFN':
        data_set = data_set[data_set['cytokine'] == 'IFN']

    data_set['Z_det'] = data_set['Z_det'] * resolution_z
    data_set['Y_det'] = data_set['Y_det'] * resolution_xy
    data_set['X_det'] = data_set['X_det'] * resolution_xy
    return data_set

def make_outline_dataset(file_path):
    #this code rewrite the XYZ positions from the outline_files of the NUCLEUS to an usable dataframe

    with open(file_path, 'r') as f:
        lines = f.readlines()

    nucleus_data = []
    current_cell_name = None
    parsing_nucleus = False
    nucleus_section = {"X_POS": [], "Y_POS": [], "Cell_Name": []}

    for line in lines:
        line = line.strip()
        if line.startswith("CELL_START"):
            current_cell_name = line.split("\t")[1]
        elif line.startswith("Nucleus_START"):
            parsing_nucleus = True
            nucleus_section = {"Cell_Name": [], "X_POS": [], "Y_POS": []}
        elif line.startswith("Nucleus_END"):
            parsing_nucleus = False
            if nucleus_section["X_POS"] and nucleus_section["Y_POS"]:
                nucleus_section["Cell_Name"] = [current_cell_name] * len(nucleus_section["X_POS"])
                nucleus_data.append(pd.DataFrame(nucleus_section))
        elif parsing_nucleus:
            if line.startswith("X_POS"):
                nucleus_section["X_POS"] = list(map(int, line.split("\t")[1:]))
            elif line.startswith("Y_POS"):
                nucleus_section["Y_POS"].extend(map(int, line.split("\t")[1:]))
            elif line.startswith("Z_POS"):
                pass

    if nucleus_data:
        outline_df = pd.concat(nucleus_data, ignore_index=True)
    else:
        outline_df = pd.DataFrame(columns=["X_POS", "Y_POS", "Cell_Name"])

    return outline_df

def make_outline_dataset_cyto(file_path):
    #this code rewrite the XYZ positions from the outline_files of the CYTOPLASM to an usable dataframe
    with open(file_path, 'r') as f:
        lines = f.readlines()

    cytoplasm_data = []
    current_cell_name = None
    parsing_cytoplasm = False
    cytoplasm_section = {"X_POS": [], "Y_POS": [], "Cell_Name": []}

    for line in lines:
        line = line.strip()
        if line.startswith("CELL_START"):
            current_cell_name = line.split("\t")[1]
            parsing_cytoplasm = True #new
            cytoplasm_section = {"Cell_Name": [], "X_POS": [], "Y_POS": []} #new
        elif line.startswith("CELL_END"):
            parsing_cytoplasm = False
            if cytoplasm_section["X_POS"] and cytoplasm_section["Y_POS"]:
                cytoplasm_section["Cell_Name"] = [current_cell_name] * len(cytoplasm_section["X_POS"])
                cytoplasm_data.append(pd.DataFrame(cytoplasm_section))
        elif parsing_cytoplasm:
            if line.startswith("X_POS"):
                cytoplasm_section["X_POS"] = list(map(int, line.split("\t")[1:]))
            elif line.startswith("Y_POS"):
                cytoplasm_section["Y_POS"].extend(map(int, line.split("\t")[1:]))
            elif line.startswith("Z_POS"):
                pass

    if cytoplasm_data:
        outline_df = pd.concat(cytoplasm_data, ignore_index=True)
    else:
        outline_df = pd.DataFrame(columns=["X_POS", "Y_POS", "Cell_Name"])

    return outline_df

def find_cell_name(file_path_tiff, final_df):
    image = imageio.v2.imread(file_path_tiff)

    #select where there are cells and find intensities
    z_indices, y_indices, x_indices = np.where(image > 0)
    intensities = image[z_indices, y_indices, x_indices]

    #make a datafram
    positions_df = pd.DataFrame({'Y_POS': y_indices, 'X_POS': x_indices, 'Intensity': intensities})

    #merge the dataframe from the outline with the intensity coordinates
    merged_df = positions_df.merge(final_df, on=['Y_POS', 'X_POS'], how='inner')

    #group the cells and intensities, match_count is the number of overlaps between the two in pixels
    overlap_counts = merged_df.groupby(['Intensity', 'Cell_Name']).size().reset_index(name='Match_count')

    #find the maximal overlap count
    best_matches = overlap_counts.loc[overlap_counts.groupby('Intensity')['Match_count'].idxmax()]
    cell_name_df = best_matches[['Cell_Name', 'Intensity', 'Match_count']]
    cell_name_df = cell_name_df.groupby("Cell_Name")["Intensity"].apply(list).reset_index()

    return cell_name_df, image

def data_prep(file_path_cellpose, file_path_outline_files):
    # this code is made to find the overlap between the nucleus cellpose file and the outline file path (only the names thus!!!)
    # if files are differently saved than explained in the GITHUB, change this code (probably change the for loop structure)

    #cellpose
    data = []
    for item in os.listdir(file_path_cellpose):
        parts = item.split('_')
        if len(parts) > 1:
            donor = parts[0][:2]  #IF DONOR>9 THAN THIS NEEDS TO BE CHANGED
            timepoint = parts[1]
            type = parts[0][2:]
            # sample = parts[2] if len(parts) > 2 else None #This was for the dataset of Valeria
            sample = parts[3] #This is for dataset Zan
            data.append({'ID': item, 'Donor': donor, 'TimePoint': timepoint, 'Type': type, 'Sample': sample})
    files_mask_df = pd.DataFrame(data)

    #outline files
    data_2 = []
    for timepoint in os.listdir(file_path_outline_files):
        for item in os.listdir(f"{file_path_outline_files}/{timepoint}/FQ_outline"):
            parts = item.split('_')
            if len(parts)>1:
                donor = parts[0][:2]
                type = parts[0][2:]
                #sample = parts[1] if len(parts) > 2 else None #This was for the dataset of Valeria
                sample = parts[3] #This is for data set Zan
                data_2.append({'ID': item, 'Donor': donor, 'TimePoint': timepoint, 'Type': type, 'Sample': sample})
    files_outline_df = pd.DataFrame(data_2)

    #merge the two dataframes
    name_overlap_df = pd.merge(files_mask_df, files_outline_df, on=["Donor", "TimePoint", "Type", "Sample"], suffixes=('_cp', '_outline'))

    return name_overlap_df

def data_prep_cyto(file_path_cellpose, file_path_outline_files):
    #this code is made to find the overlap between the cytoplasm cellpose file and the outline file path (only the names thus!!!)
    #if files are differently saved than explained in the GITHUB, change this code (probably change the for loop structure)

    #cellpose
    data = []
    for folder in os.listdir(file_path_cellpose):
        for item in os.listdir(f'{file_path_cellpose}/{folder}'):
            parts = item.split('_')
            if len(parts) > 1:
                donor = parts[0][:2]
                timepoint = folder
                type = parts[0][2:]
                #sample = parts[1] if len(parts) > 2 else None #This was for the dataset of Valeria
                sample = parts[3] #This is for dataset Zan
                data.append({'ID': item, 'Donor': donor, 'TimePoint': timepoint, 'Type': type, 'Sample': sample})
    files_mask_df = pd.DataFrame(data)

    #outline_files
    data_2 = []
    for timepoint in os.listdir(file_path_outline_files):
        for item in os.listdir(f"{file_path_outline_files}/{timepoint}/FQ_outline"):
            parts = item.split('_')
            if len(parts)>1:
                donor = parts[0][:2]
                type = parts[0][2:]
                #sample = parts[1] if len(parts) > 2 else None #This was for the dataset of Valeria
                sample = parts[3] #This is for data set Zan
                data_2.append({'ID': item, 'Donor': donor, 'TimePoint': timepoint, 'Type': type, 'Sample': sample})
    files_outline_df = pd.DataFrame(data_2)

    # merge the two dataframes
    name_overlap_df = pd.merge(files_mask_df, files_outline_df, on=["Donor", "TimePoint", "Type", "Sample"], suffixes=('_cp', '_outline'))
    return name_overlap_df



def find_belonging_data_set(data_set, name, cell_name):
    #pre_name = '_'.join(name.split('_')[:3]) #this was for data set valeria
    pre_name = '_'.join(name.split('_')[:4]) #this is for data set ZAN
    NAME = pre_name + '_' + cell_name
    data = data_set[data_set.ID == NAME]
    return data, NAME


if __name__ == '__main__':
    #input parameters
    outline_file_path= "F:/Eva/TryOuts/D1C_01_CY35_outline.txt"
    file_name = "F:/Eva/TryOuts/spots_localization_D1C_3H_01.csv"
    file_path_tiff = "F:/Eva/TryOuts/D1C_3H_01_22_38_cp_masks.tif"
    nucleus = True
    cytoplasm = False
    resolution_z = 65
    resolution_xy = 200
    name = "D1C_3H_01_22_38_cp_masks"
    cell_name = 'Cell_CP_32'

    data_set = pre_process_data_set(file_name, nucleus, cytoplasm, resolution_z, resolution_xy)
    data = find_belonging_data_set(data_set, name, cell_name)
    outline_df = make_outline_dataset(outline_file_path)
    cell_name_df = find_cell_name(file_path_tiff, outline_df)

    print(cell_name_df)

