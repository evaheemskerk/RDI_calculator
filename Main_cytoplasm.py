import pandas as pd
import numpy as np
from tqdm import tqdm
from process_data_sets import pre_process_data_set, make_outline_dataset, make_outline_dataset_cyto, find_cell_name, find_belonging_data_set, data_prep, data_prep_cyto
from measurements import polarization_index, dispersion_index, peripheral_distribution_index
from plotting import plot_RDI

def cytoplasm_RDI(file_path_outline_files, file_name, nucleus, cytoplasm, resolution_z, resolution_xy, file_path_outline, file_path_cellpose_n, file_path_cellpose_c):
    name_overlap_n = data_prep(file_path_cellpose_n, file_path_outline_files)
    name_overlap_c = data_prep_cyto(file_path_cellpose_c, file_path_outline_files)
    name_overlap_df = pd.merge(name_overlap_n, name_overlap_c, on=["Donor", "TimePoint", "Type", "Sample"],
                               suffixes=('_n', '_c'))  # ff kijken of dit goed gaat?

    index_data = pd.DataFrame()
    for i in tqdm(range(len(name_overlap_df))):
        timepoint = name_overlap_df.loc[i, 'TimePoint']

        outline_file_path = f"{file_path_outline}/{timepoint}/FQ_outline/{name_overlap_df.loc[i, 'ID_outline_n']}"
        file_path_tiff_n = f"{file_path_cellpose_n}/{name_overlap_df.loc[i, 'ID_cp_n']}"
        file_path_tiff_c = f"{file_path_cellpose_c}/{timepoint}/{name_overlap_df.loc[i, 'ID_cp_c']}"


        #find the spots
        data_spots = pre_process_data_set(file_name, nucleus, cytoplasm, resolution_z, resolution_xy)

        outline_cyto = make_outline_dataset_cyto(outline_file_path)
        cell_name_cyto, image_cyto = find_cell_name(file_path_tiff_c, outline_cyto)

        outline_nucleus = make_outline_dataset(outline_file_path)
        cell_name_nucleus, image_nucleus = find_cell_name(file_path_tiff_n, outline_nucleus)

        cell_name_combined = pd.merge(cell_name_nucleus, cell_name_cyto, on='Cell_Name', suffixes=('_nucleus', '_cyto'))

        #find cut stack
        parts = name_overlap_df.loc[i, 'ID_cp_n'].split('_')
        cut_stack_nucleus = int(parts[3])
        parts2 = name_overlap_df.loc[i, 'ID_cp_c'].split('_')
        cut_stack_cyto = int(parts2[3])

        for cellname in cell_name_combined['Cell_Name'].unique():
            data, NAME = find_belonging_data_set(data_spots, name_overlap_df.loc[i, 'ID_cp_n'], cellname)

            particles = data.iloc[:, [3, 4, 5]].values
            cytokine = data['cytokine'].unique()
            if len(particles) > 1:
                ##PI & DI
                intensity_cyto = cell_name_combined.loc[cell_name_combined["Cell_Name"] == cellname, "Intensity_cyto"].values[0]
                pixels = np.argwhere(np.isin(image_cyto, intensity_cyto))
                pixels[:, 0] = pixels[:, 0] * resolution_z + (cut_stack_cyto * 200)
                pixels[:, 1] = pixels[:, 1] * resolution_xy
                pixels[:, 2] = pixels[:, 2] * resolution_xy


                PI = polarization_index(pixels, particles)
                DI = dispersion_index(pixels, particles)

                ##PDI
                intensity_nucleus = cell_name_combined.loc[cell_name_combined["Cell_Name"] == cellname, "Intensity_nucleus"].values[0]
                pixels_n = np.argwhere(np.isin(image_nucleus, intensity_nucleus))
                pixels_n[:, 0] = pixels_n[:, 0] * resolution_z + (cut_stack_nucleus * 200)
                pixels_n[:, 1] = pixels_n[:, 1] * resolution_xy
                pixels_n[:, 2] = pixels_n[:, 2] * resolution_xy

                PDI = peripheral_distribution_index(pixels, pixels_n, particles)

                new_row = {'Cell_Name': NAME, 'Donor': name_overlap_df.loc[i, 'Donor'],
                        'TimePoint': name_overlap_df.loc[i, 'TimePoint'],
                        'Type': name_overlap_df.loc[i, 'Type'], 'Cytokines': cytokine, 'PI': PI, 'DI': DI, 'PDI': PDI}
                index_data = pd.concat([index_data, pd.DataFrame([new_row])], ignore_index=True)

    index_data['Cytokines'] = index_data['Cytokines'].astype(str)
    index_data['Type'] = index_data['Type'].astype(str)

    return index_data

if __name__ == '__main__':
    #parameters
    nucleus = False
    cytoplasm = True
    resolution_z = 200
    resolution_xy = 65

    #filepaths
    file_name = "F:/Eva/TryOuts/Spots_with_localization.csv"
    file_path_cellpose_n = "F:/Eva/HuRKO_mask/Cellpose_DONE"
    file_path_cellpose_c = "F:/Eva/CY5_HuRKO/CellPose2"
    file_path_outline_files = "F:/Eva/Outline_files"
    file_path_output = 'RDI_calculator_cytoplasm_TNF.xlsx'

    #Call main loop and save outcome
    index_data = cytoplasm_RDI(file_path_outline_files, file_name, nucleus, cytoplasm, resolution_z, resolution_xy, file_path_outline_files,
                   file_path_cellpose_n, file_path_cellpose_c)    # #save PI_data in excel
    index_data.to_excel(file_path_output, index=False)

    #plotting the data
    plot_RDI(index_data)
