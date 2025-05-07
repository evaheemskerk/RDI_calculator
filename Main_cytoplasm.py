import pandas as pd
import numpy as np
from tqdm import tqdm
from process_data_sets import pre_process_data_set, make_outline_dataset, make_outline_dataset_cyto, find_cell_name, find_belonging_data_set, data_prep, data_prep_cyto
from measurements import polarization_index, dispersion_index, peripheral_distribution_index
from plotting import plot_RDI

def cytoplasm_RDI(file_name, nucleus, cytoplasm, resolution_z, resolution_xy, file_path_outline, file_path_cellpose_n, file_path_cellpose_c, cytokine_analysis):
    name_overlap_n = data_prep(file_path_cellpose_n, file_path_outline)
    name_overlap_c = data_prep_cyto(file_path_cellpose_c, file_path_outline)
    name_overlap_df = pd.merge(name_overlap_n, name_overlap_c, on=["Donor", "TimePoint", "Type", "Sample"], suffixes=('_n', '_c'))
    index_data = pd.DataFrame()
    for i in tqdm(range(len(name_overlap_df))):
        try:
            timepoint = name_overlap_df.loc[i, 'TimePoint']

            outline_file_path = f"{file_path_outline}/{timepoint}/FQ_outline/{name_overlap_df.loc[i, 'ID_outline_n']}"
            file_path_tiff_n = f"{file_path_cellpose_n}/{name_overlap_df.loc[i, 'ID_cp_n']}"
            file_path_tiff_c = f"{file_path_cellpose_c}/{timepoint}/{name_overlap_df.loc[i, 'ID_cp_c']}"

            #find the spots
            data_spots = pre_process_data_set(file_name, nucleus, cytoplasm, resolution_z, resolution_xy, cytokine_analysis)
            outline_cyto = make_outline_dataset_cyto(outline_file_path)
            cell_name_cyto, image_cyto = find_cell_name(file_path_tiff_c, outline_cyto)

            outline_nucleus = make_outline_dataset(outline_file_path)
            cell_name_nucleus, image_nucleus = find_cell_name(file_path_tiff_n, outline_nucleus)

            cell_name_combined = pd.merge(cell_name_nucleus, cell_name_cyto, on='Cell_Name', suffixes=('_nucleus', '_cyto'))


            #find cut stack
            parts = name_overlap_df.loc[i, 'ID_cp_n'].split('_')
            #cut_stack_nucleus = int(parts[3]) # THIS WAS FOR DATASET VALERIA
            cut_stack_nucleus = int(parts[5])  # THIS is FOR DATASET Zan
            parts2 = name_overlap_df.loc[i, 'ID_cp_c'].split('_')
            #cut_stack_cyto = int(parts2[3]) #THIS WAS FOR DATASET VALERIA
            cut_stack_cyto = int(parts2[5])  # # THIS is FOR DATASET Zan

            for cellname in cell_name_combined['Cell_Name'].unique():
                data, NAME = find_belonging_data_set(data_spots, name_overlap_df.loc[i, 'ID_cp_n'], cellname)

                particles = data.iloc[:, [3, 4, 5]].values
                cytokine = data['cytokine'].unique()
                if len(particles) > 1:
                    ##PI & DI
                    intensity_cyto = cell_name_combined.loc[cell_name_combined["Cell_Name"] == cellname, "Intensity_cyto"].values[0]
                    pixels = np.argwhere(np.isin(image_cyto, intensity_cyto))
                    pixels[:, 0] = pixels[:, 0] * resolution_z + (cut_stack_cyto * resolution_z)
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

        except:
            print('Error in code, for image:', name_overlap_df.loc[i, 'ID_cp_n'])
            print('Code is continued without this image')

    index_data['Cytokines'] = index_data['Cytokines'].astype(str)
    index_data['Type'] = index_data['Type'].astype(str)

    return index_data

if __name__ == '__main__':
    #parameters
    nucleus = False
    cytoplasm = True
    resolution_z = 200
    resolution_xy = 65

    # Fill in BOTH if interested in both TNF and IFN, or fill in 'TNF' or 'IFN' when interested in a single analysis
    cytokine_analysis = 'TNF'

    #filepaths VAL
    # file_name = "F:/Eva/TryOuts/Spots_with_localization.csv"
    # file_path_cellpose_n = "F:/Eva/HuRKO_mask/Cellpose_DONE"
    # file_path_cellpose_c = "F:/Eva/CY5_HuRKO/CellPose2"
    # file_path_outline_files = "F:/Eva/Outline_files"
    # file_path_output = 'RDI_calculator_cytoplasm_TNF.xlsx'

    #filepaths ZAN
    file_name = "F:/Eva/Data_Zan/Spots_with_localization.csv"
    file_path_cellpose_n = "F:/Eva/Data_Zan/Zan_nucleus_data"
    file_path_cellpose_c = "F:/Eva/Data_Zan/Cytoplasm_images/CellPose"
    file_path_outline_files = "F:/Eva/Data_Zan/Outline_files"
    file_path_output = 'F:/Eva/Data_Zan/Results/RDI/RDI_calculator_cytoplasm.xlsx'

    #Call main loop and save outcome
    index_data = cytoplasm_RDI(file_name, nucleus, cytoplasm, resolution_z, resolution_xy, file_path_outline_files, file_path_cellpose_n, file_path_cellpose_c, cytokine_analysis)
    index_data.to_excel(file_path_output, index=False)

    #plotting the data
    plot_RDI(index_data)
