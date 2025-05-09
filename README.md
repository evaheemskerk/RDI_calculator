# Spatial distribution using RDI-calculator
This algorithm allows you to find the three spatial localization indexes in 3d, proposed by Stueland _et al._(2019) and based on formulas from Wang _et al._ (2019) and Park _et al._ (2012)

Michael Stueland, Tianhong Wang, Hye Yoon Park, and Stavroula Mili. Rdi calculator: an analysis tool to assess rna distributions in cells. Scientific reports, 9(1):8267, 2019. 

Tianhong Wang, Susan Hamilla, Maggie Cam, Helim Aranda-Espinoza, and Stavroula Mili. Extracellular matrix stiffness and cell contractility control rna localization to promote cell migration. Nature communications, 8(1):896, 2017.

Hye Yoon Park, Tatjana Trcek, Amber L Wells, Jeffrey A Chao, and Robert H Singer. An unbiased analysis method to quantify mrna localization reveals its correlation with cell motility. Cell reports, 1(2):179–184, 2012.

# Content
* main_nucleus: run when interested in the nucleus
* main_cytoplasm: run when interested in the cytoplasm 

## Installation
Before using the python algorithm, ensure you have the following libraries installed in your local environment: numpy pandas seaborn tqdm imageio scipy.stats matplotlib.pyplot
If not; use the package manager [pip](https://pip.pypa.io/en/stable/) to install these.

```bash
pip install numpy pandas seaborn tqdm imageio scipy.stats matplotlib.pyplot
```
Download all the python files in this repository. 

## Inputs 
Change the inputs at the bottom of the script. 

Input parameters:
* **nucleus** = Boolean set on True when interested in nucleus 
* **cytoplasm** = Boolean set on True when interested in cytoplasm 
* **resolution_z** = resolutioin in Z-axis
* **resolution_xy** = resolution in XY-direction
* **cytokine_analysis** = Fill in BOTH if interested in both TNF and IFN, or fill in 'TNF' or 'IFN' when interested in a single analysis
  
Input files:
* **file_name** = File path to the Spots_with_localization.csv file
* **file_path_cellpose_n** = File path to the segmentation images of the nucleus
* **file_path_cellpose_c** = File path to the segmentation images of the cytoplasm (does not exist when main_cytoplasm is run)
* **file_path_outline_files** =  File path to the outline files
* **file_path_output** = Output file path

# File handling
Be sure that all outline files are saved as follows: 
```python
file_path_outline_files\{timepoint}\FQ_outline
```
And the files are saved as:
```python
{Donor}{Type}_{Timepoint}__{Sample}
```

Be sure that the nucleus segmentation images files are all saved in one folder and the images are saved as follows: 
```python
{Donor}{Type}_{Timepoint}__{Sample}_{stain}_{firstcut}_{second_cut}
```

Be sure that the cytoplasm segmentation images files are all saved as follows: 
```python
file_path_cellpose_c\{timepoint}
```
And the files are saved as:
```python
{Donor}{Type}_{Timepoint}__{Sample}_{stain}_{firstcut}_{second_cut}
```

# Run
Run the script. 
