
<h1 align="center">
  Processing and Visualization of CFD Results of Type B Aortic Dissection  using SimVascular and Paraview 
  <br>
</h1>

<h4 align="center">This repository contains code that helps run large-scale CFD experiments on patients with Type B Aortic Dissection using SimVascular and Paraview (see <a href="#repository-structure">repository structure</a> below).  </h4>


<p align="center">
  <a href="#repository-structure">Repository Structure</a> •
  <a href="#all-scripts">All Scripts</a> •
  <a href="#environment-set-up">Environment Set-up</a> •
  <a href="#rcr_prep">rcr_prep</a> •
  <a href="#compflow">compflow</a> • 
  <a href="#analysisrange">analysis/range</a> •
  <a href="#analysisparaview">analysis/paraview</a> •
  <a href="#resources">Resources</a>
</p>

<!-- TODO UPDATE -->
![screenshot](github_summary.png)

## Repository Structure
* **rcr_prep**: Scale total resistance and capacitance, then output RCR at each outlet
* **compflow**: C++ wrapper that calculates the pressure and flow rate at the faces using SimVascular
* **analysis**
  - **range**: Comparison of model's pressure VS Patient's blood pressure
  - **paraview**: Visualization + Processing of final CFD results


## All Scripts 
<!-- * Bash scripts (can be run with bash or slurm)
  - ...sh
  - ...sh -->
* **Python scripts** (see examples below)
  - rcr_prep/get_rcr.py
  - analysis/range/range_analyze.py
  - analysis/paraview/paraview_analyze.py
* **C++ code** (see example below)
  - compflow/src/compute_flow.cxx


## Environment Set-up
* Python Version: >=3.8
* [Paraview 5.11](https://www.paraview.org/download/)
* Optional, needed for compflow: [SimVascular - 2022-07-20](https://simvascular.github.io/)
* Optional, needed for compflow: [VTK 8.2](https://vtk.org/download/)
* Code tested in Ubuntu 20.04.6 LTS
```
conda env create --name ENV_NAME --file ad_environment.yaml
conda activate ENV_NAME
```

## rcr_prep
* **Goal**: Write RCR boundary conditions based on outlet flow fractions and R/C scaling factors <!-- * Bash script: rcr_prep/get_rcr.sh -->
* **Python script**: rcr_prep/get_rcr.py (see python file for more details)
* **Input(s)**:
    - Total resistance scale factor, -R (float)
    - Total capacitance scale factor, -C (float)
    - A file with the outlet flow fractions, --flow_fractions_p (.txt)
    - A file with the order in which the outlets will be written to the output file, --outlet_order_p (.txt)
    - A path to an output directory, --out_dir (directory path)
* **Output(s)**:
    - A file with scaled RCR boundary conditions (.txt)

Example:
```
python get_rcr.py -R 0.85 -C 2 --flow_fractions_p ./outlet_flow_fractions_bm.txt --outlet_order_p ./outlet_order_p2.txt --out_dir .

Output:
R0.85-C2.0.txt
```
## compflow
* **Goal**: C++ wrapper to compute the flow and pressure at the faces. svpost from svSolver converts the data to .vtu and .vtp files, but does not run the second part available in the GUI (flow calculations).
* **C++ source code**: compflow/src/compute_flow.cxx (see code for more details)
* **Requirements**: 
    - C++ 17 (for filesystem library)
    - SimVascular (use underlying GUI function)
    - VTK 8.2 (needed by SimVascular)
    - Please follow the directions from the [svSolver](https://github.com/SimVascular/svSolver) to build VTK. 
* **Input(s)**:
    - The output folder of svpost, i.e. a directory containing the CFD results (.vtu and .vtp files), -e (path)
    - The original job directory, i.e. the directory with the mesh-complete folder under it, -j (path)
* **Output(s)**:
    - Text files that have the pressure/flow at each face. These are the same as the ones ouput by the SimVascular GUI (.txt) 
* **Info**: This section requires some familiarity with building code, and assumes that you built the svSolver from SimVascular and VTK. Feel free to reach out if you are having troubling with it.

Building with CMake:
```
1) Update SV_DIR and VTK_DIR in CMakeLists.txt to point to the SimVascular directory and the VTK directory. 
2) Make a build directory and enter it

mkdir build
cd build
   
3) Execute the build
  cmake ..
  make

4) You should now have an executable called compflow.

ls 

Output: CMakeCache.txt  CMakeFiles  cmake_install.cmake  compflow  Makefile

5) You can make a symbolic link to compflow in your bin folder to avoid needing a full path at runtime.
```


Running compflow - Example:
```
FULL_PATH/compflow -e FULL_PATH/test-converted-results -j FULL_PATH/test
```


## analysis/range
* **Goal**: Plot pressure/flow at different faces and show if we are within 5% or 10% of the patient-specific blood pressure (systolic pressure, diastolic pressure, pulse pressure, and mean pressure) <!-- * Bash script: rcr_prep/get_rcr.sh -->
* **Python script**: analysis/range/range_analyze.py (see bash script and python file for more details)
* **Input(s)**:
    - A path to the converted results folder from SimVascular or compflow. The code uses the text files that have the pressure/flow at each face. --converted_res_dir (path)
    - A file with the all patients' blood pressure, --patient_info (.tsv)
    - A path to an output directory, --out_dir (directory path)
    - The number of points before and after the current one to look for a local min/max. You can use the number of files per cycle divided by 2. --nb_points_comp (int)
* **Output(s)**:
    - Plots of the flow and pressure at all faces (.png) and a summary with the absolute and relative differences (.txt) 
* **Info**: Note that the name of the folder --converted_res_dir can be formatted in a specific way so that some parameters can be extracted from the name and not provided as inputs. See code for details.

Example:
```
python range_analyze.py --patient 2 --converted_res_dir ../test_converted_res_dir --patient_info ./patient_info.tsv --out_dir ./analysis_output  --nb_points_comp 50 
```


## analysis/paraview 


* **Goal**: Processing and Visualization of the CFD results of a patient-specific Type B Aortic Dissection model. 
  - **Slicing**:     

    1. Cuts along the centerline of the aorta using given locations or a stride
    2. Splits slices with 2 regions into a large and a small region
    3. Matches regions to the true lumen and the false lumen using a file identifying the false lumen as the small or large region
  - **Flow calculations**: 

    4. Computes at all timesteps the cross-sectional flow rate and pressure for the TL and FL separately 
  - **OSI and TAWSS calculations**: 
  
    5. Computes the OSI and TAWSS
  - **Visualization**:

    6. Shows the fluid model with the true lumen
    7. Shows the OSI and TAWSS
    8. Shows the velocity at all cross-sections at a given time step (i.e. systole) or all time steps
    8. Shows the velocity streamlines colored by velocity magnitude at a given time step (i.e. systole) or all time steps
  - **Screenshots and state file**:

    9. Saves a screenshot of all views 
    10. Writes a state file that can be open in Paraview if more modifications are needed
* **Requirements**: pvpython (Paraview)
* **Python script**: analysis/paraview/paraview_analyze.py (see config file and python file for more details)
* **Input(s)**:
  - A config file, see example.ini, -c or --config (.ini)
    
  Paths to all the following inputs should be specified in the config file.
  - A centerline of the fluid model (can be obtained from SimVascular) (.vtp) **[REQUIRED]**
  - A folder containing the CFD results at each time step in a .vtu format **[REQUIRED]**
  - A true lumen model if it should be added to the model visualization (.vtp) 
  - A file containing the location of the slices to process. It is suggested to use the automatic slicing function, then move the slices of interest to this file. (.tsv)
* **Main Output(s)**:
    - Views described above (.png) and a file with the flow rate and pressure at each cross-section for all timesteps  (.tsv) 
* **Info**: Suggested use for new data: 

  **[Cross-section dependent]** (1) Automatically slice the aorta, (2) Select slices of interest and copy them to the fixed slice file, (3) Visualize the slices, (4) Update the false_lumen_id file accordingly, (5) Process the cross-sections using loop_timestep

  **[Model without cross-sections, TAWSS, OSI, Velocity Streamlines]** (6) Can be run directly after updating the parameters in the config files

Example: 
```
FULL_PATH/pvpython paraview_analyze.py --config example.ini
```

## Resources

* [SimVascular](https://simvascular.github.io/)
* [svSolver code](https://github.com/SimVascular/svSolver) from SimVascular
* [Paraview](https://www.paraview.org/)
* Interactive Slurm Session on Zaratan
```
sinteractive --t 120 -mem=20g --cpus-per-task=8
```
