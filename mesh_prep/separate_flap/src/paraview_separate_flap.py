import os
import argparse
import glob
import shutil
from paraview.simple import *
import inspect

from paraview_utils import get_array_names, print_header, load_vtp_and_get_arrays, add_suffix
from paraview_calc import threshold_face_id, add_array, get_centroid
import paraview.servermanager as sm 
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import numpy as np
from vtk.numpy_interface import dataset_adapter as dsa

from argparse import RawTextHelpFormatter
from scipy.spatial import KDTree

#Read vtu mesh, input can be a folder with a VTU file or a VTU file
def read_vtu_mesh(file_or_folder_path,reg_name,dont_copy_old_mesh=False):
    meshes_to_rm = {}
    # Check if the input is a folder or a VTU file
    if os.path.isdir(file_or_folder_path):
        # If it's a directory
        if not os.path.exists(file_or_folder_path) or not os.path.isdir(file_or_folder_path):
            raise FileNotFoundError(f"Mesh folder '{file_or_folder_path}' not found or is not a directory.")

        # Get a list of all VTU files in the 'mesh' folder
        vtu_files = [f for f in os.listdir(file_or_folder_path) if f.endswith('.vtu')]
        
        if not vtu_files:
            raise FileNotFoundError(f"No VTU files found in '{file_or_folder_path}'.")
        else:
            print(f"{len(vtu_files)} vtu file(s) found in {file_or_folder_path}:\n {vtu_files}.")


        vtu_files_to_rm = [x for x in vtu_files if x.startswith("mesh_to_remove")]
        vtu_files = [x for x in vtu_files if not x.startswith("mesh_to_remove")]
        # Load the first VTU file using Paraview's SimpleXMLReader
        vtu_file_path = os.path.join(file_or_folder_path, vtu_files[0])

        print(f"[Main solid] Reading {vtu_file_path}")
        mesh = XMLUnstructuredGridReader(registrationName=reg_name,FileName=[os.path.abspath(vtu_file_path)])
        
        #Regions that have to be part of the wall and might need to be removed from the flap 
        #can be provided as "mesh_to_remove_SUFFIX.vtu"
        for vtu_file_path_i in vtu_files_to_rm: 
            print(f"[Mesh to rm] Reading {vtu_file_path_i}")
            reg_name_i = os.path.basename(vtu_file_path_i).replace("mesh_to_remove","to_rm").replace(".vtu","")
            meshes_to_rm[reg_name_i] = XMLUnstructuredGridReader(registrationName=reg_name_i,FileName=[os.path.abspath(os.path.join(file_or_folder_path, vtu_file_path_i))])
        print(f"{len(meshes_to_rm)} files found to change from flap to wall")

        #Copy old file
        if not dont_copy_old_mesh:
            vtu_file_path_old = vtu_file_path + "_old0"
            version = 0
            while os.path.exists(vtu_file_path_old):
                version += 1
                vtu_file_path_old = vtu_file_path_old[:-1] + str(version)
            shutil.copy2(vtu_file_path, vtu_file_path_old)

    elif os.path.isfile(file_or_folder_path) and file_or_folder_path.endswith('.vtu'):
        # If it's a single VTU file directly
        mesh = XMLUnstructuredGridReader(registrationName=reg_name,FileName=[os.path.abspath(file_or_folder_path)])

    else:
        raise ValueError(f"Invalid input {file_or_folder_path}. Please provide a folder containing a mesh folder with VTU files inside or a direct path to a VTU file.")

    # # Optionally, visualize the loaded mesh (you can add more visualization steps here)
    # Show(mesh)
    # Render()
    mesh_vtk = paraview.servermanager.Fetch(mesh)
    print(f"--{mesh_vtk.GetNumberOfCells()} elements found")
    point_data_array_names, _ = get_array_names(mesh_vtk,data_type="PointData",print_all=True)
    cell_data_array_names, _ = get_array_names(mesh_vtk,data_type="CellData",print_all=True)

    return mesh, point_data_array_names, cell_data_array_names, vtu_file_path, meshes_to_rm

#Knn using scipy
def knn_with_labels(points, labels, test_points, k=3):
    tree = KDTree(points)
    distances, nearest_ind = tree.query(test_points, k=k)
    nearest_labels = labels[nearest_ind]
    return nearest_labels, distances, nearest_ind

#Add DOMAIN_ID array and set all values to 2 
def add_default_domain(obj_src): 
    calculator = add_array(obj_src,reg_name='DOMAIN_ID',res_name='DOMAIN_ID',func='2',res_type='Int')
    return calculator

#Script for Programmable Filter
#Identify the flap as domain 3
#Only send the inside of the function to the filter
def assign_flap_domain(inputs):
    mesh_solid = inputs[0]
    mesh_solid_center =  inputs[1]
    wall_outer = inputs[2]
    wall_inner = inputs[3]

    # print(type(mesh_solid))
    try:
        get_array_names(mesh_solid,data_type="CellData",print_all=True)
    except: 
        mesh_solid = dsa.WrapDataObject(sm.Fetch(inputs[0]))
        mesh_solid_center = dsa.WrapDataObject(sm.Fetch(inputs[1]))
        wall_outer = dsa.WrapDataObject(sm.Fetch(inputs[2]))
        wall_inner = dsa.WrapDataObject(sm.Fetch(inputs[3]))
    get_array_names(mesh_solid_center,data_type="PointData",print_all=True)

    nb_neighbors = wall_outer.GetCellData()["nb_neighbors"][0]
    print(f"Using {nb_neighbors} neighbor cells to find solid.")
    mesh_solid_center_points = mesh_solid_center.GetPoints()
    wall_outer_points = wall_outer.GetPoints()
    wall_inner_points = wall_inner.GetPoints()

    #Prepare arrays for knn
    #All fixed cells are the ones from the outer all or the inner wall 
    #Test cells are all the cells in the solid 
    points = np.concatenate([wall_outer_points,wall_inner_points])
    labels = np.zeros(points.shape[0]) + 2 #Outer wall points have WALL label for now
    labels[wall_outer_points.shape[0]:] = 3 #Inner wall points have FLAP label for now
    test_points =  mesh_solid_center_points

    print(f"[OUTER] Nb labels = 2: {np.sum(labels==2)}")
    print(f"[INNER] Nb labels = 3: {np.sum(labels==3)}")
    print("\nRunning Knn...")
    nearest_labels, distances, nearest_ind = knn_with_labels(points, labels, test_points, k=nb_neighbors)
    # plot_3d_points(points, labels, test_points, nearest_labels,nearest_ind,filename="3d.png")
    #Solid must be correct since any cell that is the closest to the outer wall mush be part of the wall
    domain_id_arr = nearest_labels[:,0] 

    #If the outer wall is within nb_neighbors cells, assign it as the wall
    nb_outer_neighbors = nearest_labels.copy()
    nb_outer_neighbors[nb_outer_neighbors != 2] = 0
    nb_outer_neighbors = np.sum(nb_outer_neighbors,axis=1)/2
    domain_id_arr[nb_outer_neighbors >= 1] = 2 
    
    np_array_modified = domain_id_arr.copy()

    print(f"Nb cells in solid: {len(mesh_solid_center.GetPoints())}")
    print(f"Nb of points in wall: {np.sum(np_array_modified==2)}")
    print(f"Nb of points in flap: {np.sum(np_array_modified==3)}")

    domain_id_arr = mesh_solid.CellData["DOMAIN_ID"]

    print(f"[Original] Unique values in DOMAIN_ID: {np.unique(domain_id_arr)} | length: {len(domain_id_arr)}")
    print(f"[Updated] Unique values in DOMAIN_ID: {np.unique(np_array_modified)} | length: {len(np_array_modified)}")

    output.CellData.append(np_array_modified,"DOMAIN_ID")
    print("[assign_flap_domain] VTK array modified successfully and replaced.")

#Add a programmable filter
def add_programmable_filter(input_arrays,script,reg_name='ProgrammableFilter',copy_arrays=1):
    programmable_filter = ProgrammableFilter(registrationName=reg_name,Input=input_arrays)
    programmable_filter.CopyArrays = copy_arrays
    programmable_filter.Script = script
    programmable_filter.UpdatePipeline()

    return programmable_filter 

def clip_src(src,reg_name,clip_type='Box',clip_pos_orig=[-56.39704513549805, -245.451416015625, -200.0],clip_len_norm=[140.71016311645508, 129.28773498535156, 20.0]):
    # create a new 'Clip'
    clip_src = Clip(registrationName=reg_name, Input=src)
    clip_src.ClipType = clip_type
    # clip_src.HyperTreeGridClipper = 'Plane'
    # clip_src.Scalars = ['CELLS', 'ModelRegionID']
    # clip_src.Value = 1.0
    # clip_src.HyperTreeGridClipper.Origin = [13.958036422729492, -180.80754852294922, -250.37559700012207]

    if clip_type == 'Box':
        clip_src.ClipType.Position = clip_pos_orig
        clip_src.ClipType.Length = clip_len_norm 
    elif clip_type == 'Plane':
        clip_src.ClipType.Origin = clip_pos_orig
        clip_src.ClipType.Normal = clip_len_norm
    else:
        print(f"clip_type {clip_type} not recognized")


    return clip_src

#Make sure that the flap is one connected domain 
#Keep the largest connected region 
def keep_largest_domain(mesh_solid): 
    #Get potential flap
    flap = threshold_face_id(mesh_solid,3,"flap",array_name="DOMAIN_ID",cells=True)
    
    #Separate largest connected body of cells (real flap)
    flap_conn = Connectivity(registrationName='Connectivity_flap', Input=flap)

    #Get elements non-connected to the largest body of cells
    non_flap = threshold_face_id(flap_conn,1,"non_flap",array_name="RegionId",cells=False,method="Above Upper Threshold")
    non_flap_vtk = paraview.servermanager.Fetch(non_flap)
    nb_cells_non_flap = non_flap_vtk.GetNumberOfCells()
    print(f"{nb_cells_non_flap} cells found in non-connected flap")
    
    if nb_cells_non_flap > 0:
        #Move non-connected elements to the wall instead of the flap 
        change_cells_to_wall_src = inspect.getsource(change_cells_to_wall).split('\n', 1)[1]
        input_arrays = [mesh_solid,non_flap]

        # change_cells_to_wall(input_arrays); exit()
        mesh_solid = add_programmable_filter(input_arrays=input_arrays,script=change_cells_to_wall_src,reg_name="ProgFilt_change")

    return mesh_solid

#Find elements from the 2nd array in the 1st one and assign DOMAIN_ID 2 to them 
#Should be used as a programmable filter to update values, but can be tested as a function
#Up until the output part since it does not exist
def change_cells_to_wall(inputs):
    mesh_solid = inputs[0]
    non_flap = inputs[1]

    try:
        get_array_names(mesh_solid,data_type="CellData",print_all=True)
    except: 
        mesh_solid = dsa.WrapDataObject(sm.Fetch(inputs[0]))
        non_flap = dsa.WrapDataObject(sm.Fetch(inputs[1]))

    get_array_names(non_flap,data_type="CellData",print_all=True)
    get_array_names(non_flap,data_type="PointData",print_all=True)

    if non_flap.GetNumberOfCells() > 0:
        non_flap_ids = set(non_flap.CellData["GlobalElementID"])
        wall_idx = [ i for i,elem_id in enumerate(mesh_solid.CellData["GlobalElementID"]) if elem_id in non_flap_ids]
        print(f"Nb of points to change from flap to wall: {len(wall_idx)}")
        domain_id_arr = mesh_solid.CellData["DOMAIN_ID"]
        np_array_modified = domain_id_arr.copy()
        np_array_modified[wall_idx] = 2

        output.CellData.append(np_array_modified,"DOMAIN_ID")
        print("[Change_cells_to_wall] VTK array modified successfully and replaced.")
    else:
        print("No element found in non_flap array.")

#TODO: merge this one and the one above
#Find elements from the 2nd array in the 1st one and assign DOMAIN_ID 3 to them 
#Should be used as a programmable filter to update values, but can be tested as a function
#Up until the output part since it does not exist
def change_cells_to_flap(inputs):
    mesh_solid = inputs[0]
    flap = inputs[1]

    try:
        get_array_names(mesh_solid,data_type="CellData",print_all=True)
    except: 
        mesh_solid = dsa.WrapDataObject(sm.Fetch(inputs[0]))
        flap = dsa.WrapDataObject(sm.Fetch(inputs[1]))

    get_array_names(flap,data_type="CellData",print_all=True)
    get_array_names(flap,data_type="PointData",print_all=True)

    if flap.GetNumberOfCells() > 0:
        flap_ids = set(flap.CellData["GlobalElementID"])
        wall_idx = [ i for i,elem_id in enumerate(mesh_solid.CellData["GlobalElementID"]) if elem_id in flap_ids]
        print(f"Nb of points to change from wall to flap: {len(wall_idx)}")
        domain_id_arr = mesh_solid.CellData["DOMAIN_ID"]
        np_array_modified = domain_id_arr.copy()
        np_array_modified[wall_idx] = 3

        output.CellData.append(np_array_modified,"DOMAIN_ID")
        print("[Change_cells_to_flap] VTK array modified successfully and replaced.")
    else:
        print("No element found in flap array.")

#Clip a small part of the mesh for an example 
def clip_all_for_example(mesh_solid,solid_surfaces,clip_solid_pos=[],clip_solid_len=[200,200,10]):
    if clip_solid_pos == []:
        centroid_coords, _ = get_centroid(mesh_solid)
        clip_solid_pos = [centroid_coords[0]-80,centroid_coords[1]-80,centroid_coords[2]-30]

    mesh_solid = clip_src(mesh_solid,"clip_solid",clip_type="Box",clip_pos_orig=clip_solid_pos,
        clip_len_norm=clip_solid_len)
    solid_surfaces["wall_outer"] = clip_src(solid_surfaces["wall_outer"],"clip_wall_outer",clip_type="Box",clip_pos_orig=clip_solid_pos,
        clip_len_norm=clip_solid_len)
    solid_surfaces["wall_inner"] = clip_src(solid_surfaces["wall_inner"],"clip_wall_inner",clip_type="Box",clip_pos_orig=clip_solid_pos,
        clip_len_norm=clip_solid_len)
    return mesh_solid, solid_surfaces

#Add DOMAIN_ID to separate the flap and the wall
def add_domain(mesh_solid,new_solid_path,solid_surfaces,meshes_to_rm,nb_neighbors=45,
    clip_solid=False,clip_solid_orig=[0,0,0],clip_solid_norm=[0,0,1],save_state=False): 
    # mesh_solid, solid_surfaces = clip_all_for_example(mesh_solid,solid_surfaces)

    #(1) Add default DOMAIN_ID = 2 (wall) and (2) pass nb_neighbors as array (Might be a better way to do it)
    mesh_solid_final = add_default_domain(mesh_solid)
    solid_surfaces["wall_outer"] = add_array(solid_surfaces["wall_outer"],reg_name="wall_outer_n",res_name="nb_neighbors",func=nb_neighbors,res_type='Int')
    
    #Clip solid to only the part that needs to be processed (i.e, other part has no flap)
    if clip_solid: 
        print("Clipping solid")
        mesh_solid = clip_src(mesh_solid_final,"clip_solid",clip_type="Plane",clip_pos_orig=clip_solid_orig,
        clip_len_norm=clip_solid_norm)
    else: 
        mesh_solid = mesh_solid_final

    #Get center of solid cells
    mesh_solid_center = CellCenters(registrationName='CellCenters_solid', Input=mesh_solid)    
    mesh_solid_center.VertexCells = 1

    #=======================
    # PREP
    #=======================
    #Read solid and get info
    mesh_solid_vtk = paraview.servermanager.Fetch(mesh_solid)
    outer_wall_vtk = paraview.servermanager.Fetch(solid_surfaces["wall_outer"])
    inner_wall_vtk = paraview.servermanager.Fetch(solid_surfaces["wall_inner"])

    cell_data_array_names, _ = get_array_names(mesh_solid_vtk,data_type="CellData",print_all=True)
    nb_cells_solid = mesh_solid_vtk.GetNumberOfCells()
    nb_cells_outer_wall = outer_wall_vtk.GetNumberOfCells()
    nb_cells_inner_wall = inner_wall_vtk.GetNumberOfCells()
    nb_points_outer_wall = outer_wall_vtk.GetNumberOfPoints()
    nb_points_inner_wall = inner_wall_vtk.GetNumberOfPoints()

    print(f"--{nb_cells_solid} elements found in solid")
    print(f"--{nb_cells_outer_wall} elements and {nb_points_outer_wall} points found in solid - outer wall")
    print(f"--{nb_cells_inner_wall} elements and {nb_points_inner_wall} points  found in solid - inner wall")

    #=======================
    # PROCESS
    #=======================
    #----------------------------------------------------------------
    # Use ProgrammableFilter to replace DOMAIN_ID
    print_header("[ASSIGN FLAP DOMAIN]")
    assign_flap_domain_src = inspect.getsource(assign_flap_domain).split('\n', 1)[1]
    input_arrays = [mesh_solid,mesh_solid_center,solid_surfaces["wall_outer"],solid_surfaces["wall_inner"]]

    # assign_flap_domain(input_arrays); exit()
    mesh_solid = add_programmable_filter(input_arrays=input_arrays,script=assign_flap_domain_src)

    #----------------------------------------------------------------
    #Remove parts that might be in the flap but should be in the wall as specified by the user 
    print_header("[REMOVE WALL CELLS PROVIDED BY USER AS VTU]")
    print(f"{len(meshes_to_rm)} regions to change from flap to wall")
    for k, mesh_to_rm in meshes_to_rm.items(): 
        print(f"[FLAP_TO_WALL] {k}")
        change_cells_to_wall_src = inspect.getsource(change_cells_to_wall).split('\n', 1)[1]
        input_arrays = [mesh_solid,mesh_to_rm]
        mesh_solid = add_programmable_filter(input_arrays=input_arrays,script=change_cells_to_wall_src,reg_name="ProgFilt_change_to_wall_"+k)
    

    #----------------------------------------------------------------
    #Make sure that the flap is one connected domain, keep the largest connected region 
    print_header("[KEEP LARGEST CONNECTED DOMAIN AS THE FLAP]")
    mesh_solid = keep_largest_domain(mesh_solid)

    #----------------------------------------------------------------
    #Copy the final flap to the original mesh
    print_header("[COPY FLAP TO ORIGINAL MESH]")
    flap = threshold_face_id(mesh_solid,3,"flap_conn",array_name="DOMAIN_ID",cells=True)
    change_cells_to_flap_src = inspect.getsource(change_cells_to_flap).split('\n', 1)[1]
    input_arrays = [mesh_solid_final,flap]

    # change_cells_to_flap(input_arrays); exit()
    mesh_solid = add_programmable_filter(input_arrays=input_arrays,script=change_cells_to_flap_src,reg_name="ProgFilt_change_to_flap")

    #Print info
    domain_id = vtk_to_numpy(sm.Fetch(mesh_solid).GetCellData().GetArray('DOMAIN_ID'))
    print(np.unique(domain_id))

    #=======================
    # SAVE
    #=======================
    if save_state:
        print_header("Saving State File")
        path_state = os.path.join(new_solid_path.replace(".vtu","_FlowState.pvsm"))
        SaveState(path_state)
    
    #Save vtp file
    print("Saving {}".format(new_solid_path))
    SaveData(new_solid_path, mesh_solid)
    print_header("[TEST]")
    print("Loading {}".format(new_solid_path))
    src_obj, _, _ = load_vtp_and_get_arrays(new_solid_path,vtu=True)
    domain_id = vtk_to_numpy(sm.Fetch(src_obj).GetCellData().GetArray('DOMAIN_ID'))
    print(np.unique(domain_id)) 

#Use a config file ideally
def get_patient_info(solid_path,patient=None,nb_neighbors=45):
    clip_solid_p2 = {"nb_neighbors":70,
        "clip_solid":True,
        "clip_solid_orig": [28.372434912856573, -189.36982169050356, -192.3585782482902],
        "clip_solid_norm":[-0.29355676595892394, -0.9317686687336415, 0.21361547960312968]}
    
    clip_solid_p3F = {"nb_neighbors":45,
        "clip_solid":True,
        "clip_solid_orig": [14.047147750854492, -180.42921829223633, -68.0],
        "clip_solid_norm":[0.0, 0.0, 1.0]}
    
    clip_solid_false = {"nb_neighbors":nb_neighbors,"clip_solid":False, "clip_solid_orig": None, "clip_solid_norm":None}
    patients_info = {"2":clip_solid_p2,"3F":clip_solid_p3F,"-1":clip_solid_false}

    if patient is None:
        if "p2" in solid_path: 
            patient = "2"
        elif "p3F" in solid_path: 
            patient = "3F"

    if patient not in ["2","3F"]:
        patient = "-1"

    return patients_info[patient]


    
#Read a surface Mesh (.vtp)
def read_surface_mesh(surface_wall_path,name): 
    if surface_wall_path is None or not os.path.exists(surface_wall_path):
        err = f"{name} Surface Mesh Path: ({surface_wall_path}) not found. Cannot extract flap without {name} surface mesh (.vtp)."
        print(err)
        raise FileNotFoundError(err)
        surface_src = None
    else: 
        print(f"{name} Surface Mesh Path: {surface_wall_path}")
        surface_src = XMLPolyDataReader(registrationName=f'{name.replace(" ","_")}.vtp', FileName=[os.path.abspath(surface_wall_path)])
    return surface_src 
     
if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Identify the FLAP in the solid mesh and separate it from the WALL using the solid and the outer/inner wall surface meshes.\n' \
            '1) The DOMAIN_ID array is overwritten (or created) with ID 2 for the WALL and ID 3 for the FLAP. \n' \
            '   "_old" will be appended to the original solid mesh.\n' \
            '2) Regions provided in the solid folder as mesh_to_remove_*vtu will be fixed as the WALL.\n' \
            '3) The function get_patient_info() can be updated to provide a plane where the elements\n' \
            '   above the plane are fixed as the wall.\n\n' \
            'Tips: A) Run the code on the original mesh and visualize the result in Paraview\n' \
            '      B) Adapt the number of neighbors as needed \n' \
            '           --LOWER: The flap will be larger | HIGHER: The flap will be smaller  \n' \
            '           --Example : nb_neighbors = 2 | If the OUTER WALL is within 2 cells, the cell is classifed as a WALL cell \n' \
            # '      B) Adapt the number of neighbors as needed (minimum number of neighbor cells from the outer wall needed for a cell to be part of the WALL), \n' \
            '      C) Use (3) if you know where the first entry tear is and some cells are incorrectly classified as the flap\n' \
            '      D) For complex flaps, use Paraview to select incorrect flap regions and save them as mesh_to_remove_*vtu\n' \
            , formatter_class=RawTextHelpFormatter)
    parser.add_argument('file_or_folder_path',type=str, help='(1) Path to a folder with structure\n' \
            '   mesh/*solid*/*vtu\n' \
            '   mesh/*solid*/mesh-surfaces/solid_wall_outer.vtp\n' \
            '   mesh/*solid*/mesh-surfaces/solid_wall_inner.vtp\n' \
            '(OR 2) Path to a folder containing VTU file(s)\n' \
            '(OR 3) Path to a VTU file')
    parser.add_argument('--wall_outer_path', type=str, help='Path to a VTP file representing the outer wall surface mesh. (Not needed if file_or_folder_path expected structure is respected)')
    parser.add_argument('--wall_inner_path', type=str, help='Path to a VTP file representing the inner wall surface mesh. (Not needed if file_or_folder_path expected structure is respected)')
    parser.add_argument('-p', '--patient', type=str, help='Patient ID')
    parser.add_argument('-n', '--nb_neighbors', default=45,type=int, help='If the OUTER WALL is within nb_neighbors cells, the cell is classifed as a WALL cell')
    # parser.add_argument('-n', '--nb_neighbors', default=45,type=int, help='minimum number of neighbor cells from the outer wall needed for a cell to be part of the WALL')
    parser.add_argument('--save_state', action='store_true',help='Save the state for easy visualization (PATH_TO_VTU_SOLID_FlowState.pvsm).')
    parser.add_argument('--dont_copy_old_mesh', action='store_true',help='Will overwrite the old solid mesh and not save a copy.')
    args = parser.parse_args()

    #=======================
    # GET PATHS AND READ SOLID MESH 
    #=======================
    #In case we are given a folder already processed with the solid
    #Structure: 
    #---Solid Folder: FILE_OR_FOLDER_PATH/mesh/*solid*
    #---Solid Folder: FILE_OR_FOLDER_PATH/mesh/*solid*/mesh-surfaces/solid_wall_outer.vtp
    ## ---Fluid Folder: FILE_OR_FOLDER_PATH/mesh/*fluid* (not needed)
    #Otherwise the outer wall surface mesh should be provided separately
    if os.path.isdir(os.path.join(args.file_or_folder_path, 'mesh')):
        file_or_folder_path = glob.glob(os.path.join(args.file_or_folder_path, 'mesh',"*solid*"))[0]
        outer_wall_path = os.path.join(file_or_folder_path,"mesh-surfaces","solid_wall_outer.vtp")
        inner_wall_path = os.path.join(file_or_folder_path,"mesh-surfaces","solid_wall_inner.vtp") #not needed 
        # fluid_file_or_folder_path = glob.glob(os.path.join(args.file_or_folder_path, 'mesh',"*fluid*"))[0]
    else:
        file_or_folder_path = args.file_or_folder_path
        outer_wall_path = args.wall_outer_path
        inner_wall_path = args.inner_wall_path
    print_header("[READ SOLID MESH]")    
    mesh_solid, _, _, solid_path, meshes_to_rm = read_vtu_mesh(file_or_folder_path,"solid",args.dont_copy_old_mesh)

    #=======================
    # READ OUTER WALL AND INNER WALL SURFACE MESHES
    #=======================
    print_header("[READ OUTER WALL AND INNER WALL SURFACE MESHES]")    
    solid_surfaces = {}
    solid_surfaces["wall_outer"] = read_surface_mesh(outer_wall_path,name="Outer Wall")
    solid_surfaces["wall_inner"] = read_surface_mesh(inner_wall_path,name="Inner Wall") 

    #=======================
    # Get stored clipping info
    #=======================
    print_header("[GET stored clipping coordinates based on PATIENT]")
    print("Will remove the part above the plane if clipping is true")    
    patient_info = get_patient_info(solid_path,args.patient,args.nb_neighbors)
    print(patient_info)

    #=======================
    # ADD DOMAIN_ID to separate the flap and the wall
    #=======================
    print_header("[ADD DOMAIN]")    
    add_domain(mesh_solid,solid_path,solid_surfaces,meshes_to_rm,clip_solid=patient_info["clip_solid"],
        clip_solid_orig=patient_info["clip_solid_orig"],clip_solid_norm=patient_info["clip_solid_norm"],
        nb_neighbors=patient_info["nb_neighbors"],save_state=args.save_state)
    print(f"Successfully loaded and processed VTU mesh from '{file_or_folder_path}'")

