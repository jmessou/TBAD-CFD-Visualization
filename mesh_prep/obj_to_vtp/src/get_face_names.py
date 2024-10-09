import argparse
from argparse import RawTextHelpFormatter
import paraview
import warnings
import os
paraview.compatibility.major = 5
paraview.compatibility.minor = 11

from paraview.simple import *

from paraview_utils import load_vtp_and_get_arrays, print_header,save_screen, add_suffix, read_data, get_array_names
from paraview_calc import threshold_face_id, get_centroid, get_area

import vtk 
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np
import pandas as pd

from convert_txt_to_mdl import convert_df_to_mdl
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()


# Parse arguments
def _parse_args():
    #~/bin/pvpython convert_obj_to_vtp.py
    parser = argparse.ArgumentParser(prog="get_face_names.py",
        description="(1) Read .vtp file and deduce face names from position.\n"
            "(2) Write face names to a file.\n"
            "Run as PATH_TO_PARAVIEW/bin/pvpython get_face_name.py", formatter_class=RawTextHelpFormatter)
    parser.add_argument("-i",'--input_file_vtp', required=True,type=str, help="Input file (.vtp) "
            "with ModelFaceID cell array starting from 1")
    parser.add_argument('--model_type', type=str, default="solid",
            help="Model type (fluid or solid), 13 faces for fluid, 14 faces for solid.")
    parser.add_argument("-o",'--out_dir', type=str, default=None, help="Output directory "
            "where the file with the face names will be written.\n" 
            "Defaults to the input file folder if not provided.")
    parser.add_argument("-t",'--face_name_template', type=str, default="./face_name_template.txt", help="See template for face order. " 
            "Names can be modified, but the order should stay the same. For a fluid, rename the inner wall as wall for example.")

    return parser.parse_args()

#==========================================
#region          HIGH-LEVEL FUNCTIONS
#==========================================

#Deduce Face names from their location
def get_face_names(aorta_src,model_type="solid",out_dir=".",fname_prefix="test",face_name_template_path="face_names_template.txt",view=None):
    
    #Get face names and face info (centroid/nb_points in face/area)
    face_id, face_names_mapping_dic = get_face_id_dic(face_name_template_path)
    face_df, faces_range, bad_face_ids = get_face_info(aorta_src,view=view)
    face_df["z_section"] = None
    face_df["face_type"] = "cap"
    
    print("Dropping bad face ids...")
    face_df = face_df.drop(bad_face_ids)
    # face_df = face_df.drop([13])
    sorted_area_df = face_df.sort_values(by='area',ascending=False)

    # model_type = "fluid"
    if model_type == "fluid":
        assert len(face_df) == 13, "Received {} faces after drop instead of 13".format(len(face_df))
    else: 
        assert len(face_df) == 14, "Received {} faces after drop instead of 14".format(len(face_df))
        

    #Should always be true, otherwise there is an issue somewhere
    # Outer wall area > Inner wall area > Inlet area > All others 
    if model_type == "fluid":
        face_id["wall_inner"] = sorted_area_df.index[0] 
        face_id["wall_outer"] = -1
        inlet_index=1
    else: #TODO: Apparently not always true, find other way 
        face_id["wall_inner"], face_id["wall_outer"] = assign_wall_order(sorted_area_df) 
        face_df = face_df.drop("face",axis=1) #Column used to assign wall
        inlet_index = 2

    face_df.loc[sorted_area_df.index[:inlet_index],"face_type"] = "wall"
    face_id["inlet"] = sorted_area_df.index[inlet_index] 

    #Compute vector from other faces to the inlet 
    inlet = sorted_area_df.iloc[inlet_index]
    sorted_area_df = change_origin(sorted_area_df,inlet,"inlet") 

    #Now that we have the first 3, we remove the walls and break down the other outlets into top,middle,bottom
    sorted_z_df = sorted_area_df[inlet_index:].sort_values(by='coords_z',ascending=False)
    sorted_z_df.loc[sorted_z_df.index[:4],"z_section"] = "top"
    sorted_z_df.loc[sorted_z_df.index[4:8],"z_section"] = "middle"
    sorted_z_df.loc[sorted_z_df.index[8:],"z_section"] = "bottom"

    #Process each section separately 
    face_id = process_top(sorted_z_df[sorted_z_df["z_section"] == "top"],face_id) 
    face_id, celiac = process_middle(sorted_z_df[sorted_z_df["z_section"] == "middle"],face_id) 
    sorted_z_df = change_origin(sorted_z_df,celiac,"celiac")
    face_id = process_bottom(sorted_z_df[sorted_z_df["z_section"] == "bottom"],face_id) 
    
    
    #Rename faces based on given template, make indexes as keys by swapping keys/values,
    #add vectors from inlet/celiac, and add face names to original dataframe 
    face_id_renamed = {face_names_mapping_dic.get(k, k): v for k, v in face_id.items()}
    face_name_dic = {v: k for k, v in face_id_renamed.items()}
    face_df = change_origin(face_df,inlet,"inlet") 
    face_df = change_origin(face_df,celiac,"celiac")
    face_df['face_names'] = [face_name_dic[i] for i in face_df.index]

    #Move face names to first column
    cols = face_df.columns.tolist()
    cols = cols[-1:] + cols[:-1]  # Move the last element to the front
    face_df = face_df[cols]
    print(face_df)

    face_df = face_df.rename_axis("ModelFaceID")
    save_path = os.path.join(out_dir,fname_prefix + "_face_info.csv")
    save_path_txt = save_path.replace("_face_info.csv","_face_id.txt")
    print("Saving face info dataframe: {}".format(save_path))
    print("Saving face id/names txt: {}".format(save_path_txt))
    face_df.to_csv(save_path)
    face_df[["face_names","face_type"]].to_csv(save_path_txt,sep=" ")
    # sorted_area_df["z_section"]

    mdl_template_path = os.path.join(os.path.dirname(face_name_template_path),"template.mdl")
    
    print_header("[Convert to mdl]")
    convert_df_to_mdl(face_df,save_path.replace("_face_info.csv",".mdl"),template_path=mdl_template_path)

    # SaveState('debug_test.pvsm')

#Verify if inner wall and outer wall were correctely assigned 
#Outer wall doesn't always have a larger area 
#Flip them if needed 
def assign_wall_order(sorted_area_df):
    #Assumer outer has largest surface area 
    wall_outer_df = sorted_area_df.iloc[0]
    wall_inner_df = sorted_area_df.iloc[1]
    wall_outer = sorted_area_df.index[0] #Largest
    wall_inner = sorted_area_df.index[1] #2nd largest 
    outer_wall_coords = [wall_outer_df["coords_x"],wall_outer_df["coords_y"],wall_outer_df["coords_z"]]
    normal = [0,0,1] 

    #Slice outer wall and inner wall
    outer_slice = make_slice("outer",wall_outer_df["face"],outer_wall_coords,normal)
    inner_slice = make_slice("inner",wall_inner_df["face"],outer_wall_coords,normal)
    
    #Make 2D convex hulls and compute their areas
    outer_hull = Delaunay2D(registrationName='outer_hull', Input=outer_slice)
    inner_hull = Delaunay2D(registrationName='inner_hull', Input=inner_slice)
    outer_area = get_area(outer_hull)
    inner_area = get_area(inner_hull)

    #Switch outer/inner wall if necessary
    if outer_area < inner_area: 
        print(f"[ASSIGN_WALL_ORDER - Slice] outer_area ({outer_area}) < inner_area ({inner_area})")
        print(f"[ASSIGN_WALL_ORDER - Slice] Flipping inner/outer wall (outer wall overall surface area is smaller)")
        wall_inner_tmp = wall_inner
        wall_inner = wall_outer
        wall_outer = wall_inner_tmp 
    else:
        print(f"[ASSIGN_WALL_ORDER - Slice] inner_area ({inner_area}) < outer_area ({outer_area})")

    return wall_inner, wall_outer

#Get nb points, area, and coordinates of each face (ModelFaceID)
def get_face_info(aorta_src,view=None):
    aorta_vtk = paraview.servermanager.Fetch(aorta_src)
    faces_range = aorta_vtk.GetCellData().GetArray("ModelFaceID").GetRange(0)
    nb_faces = faces_range[1]
    
    face_locations = {}
    bad_face_ids = [] 
    area_threshold = 5 #Add as a bad face if area in below this number 
    print("Nb Faces: {}".format(nb_faces))

    #Get info about each face 
    for face_id_i in range(1,int(nb_faces)+1):
        centroid = {}
        print("Face ID: {}".format(face_id_i))
        face_i_src  = threshold_face_id(aorta_src,face_id_i,str(face_id_i))
        centroid["face"] = face_i_src

        face_i_src.UpdatePipeline()
        centroid["nb_points"] = face_i_src.GetDataInformation().GetNumberOfPoints()
        centroid_coords, centroid["area"] = get_centroid(face_i_src)
        centroid["coords_x"], centroid["coords_y"], centroid["coords_z"] = centroid_coords
        # centroid["coords"] = centroid_coords
        face_locations[face_id_i] = centroid

        # print(face_locations[face_id_i])
        print("Nb Points: {} | Area: {:.2f} | Centroid: {}".format(centroid["nb_points"],
                centroid["area"], centroid_coords))

        if centroid["area"] < area_threshold: 
            warnings.warn("Area < {} for face ID {}!".format(area_threshold,face_id_i))
            bad_face_ids += [face_id_i] 

        if test:
            SetActiveSource(face_i_src)
            save_path = "."
            Show(); save_screen(add_suffix(save_path,"_"+str(face_id_i))); Hide()
            
    face_df = pd.DataFrame.from_dict(face_locations, orient='index')
    
    return face_df, faces_range, bad_face_ids

#def 
#endregion         END - HIGH-LEVEL FUNCTIONS
#==========================================

#==========================================
#region          FACE PROCESSING FUNCTIONS
#==========================================
#Might need to make it more thorough
def process_top(sorted_z_df,face_id): 
    sorted_z_df = sorted_z_df[sorted_z_df["inlet_norm"] != 0].sort_values(by='inlet_x')
    face_id["b-trunk"] =  sorted_z_df.index[0]
    face_id["carotid"] =  sorted_z_df.index[1]
    face_id["subclavian"] = sorted_z_df.index[2]

    return face_id 

#Might need to make it more thorough
def process_middle(sorted_z_df,face_id): 
    #Get celiac trunk
    sorted_z_df = sorted_z_df.sort_values(by='coords_z',ascending=False)
    face_id["celiac"] =  sorted_z_df.index[0]
    celiac = sorted_z_df.iloc[0]

    #Get others using celiac trunk
    sorted_z_df = change_origin(sorted_z_df,celiac,"celiac")
    sorted_z_df = sorted_z_df[sorted_z_df["celiac_norm"] != 0].sort_values(by='celiac_x')
    face_id["sup-mes"] = sorted_z_df.index[1]
    face_id["r-renal"] = sorted_z_df.index[0]
    face_id["l-renal"] = sorted_z_df.index[2]

    return face_id, celiac

#Might need to make it more thorough
def process_bottom(sorted_z_df,face_id): 
    sorted_z_df = sorted_z_df.sort_values(by='celiac_y')
    sorted_z_df_ext = sorted_z_df[:2].sort_values(by='celiac_x')
    sorted_z_df_int = sorted_z_df[2:].sort_values(by='celiac_x')

    face_id["r-ext-ill"] =  sorted_z_df_ext.index[0]
    face_id["r-int-ill"] = sorted_z_df_int.index[0]
    face_id["l-int-ill"] = sorted_z_df_int.index[1]
    face_id["l-ext-ill"] = sorted_z_df_ext.index[1]

    return face_id 

#Compute vector from origin to all vectors in df, and add it as new columns
def change_origin(df,origin,col_name): 
    df[col_name+'_x'] = df["coords_x"] - origin["coords_x"]
    df[col_name+'_y'] = df["coords_y"] - origin["coords_y"]
    df[col_name+'_z'] = df["coords_z"] - origin["coords_z"]
    vectors = df[[col_name+'_x',col_name+'_y',col_name+'_z']].values
    df[col_name+'_norm'] = np.linalg.norm(vectors,axis=1)

    return df

#endregion         END - FACE PROCESSING FUNCTIONS
#==========================================

#==========================================
#region          SLICE FUNCTIONS (copied from analysis code)
#==========================================
#Make single slice, use box for ATA to avoid cutting in descending aorta
#TODO improve cutting method instead of shifting box size
def make_slice(slice_name,input_src,origin,normal,box_size=70,z_normal=True):
    if z_normal:
        normal = [0,0,1]

    #TODO: Always use box for now
    input_src = Clip(Input=input_src)
    input_src.ClipType = 'Box'
    input_src.ClipType.Position = [x-box_size/2 for x in origin] #Bottom square set at bottom quarter
    input_src.ClipType.Length = [box_size]*3

    # create a new 'Slice'
    slice_ = Slice(registrationName=slice_name, Input=input_src)
    slice_.SliceType = 'Plane'
    slice_.SliceType.Origin = origin
    slice_.SliceType.Normal = normal
    if slice_.PointData.GetNumberOfArrays() == 0:
        print("--[WARNING!-Make Slice] No point")
        print("    Returning None")
        slice_ = None
        raise ValueError(f"Empty slice: {slice_name} | Origin: {origin} | Normal: {normal} | Box size: {box_size}")
    
    slice_.UpdatePipeline()
    return slice_

#endregion         END - (copied from analysis code)
#==========================================

#==========================================
#region          UTILS FUNCTIONS
#==========================================
#Read face names from template file 
#Names can change, but the order needs to stay the same 
def get_face_id_dic(face_name_template_path="face_names_template.txt"):
    print("Using this template for face names:\n"
        "---{}".format(face_name_template_path))
    keys_in_code = ["wall_outer","wall_inner",
    "inlet","b-trunk","carotid","subclavian",
    "celiac","sup-mes","r-renal","l-renal",
    "r-ext-ill", "r-int-ill", "l-int-ill","l-ext-ill"]
    face_names = read_data(face_name_template_path)
    face_id_dic = {key: None for key in keys_in_code}
    face_name_dic = {key: val for key,val in zip(keys_in_code,face_names)}
    return face_id_dic, face_name_dic

#endregion         END - UTILS FUNCTIONS
#==========================================


#==========================================
#region          VISUALIZATION FUNCTIONS
#==========================================

def set_up_view(coords,suffix="",assign_layout=True,view_dic={"az":-35,"el":0}, analysis_time_step=0):
    # ----------------------------------------------------------------
    # setup views used in the visualization
    # ----------------------------------------------------------------
    # get the material library
    materialLibrary1 = GetMaterialLibrary()
    x_pos = (np.abs(coords[2]["max"]-coords[2]["min"]) + 20)*2

    # Create a new 'Render View'
    renderView1 = CreateView('RenderView')
    # renderView1.ViewSize = [1177, 550]
    renderView1.ViewSize = [550, 550]
    #renderView1.ViewSize = [2300, 1100]
    renderView1.AxesGrid = 'GridAxes3DActor'
    renderView1.CenterOfRotation = [coords[0]["center"], coords[1]["center"], coords[2]["center"]]
    renderView1.StereoType = 'Crystal Eyes'
    renderView1.CameraPosition = [x_pos, coords[1]["center"], coords[2]["center"]]
    renderView1.CameraFocalPoint = [coords[0]["center"], coords[1]["center"], coords[2]["center"]]
    renderView1.CameraViewUp = [0,0, 1]
    renderView1.CameraFocalDisk = 1.0
    #renderView1.CameraParallelScale = 257.3654953693423
    renderView1.BackEnd = 'OSPRay raycaster'
    renderView1.OSPRayMaterialLibrary = materialLibrary1
    
    
    #Select active view
    SetActiveView(renderView1)
    renderView1.ViewTime = analysis_time_step
    GetActiveCamera().Azimuth(view_dic["az"])
    GetActiveCamera().Elevation(view_dic["el"])
    Render()    
    # GetActiveCamera().Azimuth(150)
    if assign_layout:
        # create new layout object 'Layout #1'
        # layout1 = CreateLayout(name='Layout' + suffix)
        layout1 = CreateLayout(name=suffix)
        layout1.AssignView(0, renderView1)
        layout1.SetSize(renderView1.ViewSize)
        #layout1.SetSize(2300, 1100)
    
    return renderView1, layout1

def get_vtk(src_vtp):
    src_vtk = paraview.servermanager.Fetch(src_vtp)
    # print(dir(src_vtk))
    # import pdb; pdb.set_trace()
    # reader = vtk.vtkXMLUnstructuredGridReader()                                  
    # reader.SetFileName(src_vtp.FileName[0])    
    # reader.Update()                                                                                                                                                            
    mesh = src_vtk                           
    pointCoords = vtk_to_numpy(mesh.GetPoints().GetData())
    camera = {"center":0,"diff":0,"max":0,"min":0,"q1":0,"q3":0}
    #coords = {"x":camera.copy(),"y":camera.copy(),"z":camera.copy()}
    coords = {0:camera.copy(),1:camera.copy(),2:camera.copy()}
    for i in range(3):
        coords[i]["max"] =  np.max(pointCoords[:,i])
        coords[i]["min"] = np.min(pointCoords[:,i])
        coords[i]["center"] = (np.max(pointCoords[:,i]) + np.min(pointCoords[:,i]))/2
        coords[i]["diff"] = np.abs(np.max(pointCoords[:,i]) - np.min(pointCoords[:,i]))   
        coords[i]["q1"] = (coords[i]["min"] + coords[i]["center"])/2
        coords[i]["q3"] = (coords[i]["max"] + coords[i]["center"])/2

    # print(coords)
    return coords, pointCoords

#endregion         END - VISUALIZATION FUNCTIONS
#==========================================

# test=True
test=False
if __name__ == '__main__':
    args = _parse_args()

    if args.out_dir is None: 
        args.out_dir = os.path.dirname(args.input_file_vtp)
    input_fname = os.path.basename(args.input_file_vtp)
    src_vtp, _, _ = load_vtp_and_get_arrays(args.input_file_vtp)

    coords, _ = get_vtk(src_vtp) 
    
    if test:
        #Set up View for testing
        renderView1, layout1 = set_up_view(coords,view_dic={"az":0,"el":0},analysis_time_step=0)
        # renderView1, layout1 = set_up_view(coords)
        SetActiveView(renderView1)
        Show(src_vtp); save_screen("./test.png"); Hide()
    else:
        renderView1 = None

    # import pdb; pdb.set_trace()

    get_face_names(src_vtp,model_type=args.model_type,out_dir=args.out_dir,
        fname_prefix=input_fname.replace(".vtp",""), face_name_template_path=args.face_name_template,view=renderView1)

