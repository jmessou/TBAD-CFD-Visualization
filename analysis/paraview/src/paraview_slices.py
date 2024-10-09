#Utility functions: centerline/loading/saving/printing 

#Depends on paraview calc 

import paraview 
from paraview.simple import *
from paraview.selection import _createSelection, _collectSelectionPorts
from vtk.util.numpy_support import vtk_to_numpy
from vtk.numpy_interface import dataset_adapter as dsa

import paraview.servermanager as sm 

from paraview_calc import *
from paraview_utils import save_screen, add_suffix, _QuerySelect, _GroupTimeSteps, switch_var, make_parent_dir, map_fields, get_array_names

from paraview_display import *  

import pandas as pd 
import numpy as np 
import os 
import glob 
import warnings 
import ast 

#==========================================
#region          HIGH-LEVEL FUNCTIONS
#==========================================

#Load slices df, then make slices
#(see load_slice_info for info on slices_df)
def load_and_make_slices(input_src,slices_tsv,results_tsv="results.tsv",box_size=70,z_normal=True,
        show=True,force_normal=150,fl_id_tsv="",process_slice_bool=True,append_res=False,
        split_by_centroid_bool=True): 
    slices_df = load_slices_info(slices_tsv)
    slices, slice_names, results_df = make_slices(input_src,slices_df,results_tsv=results_tsv,
        box_size=box_size,z_normal=z_normal,show=show,force_normal=force_normal,fl_id_tsv=fl_id_tsv,
        process_slice_bool=process_slice_bool,append_res=append_res,split_by_centroid_bool=split_by_centroid_bool)

    return slices, slice_names, results_df

#Cut along the centerline with a given stride  
def slice_centerline(input_src,centerline,box_size=70,stride=100,slices_tsv="",results_tsv="results.tsv",
        z_normal=True,show=True,force_normal=150,erase_range=[0,0],fl_id_tsv="",process_slice_bool=True,
        append_res=False,split_by_centroid_bool=True):
    #Get origins/normals from centerline and store it as a dataframe
    vtk_centerline = sm.Fetch(centerline)
    slice_names, origins, normals = centerline_to_origin_normals(vtk_centerline,stride=stride)
    slices_df = centerline_to_df(slice_names, origins, normals)
    
    #Drop last slice 
    # first_slice = slices_df.iloc[0].name[0]
    last_slice = slices_df.iloc[-1].name[0]
    # print("Dropping {} and {}".format(first_slice,last_slice))
    # slices_df = slices_df.drop([first_slice,last_slice],axis=0)
    print("Dropping {}".format(last_slice))
    slices_df = slices_df.drop([last_slice],axis=0)

    #Save cross-sections as tsv files
    if slices_tsv != "":
        save_slices_tsv(slices_df,add_suffix(slices_tsv,"_{}".format(len(slices_df))))


    results_tsv = add_suffix(results_tsv,"_{}".format(int(len(slices_df)/2)))
    slices, slice_names, results_df = make_slices(input_src,slices_df["Coord"],
            box_size=box_size,results_tsv=results_tsv,z_normal=z_normal,show=show,
            force_normal=force_normal,fl_id_tsv=fl_id_tsv,process_slice_bool=process_slice_bool,
            append_res=append_res,split_by_centroid_bool=split_by_centroid_bool)
    
    return slices, slice_names, results_df

#endregion         END - HIGH-LEVEL FUNCTIONS
#==========================================

#==========================================
#region          LOADING/SAVING FUNCTIONS
#==========================================

#Load slices from file and return a df so that
# slices_df[slice_name]["Origin"] and slices_df[slice_name]["Normal"]
# are 3D tuples 
# Slice Name (ATA...)  | Type (Origin/Normal) | x | y |z
def load_slices_info(slices_tsv,box_size=70):
    slices_df = pd.read_csv(slices_tsv,sep="\t")
    slices_df["Coord"] = slices_df[["x","y","z"]].apply(tuple,axis=1)
    slices_df = slices_df.drop(["x","y","z"],axis=1)
    slices_df = slices_df.set_index(['Slice Name','Type'])
    return slices_df["Coord"]

#Save slices df to tsv file 
# Slice Name (ATA...)  | Type (Origin/Normal) | x | y |z
def save_slices_tsv(slices_df,slices_tsv):
    #Assume all the origin blocks are together (same for normals)
    #Will write in the order Origin/Normal/Origin/Normal etc... 
    middle = int(len(slices_df)/2)
    block1 = list(range(middle))
    block2 = list(range(middle,len(slices_df)))
    if slices_df.iloc[0].name[1] == "Origin":
        new_index = [item for pair in zip(block1, block2 + [0]) for item in pair]
    else:
        new_index = [item for pair in zip(block2, block1 + [0]) for item in pair]
    
    slices_df = slices_df.copy()
    slices_df['x'], slices_df['y'], slices_df['z'] = slices_df["Coord"].str
    print("Saving slices to {}".format(slices_tsv))
    #Reindex, drop tuple column, and save 
    slices_df.reset_index().reindex(new_index).drop(["Coord"],axis=1).to_csv(slices_tsv,sep="\t",index=False)


#endregion         END - LOADING/SAVING FUNCTIONS
#==========================================

#==========================================
#region          CENTERLINE FUNCTIONS
#==========================================

# Given a VTK centerline, return origin and normal points for a given stride
# Use current and next point for normals
# NOTE: Slicing at the caps can be off since normals are estimated from next point
# You should drop the first and last slice to have accurate results 
# Use cap info instead 
#Use skip_first to drop first point 
def centerline_to_origin_normals(vtk_centerline,stride=60,skip_first=True):
    #Get ID and points
    centerline_id = vtk_centerline.GetCellData().GetArray("CenterlineIds").GetValue(0)
    all_points = vtk_to_numpy(vtk_centerline.GetPoints().GetData())
    print("{} point(s) along centerline with ID {}".format(len(all_points),centerline_id))
    
    if skip_first: 
        all_points = all_points[1:]
    #Normals are from current point to next point, use previous point normal for last slice
    #TODO: Normal might be 0 if points are on the same plane, then need to use next point until normal is not 0
    all_normals = all_points[1:] - all_points[:-1]
    all_normals = np.concatenate([all_points[2:],all_points[-1:]]) - all_points[:-1]
    all_normals = np.concatenate([all_normals,all_normals[-1:]])
    slice_names = list(range(len(all_points)))

    #Keep the ones we are interested in 
    origins = all_points[::stride]
    normals = all_normals[::stride]
    slice_names = slice_names[::stride]

    #Always keep the last point no matter the stride, so we can drop it accurately
    #outside of this function
    if (all_points[-1] != origins[-1]).all():
        origins = np.concatenate([origins,all_points[-1:]])
        normals = np.concatenate([normals,all_normals[-1:]])
        slice_names = slice_names + [len(all_points)-1]

    print("Slicing with stride {} will generate {} slices".format(stride,len(origins)))
    slice_names = ["Slice_{}".format(i) for i in slice_names]

    return slice_names, origins, normals

#Take names, origins, and normals and convert them to a hierachical dataframe 
#                       Coord
# Slice Name | Type 
def centerline_to_df(slice_names, origins, normals):
    origins_df = pd.DataFrame(index=slice_names,data=zip(np_to_tuple(origins)),columns=["Coord"])
    normals_df = pd.DataFrame(index=slice_names,data=zip(np_to_tuple(normals)),columns=["Coord"])
    slices_df = pd.concat([origins_df,normals_df],keys=["Origin","Normal"]).swaplevel()
    slices_df.index.names = ["Slice Name", "Type"]

    return slices_df 

#Change numpy array to 1D list of tuples 
def np_to_tuple(data):
    return list(map(tuple, data.tolist()))

#endregion         END - CENTERLINE FUNCTIONS
#==========================================

#==========================================
#region          SLICES FUNCTIONS
#==========================================

#Generate and process all the given slices based on provided origin/normals
#(see load_slice_info for info on slices_df)
#split_by_centroid: 
#   -When false, will split using large/small regions and then correct the 
#   regions using fl_id_tsv 
#   -When true, uses the position of the centroid in fl_id_tsv to split directly between tl/fl
def make_slices(input_src,slices_df,box_size=70,results_df=None,results_tsv="results.tsv",z_normal=True,
        show=True,force_normal=150,erase_range=[0,0],fl_id_tsv="",process_slice_bool=True,append_res=False,
        split_by_centroid_bool=True,ID_ARRAY="GlobalElementID"):
    #slices = {"ATA":None,"DTA":None,"DTA2":None,"IAA":None}
    slice_names = {name[0]: None for name in slices_df.index.unique().values}
    slices = {}
    keys_to_rm = [] 
    slice_names_fl = [] 
    surface_ids_fl = []
    surface_ids_coords_fl = []
    size_tl = "LARGE" #Default, will be corrected if split_by_centroid_bool is used 
    size_fl = "SMALL" #Default, will be corrected if split_by_centroid_bool is used 
    fl_id = None

    # count = 0
    # disp = Show(input_src)
    # disp.Opacity = 0.05
    # Render()
    if split_by_centroid_bool: 
        _fl_id_tsv = get_fl_id_tsv_full_name(fl_id_tsv,results_tsv)
        small_tl, small_fl = get_fl_ids(_fl_id_tsv,rename_small_by_tl_fl=True)
        if small_tl is None and small_fl is None: 
            print(f"[WARNING!] small_tl ({small_tl}) and small_fl ({small_fl}). Turning off split_by_centroid, will split using LARGE/SMALL")
            split_by_centroid_bool = False

    # import pdb; pdb.set_trace()
    print("-----------------------------------------------")
    print(results_df)
    for slice_name in slice_names.keys(): 
        # if "IAA" not in slice_name and "DTA2" not in slice_name:
        #     continue
        # if "DTA2" not in slice_name:
        #     continue
        if split_by_centroid_bool: 
            slice_name_tl = "TL_" + slice_name
            slice_name_fl = "FL_" + slice_name
        else:
            #Convert later to tl/fl, split is done based on large and small sections
            slice_name_tl = "LARGE_" + slice_name
            slice_name_fl = "SMALL_" + slice_name
        slice_nb = slice_name.split("_")[-1]
        #TODO make it an option
        if slice_nb.isdigit() and float(slice_nb) > force_normal:
            z_normal = True 
        #Remove certain keys TODO MOVE
        if slice_nb.isdigit() and float(slice_nb) < erase_range[1] and float(slice_nb) > erase_range[0]:
            keys_to_rm += [slice_name]
            continue
        print("Cutting: {}".format(slice_name))
        #(1) Cut , (2) Split into TL and FL, (3) Process TL
        slices[slice_name] = make_slice(slice_name,input_src,origin=slices_df[slice_name]["Origin"],normal=slices_df[slice_name]["Normal"],box_size=box_size,z_normal=z_normal)
        #Happens when the point is outside of the aorta
        if slices[slice_name] is None:
            keys_to_rm += [slice_name]
            slices.pop(slice_name)
            print("Empty slice, Popping it!")
            continue 

        if split_by_centroid_bool:
            fl_id = small_fl[ID_ARRAY].get(slice_name_fl)
            fl_id_coords = small_fl[f"{ID_ARRAY} Coordinates"].get(slice_name_fl)
            tl_id = small_tl[ID_ARRAY].get(slice_name_tl)
            tl_id_coords = small_tl[f"{ID_ARRAY} Coordinates"].get(slice_name_tl)
            # slices[slice_name_tl], slices[slice_name_fl], size_tl, size_fl = split_by_id(slices[slice_name],slice_name,
            #         fl_id=fl_id,tl_id=tl_id,fl_id_coords=fl_id_coords,tl_id_coords=tl_id_coords)
            slices[slice_name_tl], slices[slice_name_fl], size_tl, size_fl = split_by_position(slices[slice_name],slice_name,
                fl_id_coords=fl_id_coords,tl_id_coords=tl_id_coords,tolerance=2)
        else: 
            slices[slice_name_tl], slices[slice_name_fl] = split_large_small(slices[slice_name],slice_name)
 
        if process_slice_bool and slices[slice_name_tl] is not None:
            results_df, extractSurface , _ , _ = process_slice(slices[slice_name_tl],slice_name_tl,results_df,slice_size=size_tl)
        
        #If no False Lumen, set empty column
        if slices[slice_name_fl] is not None:
            slice_names_fl += [slice_name_fl]
            surface_id_fl = -1 #Default
            surface_id_coords_fl = -1
            if process_slice_bool:
                results_df, _ , _ , _ = process_slice(slices[slice_name_fl],slice_name_fl,results_df,slice_size=size_fl)
                id_column = f"{ID_ARRAY} "+slice_name_fl
                id_coords_column = f"{ID_ARRAY} Coordinates "+slice_name_fl
                if id_column in results_df.columns:
                    surface_id_fl = results_df[id_column].iloc[-1]
                    surface_id_coords_fl = results_df[id_coords_column].iloc[-1]
            else: #We still record the ID of the cell centroid if it exists
                #No need to extract surface, but leaving it just in case
                extractSurface, empty_slice = extract_surf_and_check_connectivity(slice_name,slices[slice_name_fl])
                if not empty_slice: 
                    _, _, _, surface_id_fl, surface_id_coords_fl = find_closest_cell_to_centroid(slice_name_fl,extractSurface)
                # _, _, _, surface_id_fl, surface_id_coords_fl = find_closest_cell_to_centroid(slice_name_fl,slices[slice_name_fl])

            surface_ids_fl += [surface_id_fl]
            surface_ids_coords_fl += [surface_id_coords_fl]
        else: #Doing it this way for clarity
            if process_slice_bool: 
                results_df[slice_name_fl] = None
            
        print()
        # Show(slices[slice_name])
        # count += 1
        # if count % 10 == 0:
        #     Render()
        #     import pdb; pdb.set_trace()
        # save_screen("centerline.png",image_res=[800,500])

    for key in keys_to_rm:
        slice_names.pop(key)

    # for slice_n in slices.keys():
    #     print(f"{slice_n}: {get_area(slices[slice_n])}")

    #Show all slices and wall
    if show:
        show_slices(slices,slice_names,input_src)
    #Save df results 
    if process_slice_bool:
        print(results_df)
        if append_res:
            print("Reading current tsv ({}) and appending to it...".format(results_tsv))
            print("Note: TL/FL file will be overwritten as well")
            results_df_i = pd.read_csv(results_tsv,sep="\t")
            results_df = pd.concat([results_df_i,results_df],axis=0,ignore_index=True)
            results_df = results_df.drop(["Timestep"],axis=1)
            
        results_df.to_csv(results_tsv,sep="\t",index_label="Timestep",float_format="%.3f")
    # exit()

    # print(fl_id_tsv)
    # print(results_tsv)
    if not split_by_centroid_bool and fl_id_tsv != "": #conversion from large/small to tl/fl not needed otherwise
        fl_id_tsv = get_fl_id_tsv_full_name(fl_id_tsv,results_tsv)
        if not os.path.exists(fl_id_tsv):
            write_fl_id_tsv(fl_id_tsv,surface_ids_fl,surface_ids_coords_fl,slice_names_fl,overwrite=False)
        else:
            slices, slice_names, results_df = process_fl_id_tsv(fl_id_tsv,process_slice_bool, slices, slice_names, results_tsv, results_df)

    # import pdb; pdb.set_trace()
    return slices, slice_names, results_df 

#Make single slice, use box for ATA to avoid cutting in descending aorta
#TODO improve cutting method instead of shifting box size
def make_slice(slice_name,input_src,origin,normal,box_size=70,z_normal=True):
    if z_normal:
        normal = [0,0,1]

    #Clip Volume first to avoid cutting in descending aorta
    #No longer need this
    if False and slice_name == "ATA":
        input_src = Clip(registrationName="Box_ATA", Input=input_src)
        input_src.ClipType = 'Box'
        #TODO Dont' hardcode sizes, 
        input_src.ClipType.Position = [x-box_size/2 for x in origin] #Bottom square set at bottom quarter
        input_src.ClipType.Length = [box_size]*3
       
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

#Split a surface into Large and Small Areas/Lumen
#TODO: FALSE Assumption: Assign the largest area as the True Lumen
#Going with that for now, ideally you want to track it from the beginning
#Post-process 
# FALSE ASSUMPTION ==> Need to read from tsv file I guess...
#like a correction mechanism for now 
#BEST WAY: Just write the values based on areas 
#Write the center point for each integral
#Track center point to see behavior  
#mode_extra = ["exit,"whole_slide","skip" (default)]
#Whole slice sends whole slice as TL (can be used to debug)
def split_large_small(slice_src,slice_name,mode_extra="skip"):
    slice_name_large = "LARGE_" + slice_name
    slice_name_small = "SMALL_" + slice_name
    # region_id_large = 1.0 
    # region_id_small = 0.0 

    if GetActiveView() is not None:
        Show(slice_src)
    #Find connected regions
    connectivity = Connectivity(registrationName='Connectivity_' + slice_name, Input=slice_src)
    regions = sm.Fetch(connectivity).GetCellData().GetArray("RegionId") 
    if regions is not None:
        id_range = regions.GetValueRange()
    else:
        print("--[WARNING!-Connectivity] NO region found")
        print("    Returning given slice as LARGE")
        return slice_src, None
    
    if id_range[1] not in [0,1]:
        print("--[WARNING!-Connectivity] {} regions found".format(id_range[1]+1)) #
        if mode_extra == "exit": 
            print("\n--Exit...")
            exit() 
        elif mode_extra == "whole_slice":
            print("--LARGE = whole slice...")
            large = slice_src; small=None
        else:
            print("--Skip...")
            large = None; small=None
    elif id_range[1] == 0: #No SMALL 
        # print("--No False Lumen")
        large = slice_src 
        small = None 
    else:
        print("--False Lumen found")
        large = threshold_id(connectivity,0,slice_name_large)
        small = threshold_id(connectivity,1,slice_name_small)

        area_large = get_area(large)
        area_small = get_area(small)
        # print("--Area LARGE: {} | Area SMALL: {} | LARGE>SMALL: {}".format(area_large,area_small,area_large>area_small))
        #TODO: Making this assumption for now
        if area_large < area_small: #Switch
            print("[SWITCH] Connectivity 0 is for small")
            print("--Area LARGE: {} | Area SMALL: {} | LARGE>SMALL: {}".format(area_large,area_small,area_large>area_small))
            area_large, area_small = switch_var(area_large,area_small)
            large, small = switch_var(large,small)
        
        print("--Area LARGE: {} | Area SMALL: {} | LARGE>SMALL: {}".format(area_large,area_small,area_large>area_small))
        if GetActiveView() is not None:
            HideAll()
    return large, small

#Split into large/small, then assign TL/FL based on GlobalElementIDs provided 
#If fl_id is found in large, then large, if it is found in small, then small 
#If no fl_id is given, use tl_id for assignment 
#TODO: Maybe double check for tl_id as well
#WARNING: Leaving as is, don;t use this function, use position directly (more robust) 
def split_by_id(slice_src,slice_name,mode_extra="skip",fl_id=None, tl_id=None,fl_id_coords=None,tl_id_coords=None):
    print("==Spliting by ID...")
    #Split into large/small
    large, small = split_large_small(slice_src,slice_name,mode_extra)
    #Default (True Lumen only)
    tl = large; fl = small; size_tl = "LARGE"; size_fl = ""    
    
    #One check should be enough (fl in small?), but we check both ids against both regions for discrepencies
    #Use fl_id for assignment, if not provided, use tl_id
    if small is not None: # No FL otherwise
        # #Quick test 
        # tl_ids = [1540,688,1540,688]
        # fl_ids = [1540,688,688,1540]
        # for tl_id, fl_id in zip(tl_ids,fl_ids):
        print("Check LARGE slice...")
        [tl_id_in_large,fl_id_in_large], _ = is_id_in_src(large, [tl_id,fl_id])
        print("Check Small slice...")
        [tl_id_in_small,fl_id_in_small], _ = is_id_in_src(small, [tl_id,fl_id])
        # [tl_id_in_small,fl_id_in_small], _ = is_id_in_src(large, [tl_id,fl_id])

        # backup_split_method = "Using default assignment..."
        backup_split_method = "Spliting by position..."
        if fl_id_in_large and tl_id_in_large: #Both large
            print(f"[WARNING!] fl_id ({fl_id}) and tl_id ({tl_id}) are both in the LARGE slice!!! {backup_split_method}...")
            tl, fl, size_tl, size_fl = split_by_position(slice_src,slice_name,mode_extra=mode_extra,fl_id_coords=fl_id_coords,tl_id_coords=tl_id_coords,large=large,small=small)
        elif tl_id_in_small and fl_id_in_small: #Both small 
            print(f"[WARNING!] fl_id ({fl_id}) and tl_id ({tl_id}) are both in the SMALL slice!!! {backup_split_method}...")
            tl, fl, size_tl, size_fl = split_by_position(slice_src,slice_name,mode_extra=mode_extra,fl_id_coords=fl_id_coords,tl_id_coords=tl_id_coords,large=large,small=small)
        elif fl_id_in_large and fl_id_in_small: #FL: TRUE | TRUE
            print(f"[WARNING!] fl_id ({fl_id}) is in both LARGE and SMALL slices!!! {backup_split_method}...")
            tl, fl, size_tl, size_fl = split_by_position(slice_src,slice_name,mode_extra=mode_extra,fl_id_coords=fl_id_coords,tl_id_coords=tl_id_coords,large=large,small=small)
        elif tl_id_in_large and tl_id_in_small: #TL: TRUE | TRUE
            print(f"[WARNING!] tl_id ({tl_id}) is in both LARGE and SMALL slices!!! {backup_split_method}...")
        elif tl_id is None and not fl_id_in_large and not fl_id_in_small: #FL: FALSE | FALSE
            print(f"[WARNING!] fl_id ({fl_id}) is in neither LARGE and SMALL slices!!! {backup_split_method}...")            
            tl, fl, size_tl, size_fl = split_by_position(slice_src,slice_name,mode_extra=mode_extra,fl_id_coords=fl_id_coords,tl_id_coords=tl_id_coords,large=large,small=small)
        elif fl_id_in_large or tl_id_in_small: #Switch
            tl = small; fl = large; size_tl = "SMALL"; size_fl = "LARGE"
            print(f"--fl_id ({fl_id}) is in LARGE slice OR tl_id ({tl_id}) is in SMALL slice")            
        elif fl_id_in_small or tl_id_in_large: #fl_id_in_small (default when small exists)
            tl = large; fl = small; size_tl = "LARGE"; size_fl = "SMALL"    
            print(f"--fl_id ({fl_id}) is in SMALL slice OR tl_id ({tl_id}) is in LARGE slice")            
        elif not any([fl_id_in_large,fl_id_in_small,tl_id_in_large,tl_id_in_small]): #(shouldn't be reached)
            print(f"[WARNING!] The IDs could neither be found in the LARGE nor small slices--fl_id ({fl_id}) | tl_id ({tl_id}) - {backup_split_method}...")            
            tl, fl, size_tl, size_fl = split_by_position(slice_src,slice_name,mode_extra=mode_extra,fl_id_coords=fl_id_coords,tl_id_coords=tl_id_coords,large=large,small=small)
        else: #(should never be reached)
            print(f"[ELSE!--WARNING!] --fl_id ({fl_id}) | tl_id ({tl_id})")            
            print(f"[ELSE!--WARNING!] --fl_id_in_large ({fl_id_in_large}) | fl_id_in_small ({fl_id_in_small})")            
            print(f"[ELSE!--WARNING!] --tl_id_in_large ({tl_id_in_large}) | tl_id_in_small ({tl_id_in_small})")            
            raise ValueError(f"Check if/else logic. --fl_id ({fl_id}) | tl_id ({tl_id})\n")          

    return tl, fl, size_tl, size_fl

#Split into large/small, then assign TL/FL based on position of centroid and given coordinates 
#If fl_id_coords is closer to large, then large, if it is closer to small, then small 
#If no fl_id_coords is given, use tl_id_coords for assignment 
#TODO: Maybe double check for tl_id as well
def split_by_position(slice_src,slice_name,mode_extra="skip",fl_id_coords=None,tl_id_coords=None,large=None,small=None,tolerance=2):
    print("==Spliting by position...")
    if large is None or small is None:
        #Split into large/small
        large, small = split_large_small(slice_src,slice_name,mode_extra)
    #Default (True Lumen only)
    tl = large; fl = small; size_tl = "LARGE"; size_fl = ""    
    
    if small is not None: # No FL otherwise
        #Get centroid of large and small regions 
        coords_centroid_large, _ = get_centroid(large)
        coords_centroid_small, _ = get_centroid(small)
        #Check for each centroid if the TL or the FL is the closest 
        small_closer_tl, min_distance_tl = is_small_closer(tl_id_coords,coords_centroid_small,coords_centroid_large)
        small_closer_fl, min_distance_fl = is_small_closer(fl_id_coords,coords_centroid_small,coords_centroid_large)
        
        check_tolerance(min_distance_tl,tolerance=tolerance)
        check_tolerance(min_distance_fl,tolerance=tolerance)
        # #Quick test 
        # small_closer_tls = [None,None,True,True,False,True,False]
        # small_closer_fls = [None,True,None,False,True,True,False]
        # # min_distance_tl = 200; min_distance_fl = 300
        # min_distance_tl = 300; min_distance_fl = 200
        # for small_closer_tl, small_closer_fl in zip(small_closer_tls,small_closer_fls):
        if small_closer_tl is None and small_closer_fl is None: #Both None
            print(f"[WARNING!] Most likely no coordinates were provided. small_closer_tl ({small_closer_tl}) | small_closer_fl ({small_closer_fl})")            
        elif small_closer_tl is True and small_closer_fl is True: #True, True
            print("[WARNING!]  --small region is closer to both the TL point and the FL point than the LARGE region.")
            print(f"small_closer_tl ({small_closer_tl}) | small_closer_fl ({small_closer_fl})")
            tl, fl, size_tl, size_fl = assign_using_tl_fl_dist_to_region(large,small,min_distance_tl,min_distance_fl,region="small")
        elif small_closer_tl is False and small_closer_fl is False: #False, False
            print("[WARNING!]  --LARGE region is closer to both the TL point and the FL point than the small region.")
            print(f"small_closer_tl ({small_closer_tl}) | small_closer_fl ({small_closer_fl})")
            tl, fl, size_tl, size_fl = assign_using_tl_fl_dist_to_region(large,small,min_distance_tl,min_distance_fl,region="LARGE")
        elif (small_closer_tl or small_closer_tl is None) and not small_closer_fl: #True/None, None/False - SWITCH #TL is small
            tl = small; fl = large; size_tl = "SMALL"; size_fl = "LARGE"
            print(f"--SMALL_CLOSER_TL ({small_closer_tl}, dist: {min_distance_tl}) AND small_closer_fl ({small_closer_fl}, dist: {min_distance_fl})")            
        elif not small_closer_tl and (small_closer_fl or small_closer_fl is None): #None/False, True/None, TL is LARGE
            tl = large; fl = small; size_tl = "LARGE"; size_fl = "SMALL"    
            print(f"--small_closer_tl ({small_closer_tl}, dist: {min_distance_tl}) AND SMALL_CLOSER_FL ({small_closer_fl}, dist: {min_distance_fl})")            
        else: 
            print(f"[ELSE!--WARNING!] Should not happen... Debug. small_closer_tl ({small_closer_tl}) | small_closer_fl ({small_closer_fl})")            
            print(f"[ELSE!--WARNING!] Should not happen... Debug. min_distance_tl ({min_distance_tl}) | min_distance_fl ({min_distance_fl})")  
            raise ValueError(f"Check if/else logic. small_closer_tl ({small_closer_tl}) | small_closer_fl ({small_closer_fl}) \n min_distance_tl ({min_distance_tl}) | min_distance_fl ({min_distance_fl})")          

    return tl, fl, size_tl, size_fl

#Check if the point representing the true lumen OR the false lumen is closer to the LARGE or the small slice
def is_small_closer(surface_id_coords,coords_centroid_small,coords_centroid_large):
    min_distance, min_index, distances = find_closest_point_np(surface_id_coords,[coords_centroid_small,coords_centroid_large])
    if min_index is None:
        small_closer = None
    elif min_index == 0: 
        small_closer = True; #fl_closer = False;
    elif min_index == 1:
        small_closer = False; #fl_closer = True;
    else: 
        raise ValueError(f"[ELSE!--WARNING!] Should never happen, small or large must be closer")  

    return small_closer, min_distance 

#Return warning if distance is higher than tolerance 
def check_tolerance(dist,tolerance=2): 
    if dist is not None: 
        if dist > tolerance: 
            warn_msg = f"[WARNING] Distance {dist} is higher than tolerance {tolerance}"
            print(warn_msg)
            warnings.warn(warn_msg)

#If FL is closer to the given region, assign FL as that region (otherwise assign TL to that region)
def assign_using_tl_fl_dist_to_region(large,small,min_distance_tl,min_distance_fl,region): 
    if region.lower() not in ["small","large"]:
        raise ValueError(f"Bad Region: {region}! Should be small or LARGE.") 
    print(f"Using distances to {region} to assign...")
    fl_closer_to_region = min_distance_tl>min_distance_fl
    print("--min_distance_tl: {} | min_distance_fl: {} | min_distance_tl>min_distance_fl: {}".format(min_distance_tl,min_distance_fl,fl_closer_to_region))
    if (fl_closer_to_region and region.lower() == "small") or \
             (not fl_closer_to_region and region.lower() == "large"): 
        print(f"TL as LARGE AND FL as small")
        tl = large; fl = small; size_tl = "LARGE"; size_fl = "SMALL"  
    else:
        print(f"TL as small AND FL as LARGE")
        tl = small; fl = large; size_tl = "SMALL"; size_fl = "LARGE"
    return tl, fl, size_tl, size_fl 

#Extract and integrate surface 
def process_slice(slice_src,slice_name,results_df=None,press_name='Pressure [mmHg]',vel_name='Velocity [m/s]',
        slice_size=None,ID_ARRAY="GlobalElementID",disp_name = "Displacement",disp_available=False): 
    mean_press_name = 'mean({}) {}'.format(press_name,slice_name)
    mean_vel_name = 'mean({}) {}'.format(vel_name,slice_name)
    mean_vel_mag_name = 'mag(mean({})) {}'.format(vel_name,slice_name)
    mean_flow_name = 'mag(mean({}))*Area {}'.format(vel_name,slice_name)
    mean_flow_name = 'Flow {}'.format(slice_name)
    area_name = 'Area {}'.format(slice_name)
    names_to_extract_int = [mean_press_name,mean_vel_mag_name]
    names_to_extract_calc = [mean_press_name,mean_vel_mag_name]
    #===Displacement for FSI
    mean_disp_name = 'mean({}) {}'.format(disp_name,slice_name)
    mean_disp_mag_name = 'mag(mean({})) {}'.format(disp_name,slice_name)
    mag_disp_name = 'mag({}) {}'.format(disp_name,slice_name)
    mag_disp_max_name = 'max(mag({})) {}'.format(disp_name,slice_name)
    all_names = {"mean_press_name":mean_press_name, 
        "mean_vel_name":mean_vel_name, "mean_vel_mag_name":mean_vel_mag_name, 
        "mean_disp_name":mean_disp_name , "mean_disp_mag_name":mean_disp_mag_name , 
        "mag_disp_name":mag_disp_name, "mag_disp_max_name":mag_disp_max_name }
    #=== END - Displacement for FSI
    surface_id_name = '{} {}'.format(ID_ARRAY,slice_name)
    surface_id_pos_name = '{} Coordinates {}'.format(ID_ARRAY,slice_name)
    slice_size_name = 'Size {}'.format(slice_name)

    extractSurface, empty_slice = extract_surf_and_check_connectivity(slice_name,slice_src)
    if empty_slice:
        return results_df, None, None, None
    

    #Integrate surface and get the GlobalElementID of the centroid if available 
    integrateVariables, area, _, surface_id, surface_id_coords = find_closest_cell_to_centroid(slice_name,extractSurface)
    
    names_dic = {"mean_press":mean_press_name,"mean_vel_mag":mean_vel_mag_name,"mean_flow":mean_flow_name}

    #Check is displacement is available (FSI)
    if not disp_available: 
        array_names, _ = get_array_names(sm.Fetch(integrateVariables),data_type="PointData",print_all=False)
        if disp_name in array_names: 
            disp_available = True
            print(f"{disp_name} found in arrays, will compute related stats")

    if disp_available: 
        #Get displacement, flow, and pressure
        Integral_flow, first_point = compute_press_flow_disp(slice_name,integrateVariables,extractSurface, all_names,
                press_name=press_name,vel_name=vel_name, disp_name=disp_name,integral=True,calculator=True)
        names_to_extract_int += [mean_disp_mag_name]
        names_to_extract_calc += [mean_disp_mag_name,mag_disp_max_name]
        names_dic.update({"mean_disp_mag":mean_disp_mag_name,"mag_disp_max":mag_disp_max_name})
    else: 
        #Get flow and pressure
        Integral_flow, first_point = compute_press_flow(slice_name,integrateVariables,extractSurface, mean_press_name, mean_vel_name, mean_vel_mag_name,
                press_name=press_name,vel_name=vel_name, integral=True,calculator=True)  

    
    names_dic_reversed = {v: k for k, v in names_dic.items()} #To avoid changing code after filter 
    # data_dic_o = {"mean_press":None,"mean_vel_mag":None,"mean_flow":None}
    data_dic_int = {}
    data_dic_calc = {}
    data_dic_mix = {}
    #TODO Reformat Group by timesteps takes way too long, might even be better to do it for each
    #time step separately, write multiple files then group manually
    #or store values in DF, then group 
    #Integral
    # data_dic_int["mean_press"], data_dic_int["mean_vel_mag"], data_dic_int["mean_flow"], group_dsa_int = \
    #         _filter_to_results(Integral_flow,"integral",slice_name,area,mean_press_name,mean_vel_mag_name)
    data_dic_int, group_dsa_int = filter_to_results(Integral_flow,"integral",slice_name,area,
        mean_vel_mag_name,names_to_extract_int)
    data_dic_int = {names_dic_reversed.get(k, k): v for k, v in data_dic_int.items()} #To avoid changing rest of code 
    #Calc
    # data_dic_calc["mean_press"], data_dic_calc["mean_vel_mag"], data_dic_calc["mean_flow"], group_dsa_calc = \
    #         _filter_to_results(first_point,"calc",slice_name,area,mean_press_name,mean_vel_mag_name)
    data_dic_calc, group_dsa_calc = filter_to_results(first_point,"calc",slice_name,area,
        mean_vel_mag_name,names_to_extract_calc)
    data_dic_calc = {names_dic_reversed.get(k, k): v for k, v in data_dic_calc.items()} #To avoid changing rest of code 

    #Keep Flow from Integral and Pressure from Calculator
    data_dic_mix["mean_press"] = data_dic_calc["mean_press"]
    data_dic_mix["mean_vel_mag"] = data_dic_int["mean_vel_mag"]
    data_dic_mix["mean_flow"] = data_dic_int["mean_flow"]

    #More info added for FSI
    if slice_size is not None: 
        if disp_available:
            #Putting both for now to check correctness
            data_dic_mix["mean_disp_mag_int"] = data_dic_int["mean_disp_mag"]
            names_dic["mean_disp_mag_int"] = names_dic["mean_disp_mag"] + "_int"
            data_dic_mix["mean_disp_mag"] = data_dic_calc["mean_disp_mag"]
            data_dic_mix["mag_disp_max"] = data_dic_calc["mag_disp_max"]
            #data_dic_mix["mag_disp_max"] = data_dic_calc["mag_disp_max"]
            
        {surface_id_name: [surface_id]*len(data_dic_mix["mean_press"]),
        surface_id_pos_name: [surface_id_coords]*len(data_dic_mix["mean_press"])}
        # max_disp = 
        # avg_disp = 
        df_preprend_dic = {surface_id_name: [surface_id]*len(data_dic_mix["mean_press"]),
        surface_id_pos_name: [surface_id_coords]*len(data_dic_mix["mean_press"]),
        slice_size_name: [slice_size]*len(data_dic_mix["mean_press"])}
    else: 
        df_preprend_dic = {surface_id_name: [surface_id]*len(data_dic_mix["mean_press"])}

    results_df = results_to_df(data_dic_mix, names_dic,area,area_name,results_df,df_preprend_dic)
    # results_df = results_to_df(data_dic_calc, names_dic,area,area_name,results_df)
    # results_df = results_to_df(data_dic_int, names_dic,area,area_name,results_df)
    if "DTA2" in slice_name:
        # import pdb; pdb.set_trace()
        print(data_dic_mix["mean_press"])
    return results_df, extractSurface, integrateVariables, [group_dsa_int,group_dsa_calc]

#Extract surface and check connectivity 
def extract_surf_and_check_connectivity(slice_name,slice_src):
    # create a new 'Extract Surface'
    extractSurface = ExtractSurface(registrationName='ExtractSurface_' + slice_name , Input=slice_src)
    if extractSurface.PointData.GetNumberOfArrays() == 0:
        print("--[WARNING!-Connectivity] NO region found")
        print("    Returning given slice as LARGE")
        empty_slice = True 
    else: 
        empty_slice = False
    return extractSurface, empty_slice 

#Compute centroid and find closest cell to centroid 
def find_closest_cell_to_centroid(slice_name,surface_src,ID_ARRAY="GlobalElementID"): 
    integrateVariables = IntegrateVariables(registrationName='IntegrateVariables_' + slice_name, Input=surface_src)
    area = sm.Fetch(integrateVariables).GetCellData().GetArray("Area").GetValue(0)
    #Centroid specifically 
    coords_centroid = vtk_to_numpy(sm.Fetch(integrateVariables).GetPoints().GetData())[0]

    surface_id_coords, closest_cell_id, _, dist2 = find_closest_cell_to_point(surface_src,coords_centroid.tolist(),print_all=False)
    surface_ids_array = sm.Fetch(surface_src).GetCellData().GetArray(ID_ARRAY)
    if surface_ids_array is not None: 
        surface_id = vtk_to_numpy(surface_ids_array)[closest_cell_id]
    else:
        surface_id = -1 

    #First point from surface 
    # surface_id = sm.Fetch(surface_src).GetCellData().GetArray("surface_src").GetValue(0)
    #TODO: Is min good? ==> min/first ID doesn't work well
    # surface_id = np.min(vtk_to_numpy(sm.Fetch(surface_src).GetCellData().GetArray("surface_src")))
    print("    {} Area: {}".format(slice_name,area)) 
    print("    {} Surface ID: {}".format(slice_name,surface_id)) 
    print("    {} Surface ID Coords: {}".format(slice_name,surface_id_coords)) 
    print(f"{slice_name} Distance from centroid to closest cell: {dist2}")

    return integrateVariables, area, coords_centroid, surface_id, surface_id_coords

#Go from Large to Small slices to FL/TL based on the helper file false_lumen_id_[...].tsv
#Accept results_tsv as file name or dataframe 
def convert_slice_names(results_tsv,fl_id_tsv):
    if isinstance(results_tsv,str):
        results_df = pd.read_csv(results_tsv,sep="\t")
        save_tsv = True
    else: 
        results_df = results_tsv
        save_tsv = False 

    fl_id = pd.read_csv(fl_id_tsv,sep="\t")
    large_fl = fl_id.loc[fl_id["Is FL"].isnull() | (fl_id["Is FL"] == "NO")]["Slice Name"].tolist()

    # fl_id = fl_id.set_index("Slice Name")
    new_names = {}
    tl_fl_mapping = {} 

    #Go throught columns to find new names
    for col in results_df.columns:
        slice_part = col.split(" ")[-1]
        
        if slice_part in large_fl:
            new_names[col] = col.replace("SMALL","TL")
            #Replace corresponding LARGE as FL
            new_names[col.replace("SMALL","LARGE")] = col.replace("SMALL","FL")
        else: 
            #Default, replace large with TL and small with FL
            if "LARGE" in slice_part:
                new_names[col] = col.replace("LARGE","TL")
            #We avoid the corresponding SMALL (replaced above)
            elif "SMALL" in slice_part and slice_part.replace("SMALL","LARGE") not in large_fl:
                # print(new_names[col] )
                new_names[col] = col.replace("SMALL","FL")
        # if col not in "Timestep":
        #     print(new_names[col])

    results_df = results_df.rename(columns=new_names)
    
    #Extract slice_names mapping (some are repeated, but using dic)
    #Can be used to arrange slices and only plot fl or tl correctly
    for key, val in new_names.items():
        if key.lower() != "timestep" or not isinstance(results_tsv,str):
            tl_fl_mapping[key.split(" ")[-1]] = val.split(" ")[-1]
    print(tl_fl_mapping)

    if save_tsv:
        results_tsv_tl_fl = add_suffix(results_tsv,"_tl_fl")
        print("Saving {}...".format(results_tsv_tl_fl))
        results_df.to_csv(results_tsv_tl_fl,sep="\t",index=False)

        mapping_path = results_tsv.replace("results","mapping").replace(".tsv",".txt")
        print("Saving {}...".format(mapping_path))
        with open(mapping_path, 'w') as f:  
            for key, value in tl_fl_mapping.items():  
                f.write('%s: %s\n' % (key, value))

    return results_df, tl_fl_mapping 


#Wrapper - Given a list of the current slices, return a mapping to the 
# TL/FL names
def convert_slice_names_list(slice_names,fl_id_tsv):
    results_tsv = pd.DataFrame(columns=slice_names)
    results_df, tl_fl_mapping  = convert_slice_names(results_tsv,fl_id_tsv)
     
    return results_df.columns, tl_fl_mapping

#Convert all results in drectory 
def convert_slice_names_dir(result_folder):
    #Get results not converted yet 
    all_results =glob.glob(os.path.join(result_folder,"results*.tsv"))
    print(all_results)
    all_results = [x for x in all_results if "converted" not in x]
    
    for results_tsv in all_results:
        fl_id_tsv = results_tsv.replace("results","false_lumen_id")
        if os.path.exists(fl_id_tsv):
            print("Converting {}".format(fl_id_tsv))
            convert_slice_names(results_tsv,fl_id_tsv)
        else:
            print("[FAIL] Cannot find {}".format(fl_id_tsv))

#endregion         END - SLICES FUNCTIONS
#==========================================

#==========================================
#region          fl_id_tsv FUNCTIONS
#==========================================

#Get element ID of small slices that belong to the FL and those that belong to the TL
def get_fl_ids(fl_id_tsv,rename_small_by_tl_fl=True,ID_ARRAY="GlobalElementID"):
    small_tl, small_fl = [None, None]
    if os.path.exists(fl_id_tsv):
        fl_id = pd.read_csv(fl_id_tsv,sep="\t")
        coords_col = f"{ID_ARRAY} Coordinates"
        if ID_ARRAY not in fl_id.columns or coords_col not in fl_id.columns:
            print(f"No {ID_ARRAY} in columns of {fl_id_tsv}")
        elif all(fl_id[ID_ARRAY] == -1): #Maybe any? might not be interested in that slice
            print(f"All IDs in {fl_id_tsv} are -1.")
        else:
            small_tl = fl_id.loc[fl_id["Is FL"].isnull() | (fl_id["Is FL"] == "NO")].set_index("Slice Name")[[ID_ARRAY, coords_col]].to_dict()
            print(f"[LARGE FL/small TL] (NO change needed): {small_tl}")
            small_fl = fl_id.loc[fl_id["Is FL"] == "YES"].set_index("Slice Name")[[ID_ARRAY, coords_col]].to_dict()
            print(f"[small FL/LARGE TL] (change needed): {small_fl}")

            if rename_small_by_tl_fl:
                small_tl = {key: dic_str_to_list(replace_key(val,'SMALL','TL')) for key, val in small_tl.items()} 
                small_fl = {key: dic_str_to_list(replace_key(val,'SMALL','FL')) for key, val in small_fl.items()} 
    else: 
        print(f"File {fl_id_tsv} does not exist... Returning None, None")
    return small_tl, small_fl

#Replace string in key
def replace_key(my_dic,old_str,new_str): 
    return {key.replace(old_str, new_str): value for key, value in my_dic.items()}

#Turn values that are strings to list/float, etc...
def dic_str_to_list(my_dic): 
    return {key: ast.literal_eval(value) if isinstance(value, str) else value
                for key, value in my_dic.items()}
     
#Write fl_id_tsv
#User: Use the file that this function writes to identify the False Lumen manually
def write_fl_id_tsv(fl_id_tsv,surface_id_fl,surface_id_coords_fl,slice_names_fl,overwrite=False,ID_ARRAY="GlobalElementID"):
    if not os.path.exists(fl_id_tsv) or overwrite:
        surface_id_coords_fl = [str(s) for s in surface_id_coords_fl]
        print("Saving {}...".format(fl_id_tsv))
        # fl_id = pd.DataFrame({"Slice Name": slice_names_fl, "Is FL":["YES"]*len(slice_names_fl)})
        fl_id = pd.DataFrame({ID_ARRAY:surface_id_fl,f"{ID_ARRAY} Coordinates": surface_id_coords_fl,"Slice Name": slice_names_fl, "Is FL":["YES"]*len(slice_names_fl)})
        fl_id.to_csv(fl_id_tsv,sep="\t",index=False)

#Determine suffix of fl_id_tsv and return the updated path 
def get_fl_id_tsv_full_name(fl_id_tsv,results_tsv): 
    if fl_id_tsv != "":
        if "fixed" in results_tsv:
            fl_id_tsv = add_suffix(fl_id_tsv,"_fixed")
        else: 
            fl_id_tsv = add_suffix(fl_id_tsv,"_{}".format(len(slice_names)))
    return fl_id_tsv 

#Read fl_id_tsv and convert from large/small to tl/fl
def process_fl_id_tsv(fl_id_tsv,process_slice_bool, slices, slice_names, results_tsv, results_df):
    if os.path.exists(fl_id_tsv):
        print("Using {}...".format(fl_id_tsv))
        if not process_slice_bool: 
            all_slices_names, tl_fl_mapping = convert_slice_names_list(slices.keys(),fl_id_tsv)
        else: #Comment out the else case if you want to keep LARGE/SMALL
            results_df, tl_fl_mapping = convert_slice_names(results_tsv,fl_id_tsv)
        slices = map_fields(slices,tl_fl_mapping)

    return slices, slice_names, results_df 

#endregion         END - fl_id_tsv FUNCTIONS
#==========================================

#==========================================
#region          DISPLAY FUNCTIONS
#==========================================

#Show all the slices
#Make an empty view before  
#filter_slices_shown: Only show slices that have filter_slices_shown in their name
def show_slices(slices,slice_names,input_src=None,input_src_opacity=0.05,bar_range=[0,1],
    view_time=3,myview=None,mode="visible",input_src_color=None,filter_slices_shown="",
    tl_opacity=None,preset="rainbow", color_by_array='Velocity [m/s]',title='Velocity Magnitude [m/s]',
    other_regions_to_show=None): 
    if myview is not None:
        SetActiveView(myview) #Eventually make format_colorbar accept a view
    HideAll()
    myview = GetActiveView()
    myview.ViewTime = view_time
    print("View Time: {}".format(myview.ViewTime))
    if tl_opacity is not None: 
        separate_lumen = True #Show both lumen separately
    else: 
        separate_lumen = False #Show both lumen separately

    #Only show slices that have filter_slices_shown in their name
    if filter_slices_shown != "":
        separate_lumen = True
        print(f"Filtering shown slices by {filter_slices_shown}...") 
        slices_filtered = {key: value for key, value in slices.items() if filter_slices_shown in key}
        if len(slices_filtered) == 0:
            print(f"[WARNING!] Filtering by ...{filter_slices_shown}... might have removed all slices.")
            print(f"Available names: {slices.keys()}")
        else:
            print(f"Slices left: {slices_filtered.keys()}")
        slices = slices_filtered
   
    use_tl_fl = False 
    for key in slices.keys():
        if "TL_" in key or "FL_" in key:
            use_tl_fl = True 
            break
    disp = None
    for slice_name in slice_names.keys(): 
        if use_tl_fl:
            slice_name_tl = "TL_" + slice_name
            slice_name_fl = "FL_" + slice_name
        else:
            slice_name_tl = "LARGE_" + slice_name
            slice_name_fl = "SMALL_" + slice_name

        #mode fixed, bar range not used
        if separate_lumen:
            if slices[slice_name_tl] != None and slices[slice_name_tl] == slices[slice_name_fl]: 
                warm_msg = "[!!!WARNING!!!] show_slices() - Both slices are the same. Debugging needed..."
                print(warm_msg)
                warnings.warn(warm_msg)
            if slice_name_tl in slices and slices[slice_name_tl] is not None:
                print(f"{slice_name_tl}: {get_area(slices[slice_name_tl])}")
                disp,disp_LUT = format_colorbar(slices[slice_name_tl],bar_range=bar_range,mode="ignore",
                        preset=preset,color_by_array=color_by_array,title=title)
                if tl_opacity is not None:
                    disp.Opacity = tl_opacity
                Render()
            if slice_name_fl in slices and slices[slice_name_fl] is not None:
                print(f"{slice_name_fl}: {get_area(slices[slice_name_fl])}")
                disp,disp_LUT = format_colorbar(slices[slice_name_fl],bar_range=bar_range,mode="ignore",
                        preset=preset,color_by_array=color_by_array,title=title)
            
        else: #Whole slice 
            if slice_name in slices and slices[slice_name] is not None:
                # disp = format_colorbar(slices[slice_name],bar_range=bar_range) 
                disp,disp_LUT = format_colorbar(slices[slice_name],bar_range=bar_range,mode="ignore",
                        preset=preset,color_by_array=color_by_array,title=title)     
        # save_screen("centerline.png",image_res=[800,500])
        # import pdb; pdb.set_trace()

    #No need to rescale every time in for loop when not debugging 
    # Render()
    if disp is not None: 
        set_scale_range(disp,mode,bar_range,disp_LUT)
        Render()
        if input_src is not None:
            disp, _ = format_colorbar(input_src,opacity=input_src_opacity,mode="ignore",separate=True,show_bar=False,bar_range=bar_range,input_src_color=input_src_color,
                        preset=preset,color_by_array=color_by_array,title=title)
            # set_solid_col(disp,opacity=0.2,color="grey")
            Render()
            # save_screen("centerline.png",image_res=[800,500])
            show_other_regions(other_regions_to_show,bar_range=bar_range,mode="ignore",separate=True,preset=preset,color_by_array=color_by_array,title=title)

#Show partial regions as provided, mode should be ignore to avoid updating the color bar 
#Separate does not need to be true if no change is made to the color bar after 
#Main goal: show regions with different opacity/color(solid) 
def show_other_regions(other_regions_to_show,bar_range=[0,1],mode="ignore",separate=True,preset="rainbow", color_by_array='Velocity [m/s]',title='Velocity Magnitude [m/s]'):
    if other_regions_to_show is not None:
        for key,region_to_show in other_regions_to_show.items():
            print(f"Adding to visualization: {key} | {region_to_show.__str__()}")
            disp, _ = format_colorbar(region_to_show.src,opacity=region_to_show.opacity,mode=mode,separate=separate,show_bar=False,bar_range=bar_range,input_src_color=region_to_show.color,
                        preset=preset,color_by_array=color_by_array,title=title)
            # set_solid_col(disp,opacity=0.2,color="grey")
            Render()
    

#Save all time views
def save_all_times(save_path,slices,slice_names,input_src=None,input_src_opacity=0.05,mode="visible",bar_range=[0,1],image_res=[800,500],
    all_times=None,myview=None,input_src_color=None,filter_slices_shown="",tl_opacity=None,preset="rainbow", color_by_array='Velocity [m/s]',title='Velocity Magnitude [m/s]',
    other_regions_to_show=None): 
    if all_times is None: 
        all_times = input_src.TimestepValues
    save_path = make_parent_dir(save_path,"slices_all_times")
    print("TimestepValues: {}".format(all_times))
    for view_time in all_times: 
        save_path_i = add_suffix(save_path,"_{}".format(view_time).replace(".","_"))
        #TODO REMOVE - overwrite eventually
        if glob.glob("{}*".format(save_path_i.rsplit(".",1)[0])):
            break
        print("Saving {}".format(save_path_i))
        show_slices(slices,slice_names,input_src=input_src,input_src_opacity=input_src_opacity,mode=mode,bar_range=bar_range,
            view_time=view_time,myview=myview,input_src_color=input_src_color,filter_slices_shown=filter_slices_shown,tl_opacity=tl_opacity, preset=preset,color_by_array=color_by_array,title=title,
                other_regions_to_show=other_regions_to_show)
        save_screen(save_path_i,image_res=image_res,myview=myview)

#Save all time views
def save_all_times2(save_path,input_src=None,image_res=[800,500],all_times=None): 
    if all_times is None: 
        all_times = input_src.TimestepValues
    save_path = make_parent_dir(save_path,"all_times")
    print("TimestepValues: {}".format(all_times))
    for view_time in all_times: 
        GetActiveView().ViewTime = view_time
        save_path_i = add_suffix(save_path,"_{}".format(view_time).replace(".","_"))
        print("Saving {}".format(save_path_i))
        save_screen(save_path_i,image_res=image_res)

def find_viewed_vtu(src):
    vtu_cur = "Not Found"
    index_time = -1
    tstep = GetActiveView().ViewTime
    for i,time_i in enumerate(src.TimestepValues):
        if tstep == time_i:
            index_time = i
    if index_time != -1:
        vtu_cur = src.FileName[index_time]

    return vtu_cur

#endregion         END - DISPLAY FUNCTIONS
#==========================================

#==========================================
#region          PRESSURE/FLOW FUNCTIONS
#==========================================

#Takes source, extract pressure and flow info 
#names_to_extract: Pass all names that need to be extracted as a list
#mean_vel_mag_name still needs to be specied separately so mean_flow can be added 
# def filter_to_results(Input,src_type,slice_name,area,mean_press_name,mean_vel_mag_name,names_to_extract):
def filter_to_results(Input,src_type,slice_name,area,mean_vel_mag_name,names_to_extract):
    #Group Results by Time Step 
    groupTimeSteps = _GroupTimeSteps(registrationName='GroupTimeSteps_{}_{}'.format(src_type,slice_name), Input=Input)

    #Use numpy interface _.PointData[...][0] is a VTKCompositeDataArray
    # Access array with .Arrays
    group_dsa = dsa.WrapDataObject(sm.Fetch(groupTimeSteps))
    all_arrays = {}
    for name in names_to_extract: 
        all_arrays[name] = group_dsa.PointData[name][0]
    if src_type.lower() == "integral":
        #Need to divide all values by the area, except the flow 
        for name in all_arrays.keys(): 
            if name != mean_vel_mag_name:
                all_arrays[name] = (all_arrays[name]/area)
            else: 
                mean_flow = all_arrays[name]
    else:
        # mean_flow = mean_vel_mag*area name != mean_vel_mag_name is never true
        mean_flow = (all_arrays[mean_vel_mag_name]*area)
    #Letting this throw an error if 
    all_arrays["mean_flow"] = mean_flow

    # Access array with .Arrays
    for name in all_arrays.keys(): 
        all_arrays[name] = all_arrays[name].Arrays


    # print(group_dsa.PointData[mean_vel_mag_name][0])
    # print(mean_press)
    # print(mean_flow)

    return all_arrays, group_dsa
    # return mean_press.Arrays, mean_vel_mag.Arrays, mean_flow.Arrays, group_dsa

#Takes source, extract pressure and flow info 
def _filter_to_results(Input,src_type,slice_name,area,mean_press_name,mean_vel_mag_name):
    #Group Results by Time Step 
    groupTimeSteps = _GroupTimeSteps(registrationName='GroupTimeSteps_{}_{}'.format(src_type,slice_name), Input=Input)

    #Use numpy interface _.PointData[...][0] is a VTKCompositeDataArray
    # Access array with .Arrays
    group_dsa = dsa.WrapDataObject(sm.Fetch(groupTimeSteps))
    mean_press = group_dsa.PointData[mean_press_name][0]
    mean_vel_mag= group_dsa.PointData[mean_vel_mag_name][0]
    if src_type.lower() == "integral":
        mean_flow = mean_vel_mag
        mean_press = mean_press/area
    else:
        mean_flow = mean_vel_mag*area 

    # print(group_dsa.PointData[mean_vel_mag_name][0])
    # print(mean_press)
    # print(mean_flow)

    return mean_press.Arrays, mean_vel_mag.Arrays, mean_flow.Arrays, group_dsa

#Turn pressure and flow into a df
def results_to_df(data,names,area,area_name,results_df=None,df_preprend_dic={}):
    results_df_i_dic_default = {area_name: [area]*len(data["mean_press"])}
    results_df_i_dic_data = {names[key]: data[key] for key in names.keys()}
    results_df_i_dic_default.update(results_df_i_dic_data)
    # names["mean_press"]: data["mean_press"],
    # names["mean_vel_mag"]: data["mean_vel_mag"],
    # names["mean_flow"]: data["mean_flow"]
    if len(df_preprend_dic) > 0:
        df_preprend_dic.update(results_df_i_dic_default)
        results_df_i = pd.DataFrame(df_preprend_dic)
    else: 
        results_df_i = pd.DataFrame(results_df_i_dic_default)

    if results_df is None:
        results_df = results_df_i 
    else:
        results_df = pd.concat([results_df,results_df_i],axis=1)

    return results_df 

#Compute pressure, mag(mean(velocity)), max(mag(displacement)), and  mag(mean(displacement)) 
#using integral and/or calculators
def compute_press_flow_disp(slice_name,integrateVariables,extractSurface, all_names,
            press_name='Pressure [mmHg]',vel_name='Velocity [m/s]', disp_name="Displacement",
            integral=True,calculator=True):
    #Moved to dic since the list is getting long 
    # all_names = {"mean_press_name":, "mean_vel_name":, "mean_vel_mag_name":,
    #     "mean_disp_name": ,"mean_disp_mag_name": , "mag_disp_name": , "mag_disp_max_name": }
    #========INTEGRAL
    if integral:
        #Displacement from python 
        Integral_disp = compute_disp_integral(slice_name,integrateVariables, all_names["mean_disp_mag_name"], disp_name=disp_name)

        #Pressure/Flow from integral
        Integral_flow = compute_press_flow_integral(slice_name,Integral_disp, all_names["mean_press_name"], all_names["mean_vel_mag_name"],
        press_name=press_name,vel_name=vel_name)

        # print_debug_info(Integral_flow,"Integral",mean_press_name, mean_vel_mag_name,
        #     press_name='Pressure [mmHg]',vel_name='Velocity [m/s]')
    else: 
        Integral_flow = None
        
    #========CALCULATOR
    if calculator:
        #Displacement from python 
        Calculator_disp = compute_disp_python(slice_name,extractSurface, all_names["mag_disp_name"], 
            all_names["mag_disp_max_name"], all_names["mean_disp_name"], all_names["mean_disp_mag_name"],
            disp_name=disp_name)
        #Pressure/Flow from python
        Calculator_flow = compute_press_flow_python(slice_name,Calculator_disp, all_names["mean_press_name"], all_names["mean_vel_name"], all_names["mean_vel_mag_name"],
            press_name=press_name,vel_name=vel_name)
        

        #Extract the first point (averages are the same across vector)
        _QuerySelect(QueryString="id == 0", FieldType="POINT", Source=Calculator_flow) #Select Centerline
        first_point = ExtractSelection()

        # print_debug_info(first_point,"Calculator",mean_press_name, mean_vel_mag_name,
        #     press_name='Pressure [mmHg]',vel_name='Velocity [m/s]')
    else: 
        Calculator_flow = None

    return Integral_flow, first_point

#Compute pressure and mag(mean(velocity)) using integral and/or calculators
def compute_press_flow(slice_name,integrateVariables,extractSurface, mean_press_name, mean_vel_name, mean_vel_mag_name,
            press_name='Pressure [mmHg]',vel_name='Velocity [m/s]', integral=True,calculator=True):

    #========INTEGRAL
    if integral:
        #Pressure/Flow from integral
        Integral_flow = compute_press_flow_integral(slice_name,integrateVariables, mean_press_name, mean_vel_mag_name,
        press_name=press_name,vel_name=vel_name)

        # print_debug_info(Integral_flow,"Integral",mean_press_name, mean_vel_mag_name,
        #     press_name='Pressure [mmHg]',vel_name='Velocity [m/s]')
    else: 
        Integral_flow = None
        
    #========CALCULATOR
    if calculator:
    #Pressure/Flow from python
        Calculator_flow = compute_press_flow_python(slice_name,extractSurface, mean_press_name, mean_vel_name, mean_vel_mag_name,
            press_name=press_name,vel_name=vel_name)

        #Extract the first point (averages are the same across vector)
        _QuerySelect(QueryString="id == 0", FieldType="POINT", Source=Calculator_flow) #Select Centerline
        first_point = ExtractSelection()

        # print_debug_info(first_point,"Calculator",mean_press_name, mean_vel_mag_name,
        #     press_name='Pressure [mmHg]',vel_name='Velocity [m/s]')
    else: 
        Calculator_flow = None

    return Integral_flow, first_point


def print_debug_info(Input,type_src,mean_press_name, mean_vel_mag_name,
            press_name='Pressure [mmHg]',vel_name='Velocity [m/s]'):
    print("Using {}:".format(type_src))
    # print(sm.Fetch(Input).GetPointData().GetArray(vel_name))
    print(vtk_to_numpy(sm.Fetch(Input).GetPointData().GetArray(mean_press_name)))
    print(vtk_to_numpy(sm.Fetch(Input).GetPointData().GetArray(vel_name)))
    print(vtk_to_numpy(sm.Fetch(Input).GetPointData().GetArray(mean_vel_mag_name)))
    # print(sm.Fetch(Input.GetPointData().GetArray(vel_name).GetMaxNorm()))
    print(vtk_to_numpy(sm.Fetch(Input).GetPointData().GetArray("velocity")))

#endregion         END - PRESSURE/FLOW
#==========================================