#Utility functions: centerline/loading/saving/printing 

#Depends on paraview calc 

import paraview 
from paraview.simple import *
from paraview.selection import _createSelection, _collectSelectionPorts
from vtk.util.numpy_support import vtk_to_numpy
from vtk.numpy_interface import dataset_adapter as dsa

import paraview.servermanager as sm 

from paraview_calc import *
from paraview_utils import save_screen, add_suffix, _QuerySelect, _GroupTimeSteps, switch_var, make_parent_dir, map_fields
from paraview_display import *  

import pandas as pd 
import numpy as np 
import os 
import glob 

#==========================================
#region          HIGH-LEVEL FUNCTIONS
#==========================================

#Load slices df, then make slices
#(see load_slice_info for info on slices_df)
def load_and_make_slices(input_src,slices_tsv,results_tsv="results.tsv",box_size=70,z_normal=True,
        show=True,force_normal=150,fl_id_tsv="",process_slice_bool=True,append_res=False): 
    slices_df = load_slices_info(slices_tsv)
    slices, slice_names, results_df = make_slices(input_src,slices_df,results_tsv=results_tsv,
        box_size=box_size,z_normal=z_normal,show=show,force_normal=force_normal,fl_id_tsv=fl_id_tsv,
        process_slice_bool=process_slice_bool,append_res=append_res)

    return slices, slice_names, results_df

#Cut along the centerline with a given stride  
def slice_centerline(input_src,centerline,box_size=70,stride=100,slices_tsv="",results_tsv="results.tsv",
        z_normal=True,show=True,force_normal=150,erase_range=[0,0],fl_id_tsv="",process_slice_bool=True,
        append_res=False):
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
            append_res=append_res)
    
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
def make_slices(input_src,slices_df,box_size=70,results_df=None,results_tsv="results.tsv",z_normal=True,
        show=True,force_normal=150,erase_range=[0,0],fl_id_tsv="",process_slice_bool=True,append_res=False):
    #slices = {"ATA":None,"DTA":None,"DTA2":None,"IAA":None}
    slice_names = {name[0]: None for name in slices_df.index.unique().values}
    slices = {}
    keys_to_rm = [] 
    slice_names_fl = [] 
    # count = 0
    # disp = Show(input_src)
    # disp.Opacity = 0.05
    # Render()
    
    # import pdb; pdb.set_trace()
    print("-----------------------------------------------")
    print(results_df)
    for slice_name in slice_names.keys(): 
        # if "IAA" not in slice_name and "DTA2" not in slice_name:
        #     continue
        # if "DTA2" not in slice_name:
        #     continue
        # slice_name_tl = "TL_" + slice_name
        # slice_name_fl = "FL_" + slice_name
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

        slices[slice_name_tl], slices[slice_name_fl] = split_large_small(slices[slice_name],slice_name)
        
        if process_slice_bool and slices[slice_name_tl] is not None:
            results_df, extractSurface , _ , _ = process_slice(slices[slice_name_tl],slice_name_tl,results_df)

        #If no False Lumen, set empty column
        if slices[slice_name_fl] is not None:
            slice_names_fl += [slice_name_fl]
            if process_slice_bool:
                results_df, _ , _ , _ = process_slice(slices[slice_name_fl],slice_name_fl,results_df)
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
    #Use this file to identify the False Lumen
    #Will be read later to convert LARGE/SMALL Areas to TL/FL
    #Does not overwrite file (remove it manually)
    # print(fl_id_tsv)
    # print(results_tsv)
    if fl_id_tsv != "":
        if "fixed" in results_tsv:
            fl_id_tsv = add_suffix(fl_id_tsv,"_fixed")
        else: 
            fl_id_tsv = add_suffix(fl_id_tsv,"_{}".format(len(slice_names)))
        if not os.path.exists(fl_id_tsv):
            print("Saving {}...".format(fl_id_tsv))
            fl_id = pd.DataFrame({"Slice Name": slice_names_fl, "Is FL":["YES"]*len(slice_names_fl)})
            fl_id.to_csv(fl_id_tsv,sep="\t",index=False)

        else:
            print("Using {}...".format(fl_id_tsv))
            if not process_slice_bool: 
                all_slices_names, tl_fl_mapping = convert_slice_names_list(slices.keys(),fl_id_tsv)
            else: #Comment out the else case if you want to keep LARGE/SMALL
                results_df, tl_fl_mapping = convert_slice_names(results_tsv,fl_id_tsv)
            slices = map_fields(slices,tl_fl_mapping)

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

        #TODO: Making this assumption for now
        if area_large < area_small: #Switch
            area_large, area_small = switch_var(area_large,area_small)
            large, small = switch_var(large,small)

        print("--Area LARGE: {} | Area SMALL: {} | LARGE>SMALL: {}".format(area_large,area_small,area_large>area_small))
        
    
    # print(sm.Fetch(large))
    return large, small

#Extract and integrate surface 
def process_slice(slice_src,slice_name,results_df=None,press_name='Pressure [mmHg]',vel_name='Velocity [m/s]'): 
    mean_press_name = 'mean({}) {}'.format(press_name,slice_name)
    mean_vel_name = 'mean({}) {}'.format(vel_name,slice_name)
    mean_vel_mag_name = 'mag(mean({})) {}'.format(vel_name,slice_name)
    mean_flow_name = 'mag(mean({}))*Area {}'.format(vel_name,slice_name)
    mean_flow_name = 'Flow {}'.format(slice_name)
    area_name = 'Area {}'.format(slice_name)

    # create a new 'Extract Surface'
    extractSurface = ExtractSurface(registrationName='ExtractSurface_' + slice_name , Input=slice_src)
    if extractSurface.PointData.GetNumberOfArrays() == 0:
        print("--[WARNING!-Connectivity] NO region found")
        print("    Returning given slice as LARGE")
        return results_df, None, None, None

    # Get Area
    integrateVariables = IntegrateVariables(registrationName='IntegrateVariables_' + slice_name, Input=extractSurface)
    area = sm.Fetch(integrateVariables).GetCellData().GetArray("Area").GetValue(0)
    print("    {} Area: {}".format(slice_name,area)) 

    #Get flow and pressure
    Integral_flow, first_point = compute_press_flow(slice_name,integrateVariables,extractSurface, mean_press_name, mean_vel_name, mean_vel_mag_name,
            press_name=press_name,vel_name=vel_name, integral=True,calculator=True)
    
    
    names_dic = {"mean_press":mean_press_name,"mean_vel_mag":mean_vel_mag_name,"mean_flow":mean_flow_name}
    
    # data_dic_o = {"mean_press":None,"mean_vel_mag":None,"mean_flow":None}
    data_dic_int = {}
    data_dic_calc = {}
    data_dic_mix = {}
    #TODO Reformat Group by timesteps takes way too long, might even be better to do it for each
    #time step separately, write multiple files then group manually
    #or store values in DF, then group 
    #Integral
    data_dic_int["mean_press"], data_dic_int["mean_vel_mag"], data_dic_int["mean_flow"], group_dsa_int = \
            filter_to_results(Integral_flow,"integral",slice_name,area,mean_press_name,mean_vel_mag_name)
    #Calc
    data_dic_calc["mean_press"], data_dic_calc["mean_vel_mag"], data_dic_calc["mean_flow"], group_dsa_calc = \
            filter_to_results(first_point,"calc",slice_name,area,mean_press_name,mean_vel_mag_name)
    
    #Keep Flow from Integral and Pressure from Calculator
    data_dic_mix["mean_press"] = data_dic_calc["mean_press"]
    data_dic_mix["mean_vel_mag"] = data_dic_int["mean_vel_mag"]
    data_dic_mix["mean_flow"] = data_dic_int["mean_flow"]

    results_df = results_to_df(data_dic_mix, names_dic,area,area_name,results_df)
    # results_df = results_to_df(data_dic_calc, names_dic,area,area_name,results_df)
    # results_df = results_to_df(data_dic_int, names_dic,area,area_name,results_df)
    if "DTA2" in slice_name:
        # import pdb; pdb.set_trace()
        print(data_dic_mix["mean_press"])
    return results_df, extractSurface, integrateVariables, [group_dsa_int,group_dsa_calc]

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
#region          DISPLAY FUNCTIONS
#==========================================

#Show all the slices
#Make an empty view before  
def show_slices(slices,slice_names,input_src=None,input_src_opacity=0.05,bar_range=[0,1],view_time=3,myview=None,mode="visible"): 
    if myview is not None:
        SetActiveView(myview) #Eventually make format_colorbar accept a view
    HideAll()
    myview = GetActiveView()
    myview.ViewTime = view_time
    print("View Time: {}".format(myview.ViewTime))
    use_tl = False 
    for key in slices.keys():
        if "TL_" in key:
            use_tl = True 
            break
    for slice_name in slice_names.keys(): 
        if use_tl:
            slice_name_tl = "TL_" + slice_name
            slice_name_fl = "FL_" + slice_name
        else:
            slice_name_tl = "LARGE_" + slice_name
            slice_name_fl = "SMALL_" + slice_name

        #mode fixed, bar range not used
        # disp = format_colorbar(slices[slice_name],bar_range=bar_range) 
        disp,disp_LUT = format_colorbar(slices[slice_name],bar_range=bar_range,mode="ignore")
        # disp = format_colorbar(slices[slice_name_tl],bar_range=bar_range,mode="ignore")
        # disp.Opacity = 0.2
        # if slices[slice_name_fl] is not None:
        #     disp = format_colorbar(slices[slice_name_fl],bar_range=bar_range,mode="ignore")
     
        # save_screen("centerline.png",image_res=[800,500])
        # import pdb; pdb.set_trace()

    #No need to rescale every time in for loop when not debugging 
    # Render()
    set_scale_range(disp,mode,bar_range,disp_LUT)
    Render()
    if input_src is not None:
        disp, _ = format_colorbar(input_src,opacity=input_src_opacity,mode="ignore",separate=True,show_bar=False,bar_range=bar_range)
        # set_solid_col(disp,opacity=0.2,color="grey")
        Render()
        # save_screen("centerline.png",image_res=[800,500])

#Save all time views
def save_all_times(save_path,slices,slice_names,input_src=None,input_src_opacity=0.05,mode="visible",bar_range=[0,1],image_res=[800,500],all_times=None,myview=None): 
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
        show_slices(slices,slice_names,input_src=input_src,input_src_opacity=input_src_opacity,mode=mode,bar_range=bar_range,view_time=view_time,myview=myview)
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
def filter_to_results(Input,src_type,slice_name,area,mean_press_name,mean_vel_mag_name):
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
def results_to_df(data,names,area,area_name,results_df=None):
    results_df_i = pd.DataFrame({area_name: [area]*len(data["mean_press"]),
                        names["mean_press"]: data["mean_press"],
                        names["mean_vel_mag"]: data["mean_vel_mag"],
                        names["mean_flow"]: data["mean_flow"]})

    if results_df is None:
        results_df = results_df_i 
    else:
        results_df = pd.concat([results_df,results_df_i],axis=1)

    return results_df 

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