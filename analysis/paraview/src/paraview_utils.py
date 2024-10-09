#Utility functions: centerline/loading/saving/printing 

#Depends on paraview calc 

import paraview 
from paraview.simple import *
from paraview.selection import _createSelection, _collectSelectionPorts

from paraview_calc import *

import glob 
import os 
import warnings 

#==========================================
#region          CENTERLINE FUNCTIONS
#==========================================

#Similar to paraview's but doesn't need a view to be created 
#Will add the selection to the Source 
def _QuerySelect(QueryString='', FieldType='POINT', Source=None, InsideOut=False):
    """Selection by query expression.
    - QueryString - string with NumPy-like boolean expression defining which attributes are selected
    - FieldType - atttribute to select, e.g., 'POINT' or 'CELL'
    - Source - if not set, then the selection will be on the active source
    - InsideOut - Invert the selection so that attributes that do not satisfy the expression are
      selected instead of elements that do
    """
    if not Source:
        Source = paraview.simple.GetActiveSource()

    if GetActiveView() is not None:
        repr = paraview.simple.GetRepresentation(Source)
        #paraview.simple.GetRepresentation(GetActiveSource())

    import paraview.vtk as vtk
    if GetActiveView() is not None:
        reprCollection = vtk.vtkCollection()
        reprCollection.AddItem(repr.SMProxy)

    # convert FieldType to ElementType. Eventually, all public API should change
    # to accepting ElementType but we'll do that after all selection sources use
    # ElementType consistently.
    ftype = vtk.vtkSelectionNode.GetFieldTypeFromString(FieldType)
    if ftype == vtk.vtkSelectionNode.NUM_FIELD_TYPES:
        raise ValueError("Invalid FieldType '%s'" % FieldType)
    ElementType = vtk.vtkSelectionNode.ConvertSelectionFieldToAttributeType(ftype)

    selection = _createSelection('SelectionQuerySource', ElementType=ElementType,
                                 QueryString=QueryString, InsideOut=InsideOut)
    if selection:
        if GetActiveView() is not None:
            selectionCollection = vtk.vtkCollection()
            selectionCollection.AddItem(selection.SMProxy)
            _collectSelectionPorts(reprCollection, selectionCollection)
        else:
            Source.SMProxy.SetSelectionInput(0, selection.SMProxy, 0)

    Source.UpdateVTKObjects()

#Find the longest centerline based on centerline_src (raw vtp from Simvascular)
#QuerySelect needs a view to be set, using _QuerySelect instead
#If the function is modified to save views, the view needs to be created outside of 
# the function first 
def get_longest_centerline(centerline_src,print_all=False):  

    centerline_len = 0
    centerline_id = -1
    centerline = None
    
    # Get the range of ids | Assume list goes from 0 to the max
    range_id = paraview.servermanager.Fetch(centerline_src).GetCellData().GetArray("CenterlineIds").GetRange(0)
    print("Centerline IDs Range: {}".format(range_id))
    
    # renderView1, layout1 = set_up_view(coords)
    #Extract a centerline and check if if it is the longest one
    for centerline_id_i in range(int(range_id[1]) + 1):
        # select_centerline(centerline_src,centerline_id_i,print_all=True,save_path="selection.png")
        select_centerline(centerline_src,centerline_id_i,print_all=print_all,save_path="")
        centerline_i = ExtractSelection() #Extract centerline

        #Get centerline length 
        integrateVariables_i = IntegrateVariables(registrationName='IntegrateVariables_center_' + str(centerline_id_i), Input=centerline_i)
        integrateVariables_vtk_i = paraview.servermanager.Fetch(integrateVariables_i) #Get the VTK Object
        centerline_len_i = integrateVariables_vtk_i.GetCellData().GetArray("Length").GetValue(0)
       
        if print_all:
            print("Length: {:.2f}".format(centerline_len_i))

        #Update centerline if it is longer 
        if centerline_len_i >= centerline_len: 
            centerline_len = centerline_len_i
            centerline = centerline_i 
            centerline_id = centerline_id_i

        SetActiveSource(centerline_src)
        ClearSelection() #Not needed, but keeping

        # import pdb; pdb.set_trace()

    #Final centerline
    print("Final centerline ID: {} | Length: {:.2f}".format(centerline_id_i,centerline_len_i))
    if GetActiveView() is not None:
        select_centerline(centerline_src,centerline_id,print_all=True,save_path="selection.png")
    
    return centerline, centerline_id

#Query and select a centerline based on a given ID
#Saves view if save_path is provided (view needs to be set before)
def select_centerline(centerline_src,centerline_id,print_all=True,save_path=""):
    query= "CenterlineIds == " + str(centerline_id) 
    if print_all:
        print("Query: {}".format(query))
    _QuerySelect(QueryString=query, FieldType="CELL", Source=centerline_src) #Select Centerline
    if save_path != "":
        SetActiveSource(centerline_src)
        Show(); save_screen(save_path); Hide()

#endregion          END - CENTERLINE FUNCTIONS
#==========================================

#==========================================
#region          LOADING FUNCTIONS
#==========================================

def load_from_scratch(vtu_folder,centerline_path=None,file_range=(0,None),file_stride=1,nb_files=None):
    print_header("Reading VTU folder")

    if os.path.isdir(vtu_folder):
        vtu_files=glob.glob(os.path.join(vtu_folder,"*.vtu"))
        vtu_files.sort()
        max_step = len(vtu_files)-1
        print("\nFound {} files in {}...".format(max_step+1,vtu_folder))
        if nb_files is not None and nb_files != max_step+1:
            file_stride = int((max_step)/(nb_files-1))
            print("Correcting file_stride to {} to have {} files".format(file_stride,nb_files))
        #Can only select files within a certain range (pass file_range based on number of cycles)
        vtu_files = restrict_files(vtu_files,file_range,file_stride=file_stride)
        print("Kept {} files".format(len(vtu_files)))
        print("First file: {}".format(os.path.basename(vtu_files[0])))
        print("Last file: {}".format(os.path.basename(vtu_files[-1])))
    else: 
        print("{} not found".format(vtu_folder))
        warnings.warn("{} not found".format(vtu_folder))
        return None, None, None

    # File reader 
    all_results_00 = XMLUnstructuredGridReader(registrationName='all_results_00*', FileName=vtu_files)
   
    #===ADDED #svFSI notation renamed with svsolver notation
    stats = {"all":None,"fluid":None,"solid":None,"flap":None}
    #Return dictionary of the names of the different fields used as written in the mesh
    array_names, fsi = get_field_names(all_results_00)

    if fsi:
        fluid_id = 1.0
        solid_id = 2.0
        wall_original_domain_id = 2
        flap_original_domain_id = 3
        flap_id = solid_id**flap_original_domain_id #Looks like this is how svfsi assigns the DOMAIN_ID
        wall_id = solid_id**wall_original_domain_id #Looks like this is how svfsi assigns the DOMAIN_ID
        all_results_00.CellArrayStatus = ['Domain_ID', 'Proc_ID', 'Mesh_ID', 'E_VonMises']
        all_results_00.PointArrayStatus = ['Velocity', 'Pressure', 'Displacement', 'VonMises_stress', 'WSS']
        all_results_00.TimeArray = 'None'
        #svFSI discards GlobalElementID, let's generate new ones so we can propagate FL/TL to other timesteps
        #Array name: GlobalCellIds, copied to GlobalElementID so type is int and it pass through filters
        all_results_00 = GenerateGlobalIds(registrationName='tmp_GenerateGlobalIds', Input=all_results_00)
        all_results_00 = add_array(all_results_00,'all_results_00_gids',"GlobalElementID","GlobalCellIds",res_type='Int',set_active=True)
        surface_first_id = print(f"GlobalElementID Range: {sm.Fetch(all_results_00).GetCellData().GetArray('GlobalElementID').GetValueRange()}")
        
        #Separate fluid and solid
        fluid_src = threshold_face_id(all_results_00,fluid_id,"fluid",array_name='Mesh_ID',cells=True)
        solid_src = threshold_face_id(all_results_00,solid_id,"solid",array_name='Mesh_ID',cells=True)
        fluid_src = WarpByVector(registrationName='fluid_warped', Input=fluid_src)
        fluid_src.Vectors = ['POINTS', 'Displacement']
        solid_src = WarpByVector(registrationName='solid_warped', Input=solid_src)
        solid_src.Vectors = ['POINTS', 'Displacement']
        
        #Convert velocity to m/s and pressure to mmHg
        #Compute TAWSS and OSI | Add centerline 
        stats["all"] = add_stats(all_results_00,max_step,array_names,centerline_path=centerline_path)
        stats["fluid"] = add_stats(fluid_src,max_step,array_names,centerline_path=None,stats_to_copy=stats["all"])
        stats["solid"] = add_stats(solid_src,max_step,array_names,centerline_path=None,stats_to_copy=stats["all"])
        if stats["solid"] is not None: #Extract Flap and Wall
            # get_array_names(sm.Fetch(solid_src),data_type="CellData",print_all=True)
            # get_array_names(sm.Fetch(solid_src),data_type="PointData",print_all=True)
            stats["solid"]["flap"] = extract_region(stats["solid"],region_name="flap",region_id=flap_id,array_name='Domain_ID')
            stats["solid"]["wall"] = extract_region(stats["solid"],region_name="wall",region_id=wall_id,array_name='Domain_ID')
        # print(dir(all_results_00))
        # import pdb; pdb.set_trace()
    #===END ADDED
    else:
        all_results_00.CellArrayStatus = ['GlobalElementID']
        all_results_00.PointArrayStatus = ['GlobalNodeID', 'pressure', 'velocity', 'vinplane_traction', 'vWSS', 'timeDeriv', 'average_speed', 'average_pressure']
        all_results_00.TimeArray = 'None'

        #Convert velocity to m/s and pressure to mmHg
        #Compute TAWSS and OSI | Add centerline 
        stats["fluid"] = add_stats(all_results_00,max_step,array_names,centerline_path=centerline_path)
        stats["all"] = stats["fluid"]
    return all_results_00, stats, array_names

#Extract the flap or the wall for all our main filters 
#Other regions can be extracted as well depending on DOMAIN_ID
def extract_region(solid_dic,region_name="flap",region_id=8,array_name='Domain_ID',cells=True): 
    region_dic = {}
    main_filters = ["pressuremmHg","osi","tAWSS"]
    main_filters_str = [f"{region_name}",f"osi",f"tAWSS"]
    array_names = [f"{array_name}",f"{array_name}_average",f"{array_name}_average"]
    update_pipelines = [True, False, False] #Don't update it for osi and TAWSS until we need it (reads all the files)
    for _filter_str, _filter, array_name, update_pipeline in zip(main_filters_str,main_filters,array_names, update_pipelines):
        region_dic[_filter_str] = threshold_face_id(solid_dic[_filter],region_id,_filter_str,array_name=array_name,cells=cells,update_pipeline=update_pipeline)
    return region_dic

#Convert velocity to m/s and pressure to mmHg
#Compute TAWSS and OSI 
def add_stats(src,max_step,array_names,centerline_path=None,stats_to_copy=None):
    #Convert velocity to m/s and pressure to mmHg
    #Compute TAWSS and OSI 
    stats = {}
    stats["all_results_00"] = src
    stats["velocityms"], stats["pressuremmHg"] = convert_velocity_pressure(src,names=array_names)
    stats["tAWSS"], stats["osi"] =  compute_wss_osi(src,max_step,first_step=1,names=array_names)
    
    #Add centerline if provided 
    if centerline_path is not None and os.path.isfile(centerline_path):
        print_header("Getting Longest Centerline")

        stats["centerline"] = XMLPolyDataReader(registrationName=os.path.basename(centerline_path), FileName=[centerline_path])
        stats["centerline"].CellArrayStatus = ['CenterlineIds', 'TractIds', 'Blanking', 'GroupIds']
        stats["centerline"].PointArrayStatus = ['MaximumInscribedSphereRadius']
        stats["centerline"].TimeArray = 'None'

        stats["longest_centerline"],stats["longest_centerline_id"] = get_longest_centerline(stats["centerline"])
    elif stats_to_copy is not None: 
        stats["centerline"] = stats_to_copy["centerline"]
        stats["longest_centerline"] = stats_to_copy["longest_centerline"]
        stats["longest_centerline_id"] = stats_to_copy["longest_centerline_id"]

    return stats 

#Return dictionary of the names of the different fields used as written in the mesh
#(only need to replace default ones when svFSI is used)
def get_field_names(src):
    fsi = False
    #Expected array names
    array_names = {"velocity":"velocity","pressure":"pressure","vWSS":"vWSS","displacement":"displacement"}
    #Array names to rename if svFSI and FSI used 
    arrays_to_rename = {"Velocity":"velocity","Pressure":"pressure","Displacement":"displacement","WSS":"vWSS"}

    src_vtk = paraview.servermanager.Fetch(src)
    point_data_array_names, _ = get_array_names(src_vtk,data_type="PointData",print_all=False)
    
    for key, val in arrays_to_rename.items():
        if key in point_data_array_names:
            print("Using {} as {}".format(val,key))
            array_names[val] = key
            # src_vtk.GetPointData().GetArray(key).SetName(val)
            fsi=True

    # point_data_array_names, _ = get_array_names(src_vtk,data_type="PointData",print_all=True)
    # cell_data_array_names, _ = get_array_names(src_vtk,data_type="CellData",print_all=True)
    
    return array_names, fsi

#Shorten list of files based on given range 
#The range is based on the file name 
#If you know the the indexes, no need for a special function
def restrict_files(vtu_files, file_range=(0, None),file_stride=1):

    start_idx = -1
    end_idx = -1
    if file_range[1] is None:
        file_range[1] =""
        end_idx = len(vtu_files)
    
    if file_range[0] > file_range[1]:
        print("[ERROR!] file_range MIN > file_range MAX: {}"
        "    Returning full range...".format(file_range))
        return vtu_files 

    #Get the number of digits for the format so we can pad the min/max 
    # that we are looking for accordingly i.e 10 ===> 00010 for 5 digts
    nb_digits = 0
    for i in vtu_files[-1].rsplit(".")[-2][::-1]: 
        if i == "_":
            break
        nb_digits += 1 


    for i, vtu_file in enumerate(vtu_files): 
        if "{idx:0{nb_dig}}".format(nb_dig=nb_digits,idx=file_range[0]) in vtu_file:
            start_idx = i  
        if "{idx:0{nb_dig}}".format(nb_dig=nb_digits,idx=file_range[1]) in vtu_file:
            end_idx =  i 
  
    if start_idx == -1 or end_idx == -1:
        print("[ERROR!] Could not find min or max file, given range: {} \n"
                "    Returning full range...".format(file_range))
        return vtu_files

    return vtu_files[start_idx:end_idx+1][::file_stride]

#Can use this for a unit test
#Open vtp and check which arrays are available
def load_vtp_and_get_arrays(fpath,vtu=False):
    vtp_fname = os.path.basename(fpath)
    if vtu: 
        vtp_src = XMLUnstructuredGridReader(registrationName=vtp_fname, FileName=fpath)
    else:
        vtp_src = XMLPolyDataReader(registrationName=vtp_fname, FileName=fpath)
    src_vtk = paraview.servermanager.Fetch(vtp_src)
    point_data_array_names, _ = get_array_names(src_vtk,data_type="PointData")
    cell_data_array_names, _ = get_array_names(src_vtk,data_type="CellData")

    print(point_data_array_names)
    print(cell_data_array_names)

    #vtp_src = XMLPolyDataReader(registrationName=vtp_fname, FileName=fpath)
    return vtp_src, point_data_array_names, cell_data_array_names
    
#endregion         END - LOADING FUNCTIONS
#==========================================

#==========================================
#region          SAVING FUNCTIONS
#==========================================

#Save current view or given view as png 
#image_res: Resolution (tuple)
#background: Disctionary where key is appended to the path and the value represents 
#   the paraview background option
#palette is a suffix added to the path
def save_screen(path_stream,image_res=(500, 500),backgrounds={"default":None},palette="",myview=None): 
    if myview is None:
        myview = GetActiveView()
    myview.OrientationAxesVisibility = 0
    Render()
    #SaveScreenshot(path_stream, myview,ImageResolution=image_res,OverrideColorPalette=background)
    print(f"Main path: {path_stream}")
    for back_key, background in backgrounds.items(): 
        path_screen = add_suffix(path_stream,"_" + palette + "_" + back_key)
        if back_key == "trans":
            SaveScreenshot(path_screen, myview,ImageResolution=image_res,OverrideColorPalette="WhiteBackground",TransparentBackground=True)
        else: 
            SaveScreenshot(path_screen, myview,ImageResolution=image_res,OverrideColorPalette=background)
    myview.OrientationAxesVisibility = 1

#Add suffix to a path
def add_suffix(path,suffix):
    path_split = path.rsplit(".",1)
    path = path_split[0] + suffix + "." + path_split[1]
    
    return path

#Make parent directory from a given path my_dir ==> parent/my_dir
def make_parent_dir(save_path,parent):
    save_path = os.path.join(os.path.dirname(save_path),parent,os.path.basename(save_path))
    os.makedirs(os.path.dirname(save_path),exist_ok=True)
    return save_path 

#endregion         END - SAVING FUNCTIONS
#==========================================

#==========================================
#region          ARRAY FUNCTIONS
#==========================================

#Get array names
#type=="PointData" OR type=="CellData"
def get_array_names(src_vtk,data_type="PointData",print_all=False):
    array_names = [] 

    if data_type=="PointData":
        data = src_vtk.GetPointData()
    elif data_type=="CellData":
        data = src_vtk.GetCellData() 
    else: 
        return array_names, None

    nb_arrays = data.GetNumberOfArrays()

    for i in range(nb_arrays):
        array_name = data.GetArrayName(i)
        array_names.append(array_name)

    if print_all:
        print("{} arrays: {}".format(data_type,array_names))

    return array_names, nb_arrays

#endregion         END - ARRAY FUNCTIONS
#==========================================

#==========================================
#region          PRINTING/HELPER FUNCTIONS
#==========================================

# _______________________
#
# |||||---EXAMPLE---|||||
# _______________________
def print_header(header):
	bar = (len(header)+16)*"_"
	print("\n{}\n".format(bar))
	print("|||||---{}---|||||".format(header))
	print("{}\n".format(bar))

#Switch the content of a and b 
def switch_var(a,b):
    tmp = b
    b = a 
    a = tmp 
    return a,b 

#Rename the keys in a dictionary using a mapping dictionary 
def map_fields(obj,mapping_dic):
    new_obj = {}
    for key in obj.keys():
        if key in mapping_dic:
            new_key = mapping_dic[key]
        else:
            new_key = key
        new_obj[new_key] = obj[key]
    return new_obj

#endregion          PRINTING/HELPER FUNCTIONS
#==========================================

#==========================================
#region          I/O
#==========================================
#Read folder or list
def read_data(data_path,extension=".dcm",pattern="*",recursive=False,no_error=False):
    if type(data_path) != type("str"):
        raise TypeError("data_path is not a string: {}".format(type(data_path)))
    elif os.path.isdir(data_path):
        flist = [path for path in glob.glob(f'{data_path}/{pattern}{extension}', recursive=recursive)]
    elif os.path.isfile(data_path):
        flist = [line.strip() for line in open(data_path, 'r')]
    elif no_error:
        flist = []
    else:
        raise ValueError("{} is neither a file or a directory".format(data_path))

    return flist

#Write a list to a file
def list_to_file(data_list,filepath):
    with open(filepath, 'w') as file_handler:
        for item in data_list:
            file_handler.write("{}\n".format(item))
#endregion          I/O
#==========================================

#==========================================
#region          TIMESTEP FUNCTIONS
#==========================================

# _______________________
#
# |||||---EXAMPLE---|||||
# _______________________

#Wrapper for GroupTimeSteps
#Will first make a block then call GroupTimeSteps
#Otherwise sm can't fetch, 
#based on quick search: issue comes from block name change? 
#Same issue when the python shell is used in Paraview
#but not in Paraview itself
#Note: Will need to call GetBlock twice: GetBlock(N).GetBlock(0)
def _GroupTimeSteps(registrationName, Input,block_names="timestep"):
	#Make Group
    groupDatasets = GroupDatasets(registrationName="Block_" + registrationName, Input=Input)
    groupDatasets.BlockNames = [block_names] 
    
    #Now calling GroupTimeSteps will work 
    return GroupTimeSteps(registrationName=registrationName, Input=groupDatasets)

#endregion          TIMESTEP FUNCTIONS
#==========================================