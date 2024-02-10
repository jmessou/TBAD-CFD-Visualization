#Utility functions: centerline/loading/saving/printing 

#Depends on paraview calc 

import paraview 
from paraview.simple import *
from paraview.selection import _createSelection, _collectSelectionPorts

from paraview_calc import *

import glob 
import os 

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
        return None, None

    # create a new 'XML Unstructured Grid Reader'
    all_results_00 = XMLUnstructuredGridReader(registrationName='all_results_00*', FileName=vtu_files)
    all_results_00.CellArrayStatus = ['GlobalElementID']
    all_results_00.PointArrayStatus = ['GlobalNodeID', 'pressure', 'velocity', 'vinplane_traction', 'vWSS', 'timeDeriv', 'average_speed', 'average_pressure']
    all_results_00.TimeArray = 'None'

    #Convert velocity to m/s and pressure to mmHg
    #Compute TAWSS and OSI 
    stats = {}
    stats["all_results_00"] = all_results_00
    stats["velocityms"], stats["pressuremmHg"] = convert_velocity_pressure(all_results_00)
    stats["tAWSS"], stats["osi"] =  compute_wss_osi(all_results_00,max_step,first_step=1)
    
    #Add centerline if provided 
    if centerline_path is not None and os.path.isfile(centerline_path):
        print_header("Getting Longest Centerline")

        stats["centerline"] = XMLPolyDataReader(registrationName='p1_fluid_v1_centerline.vtp', FileName=[centerline_path])
        stats["centerline"].CellArrayStatus = ['CenterlineIds', 'TractIds', 'Blanking', 'GroupIds']
        stats["centerline"].PointArrayStatus = ['MaximumInscribedSphereRadius']
        stats["centerline"].TimeArray = 'None'

        stats["longest_centerline"],stats["longest_centerline_id"] = get_longest_centerline(stats["centerline"])

    return all_results_00, stats

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