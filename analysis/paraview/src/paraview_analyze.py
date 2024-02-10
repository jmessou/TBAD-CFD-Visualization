import paraview
from paraview.simple import *
from paraview.simple import _DisableFirstRenderCameraReset

from paraview.vtk.numpy_interface import dataset_adapter as dsa

from paraview_utils import *
from paraview_slices import *
from paraview_display import * 

#### disable automatic camera reset on 'Show'
_DisableFirstRenderCameraReset()

import vtk 
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np 

import pdb 
import os 

#from paraview_local import *

import argparse
import subprocess
from _paraview_config import ConfigAll



def set_up_view(coords,suffix="",assign_layout=True,view_dic={"az":-35,"el":0}, analysis_time_step=6):
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

 
 

def get_vtk(pv_vtu_src):
    reader = vtk.vtkXMLUnstructuredGridReader()                                  
    reader.SetFileName(pv_vtu_src.FileName[0])                                                 
    reader.Update()                                                                                                                                                            
    mesh = reader.GetOutput()                                                                                                     
    cells = mesh.GetCells()                                                      
    nCells = cells.GetNumberOfCells()                                                                                                       
    cellConns = vtk_to_numpy(cells.GetConnectivityArray())                       
    cellOffsets = vtk_to_numpy(cells.GetOffsetsArray())                          
    cellTypes = vtk_to_numpy(mesh.GetCellTypesArray())                           
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

# Clean pipeline
def clean_pipe():
    for x in GetLayouts().values():
        Delete(x)
    for x in GetSources().values():
        Delete(x[0])


def process_paraview(config_file='example.ini',patient=None,use_diff_max=False,stream_nb_pts=None,def_view=None):
    print_header("Read Config File...")
    Config1 = ConfigAll(config_file) 
    file_range, loop_timestep, _, nb_files, file_stride =  Config1.get_analysis()
    
    #Get Patient from ini file if None
    if patient is None: 
        patient = Config1.get_general()

    #Loop through each file separately 
    if loop_timestep and file_range is not None and Config1.fc.process_slice_bool:
        append_res = False
        nb_files = None
        print_header("[loop_timestep] Processing each file separately")
        for timestep in range(file_range[0],file_range[1],file_stride):
            file_range_i = (timestep,timestep)
            print("Range: {}".format(file_range_i))
            pv_analyze(patient,Config1,use_diff_max=use_diff_max,stream_nb_pts=stream_nb_pts,file_range=file_range_i,nb_files=nb_files,append_res=append_res,def_view=def_view)
            ResetSession()
            Connect()
            append_res = True
        print("[REMINDER] loop_timestep = True, all visualization was disabled.")
    else: #Process all files together (necessary for OSI/TAWSS, should be used for other visualization as well)
        pv_analyze(patient,Config1,use_diff_max=use_diff_max,stream_nb_pts=stream_nb_pts,file_range=file_range,nb_files=nb_files,def_view=def_view)
        ResetSession()
        Connect()
    #return
    #pv_analyze(patient,palette="cool",use_diff_max=use_diff_max)


#TO REMOVE
def print_vtk(source):
    mesh = paraview.servermanager.Fetch(source) 
    print(mesh)
    #import pdb; pdb.set_trace()  
    #mesh_np = vtk_to_numpy(mesh.GetCellData().GetArray("CenterlineIds"))            
    print(mesh.GetCellData().GetArray("Length").GetValue(0)) 
    print(mesh.GetCellData().GetArray("CenterlineIds"))                                                                                     
    import pdb; pdb.set_trace()
    cells = mesh.GetCells()                                                      
    nCells = cells.GetNumberOfCells()                                                                                                       
    cellConns = vtk_to_numpy(cells.GetConnectivityArray())                       
    cellOffsets = vtk_to_numpy(cells.GetOffsetsArray())                          
    cellTypes = vtk_to_numpy(mesh.GetCellTypesArray())                           
    pointCoords = vtk_to_numpy(mesh.GetPoints().GetData())
    print(cellTypes)
    print(len(pointCoords))
    import pdb; pdb.set_trace()


def load_true_lumen(tl_path):
    print("\nLoading {}".format(tl_path))
    # create a new 'XML PolyData Reader'
    truelumenvtp = XMLPolyDataReader(registrationName='true-lumen.vtp', FileName=[tl_path])
    truelumenvtp.CellArrayStatus = ['BoundaryCells', 'CapID', 'BadTriangle', 'FreeEdge', 'BooleanRegion', 'GlobalBoundaryCells', 'ModelFaceID', 'Normals', 'ActiveCells']
    truelumenvtp.PointArrayStatus = ['BoundaryPoints', 'GlobalBoundaryPoints']

    return truelumenvtp


#Find index of current time step - might be a better way
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

#Clip entry tear forthe given source and add to view if provided 
def clip_box(src,coords,opacity=None,myview=None,box_type="entry"): 
    entry = Clip(registrationName=box_type, Input=src)
    entry.ClipType = 'Box'
    # entry.Scalars = ['POINTS', 'pressure']
    
    # init the 'Box' selected for 'ClipType'
    print("Box type: {}".format(box_type))
    if box_type == "exit":
        entry.ClipType.Position = [coords[0]["min"], coords[1]["min"], coords[2]["q1"]+0.42*(coords[2]["diff"]/4)] #Bottom square set at bottom quarter
        entry.ClipType.Length = [coords[0]["diff"]+2, coords[1]["diff"]+2, coords[2]["diff"]/4]
    elif box_type == "fenestration":
        entry.ClipType.Position = [coords[0]["min"], coords[1]["min"], coords[2]["q3"]*1.6] #Bottom square set at TBD
        entry.ClipType.Length = [coords[0]["diff"]+2, coords[1]["diff"]+2, coords[2]["diff"]/4*1.6]
    elif box_type == "thrombus":
        entry.ClipType.Position = [coords[0]["min"]+0.625*coords[0]["diff"], coords[1]["min"], coords[2]["q3"]+0.05*(coords[2]["diff"]/4)] #Bottom square set at top quarter
        entry.ClipType.Length = [coords[0]["diff"]+2, coords[1]["diff"]+2, coords[2]["diff"]/4*0.75]
    else:
        entry.ClipType.Position = [coords[0]["min"], coords[1]["min"], coords[2]["q3"]] #Bottom square set at top quarter
        entry.ClipType.Length = [coords[0]["diff"]+2, coords[1]["diff"]+2, coords[2]["diff"]/4]
    print("Box length: {}".format(entry.ClipType.Length))
    if myview is not None: 
        # Add entry tear to view 2
        entryDisplay = show_and_color(entry,myview,opacity,show_bar=False)
    return entry, entryDisplay

def do_slice_save(fc,view_dic,input_src_slices,patient,stats,slices_out_dic,slices_load_out_dic,coords,
        all_times,out_dir=".",image_res=[800,500],analysis_time_step=6,wall_opacity=0.05,
        cb_mode="visible",bar_range=[0,0.1],backgrounds=None,suffix=""):
    if cb_mode == "bar_range":
        suffix += "_" + str(bar_range[1]).replace(".","_")
    else:
        suffix += ""
    slices_out = os.path.join(out_dir,"slices_p{}{}.png".format(patient,suffix))
    slices_out_load = os.path.join(out_dir,"slices_p{}_fixed{}.png".format(patient,suffix))
    centerline_out = os.path.join(out_dir,"centerline_p{}.png".format(patient))

    if fc.save_slices_bool:
        print_header("Saving Slices")
        if fc.slice_centerline_bool: 
            renderView11, layout11 = set_up_view(coords,"Slices",view_dic=view_dic[patient[0]],analysis_time_step=analysis_time_step)
            if fc.save_all_times_bool:
                save_all_times(slices_out,slices_out_dic["slices"],slices_out_dic["slice_names"],input_src=input_src_slices,image_res=image_res,mode=cb_mode,bar_range=bar_range,all_times=all_times,myview=renderView11)
            #Main Timestep (systole)
            show_slices(slices_out_dic["slices"],slices_out_dic["slice_names"],input_src=input_src_slices,input_src_opacity=wall_opacity,mode=cb_mode,bar_range=bar_range,view_time=analysis_time_step,myview=renderView11)
            save_screen(add_suffix(slices_out,"_step{}".format(analysis_time_step)),backgrounds=backgrounds,image_res=image_res,myview=renderView11)
        
        if fc.load_and_slice_centerline_bool:
            renderView12, layout12 = set_up_view(coords,"Slices_load",view_dic=view_dic[patient[0]],analysis_time_step=analysis_time_step)
            if fc.save_all_times_bool:
                save_all_times(slices_out_load,slices_load_out_dic["slices"],slices_load_out_dic["slice_names"],input_src=input_src_slices,image_res=image_res,mode=cb_mode,bar_range=bar_range,all_times=all_times,myview=renderView12)
            #Main Timestep (systole)
            show_slices(slices_load_out_dic["slices"],slices_load_out_dic["slice_names"],input_src=input_src_slices,input_src_opacity=wall_opacity,mode=cb_mode,bar_range=bar_range,view_time=analysis_time_step,myview=renderView12)
            save_screen(add_suffix(slices_out_load,"_step{}".format(analysis_time_step)),backgrounds=backgrounds,image_res=image_res,myview=renderView12)
            
    #Show and save centerline
    if fc.save_centerline_bool: 
        print_header("Saving CenterLine")
        renderView10, layout10 = set_up_view(coords,"Centerline")
        SetActiveView(renderView10)
        # Show(stats["centerline"]); save_screen(centerline_out,myview=renderView10); Hide(stats["centerline"]);
        # Show(stats["longest_centerline"]); save_screen(add_suffix(centerline_out,"_longest"),myview=renderView10); Hide(stats["longest_centerline"])
        select_centerline(stats["centerline"],stats["longest_centerline_id"],print_all=True,save_path=add_suffix(centerline_out,"_longest"))
    
def do_slice_work(fc,input_src_slices,patient,stats,slices_tsv,out_dir="/out",analysis_time_step=6,wall_opacity=0.05,bar_range=[0,0.9],append_res=False,stride=10):
    box_size = 50
    z_normal = False
    force_normal=10000 # Will set the normal to the zaxis for every slice after this point
    # show = True
    show = False
    # stride = 100
    erase_range = [50,150]
    
    os.makedirs(out_dir,exist_ok=True)
    results_tsv=os.path.join(out_dir,"results_p{}.tsv".format(patient))
    #load_and_slice will NOT overwrite this file, remove it manually
    fl_id_tsv=os.path.join(out_dir,"false_lumen_id_p{}.tsv".format(patient))

    if show: #Need a view
        coords, _ = get_vtk(stats["all_results_00"]) 
        renderView1, layout1 = set_up_view(coords,analysis_time_step=analysis_time_step)
        camera=GetActiveCamera()
        camera.Elevation(45)
        GetActiveView().ViewTime = analysis_time_step 

    slices_out_dic = {"slices":None, "slice_names":None, "results_df":None}
    slices_load_out_dic = slices_out_dic.copy()

    #=============================
    #       Make and Analyze Slices (no GUI needed if show is False)
    #       If show is True, a view needs to be already created
    #=============================
    if fc.slice_centerline_bool:
        print_header("Slice CenterLine")
        results_tsv_i = results_tsv
        print("Results tsv: {}".format(results_tsv_i))
        slices_out_dic["slices"], slices_out_dic["slice_names"], slices_out_dic["results_df"] = slice_centerline(input_src_slices,stats["longest_centerline"],slices_tsv=add_suffix(slices_tsv,"_all"),results_tsv=results_tsv_i,
                stride=stride,box_size=box_size,z_normal=z_normal,show=show,force_normal=force_normal,erase_range=erase_range,fl_id_tsv=fl_id_tsv,process_slice_bool=fc.process_slice_bool,append_res=append_res)

    if fc.load_and_slice_centerline_bool:
        print_header("Load and Slice CenterLine")
        results_tsv_i = add_suffix(results_tsv,"_fixed")
        print("Results tsv: {}".format(results_tsv_i))
        slices_load_out_dic["slices"], slices_load_out_dic["slice_names"], slices_load_out_dic["results_df"] = load_and_make_slices(input_src_slices,slices_tsv,results_tsv=results_tsv_i,box_size=box_size,
                z_normal=z_normal,show=show,force_normal=force_normal,fl_id_tsv=fl_id_tsv,process_slice_bool=fc.process_slice_bool,append_res=append_res)

    # import pdb; pdb.set_trace()
    # if load_state(patient) != 0 :
    #     return 

    #=============================
    #       Convert Results from LARGE/SMALL to TL/FL
    #=============================
    if fc.convert_results_bool: 
        print_header("Convert Results")
        convert_slice_names_dir(os.path.dirname(add_suffix(results_tsv,"_fixed")))
        #convert_slice_names(add_suffix(results_tsv,"_fixed"),add_suffix(fl_id_tsv,"_fixed"))

    return slices_out_dic, slices_load_out_dic 

#Show OSI and TAWSS
def show_osi_tawss(coords,stats,cb_mode_osi,cb_mode_tawss,bar_range_tawss,view_dic,patient,debug=False):
    # #=============================
    # #       OSI - View 15
    # #=============================
    renderView15, _= set_up_view(coords,"15_OSI",view_dic=view_dic[patient[0]])
    format_colorbar(stats["osi"],mode=cb_mode_osi,color_by_array="OSI",title="OSI",debug=debug)
    ResetCamera()  

    # #=============================
    # #       OSI - View 15_B
    # #=============================
    renderView15_B, _= set_up_view(coords,"15_OSI_B",view_dic=view_dic[patient[0]])
    format_colorbar(stats["osi"],mode=cb_mode_osi,color_by_array="OSI",title="OSI",debug=debug)
    ResetCamera()  
    GetActiveCamera().Azimuth(180)
    Render()

    #=============================
    #       TAWSS - View 16
    #=============================
    renderView16, _= set_up_view(coords,"16_TAWSS",view_dic=view_dic[patient[0]])
    format_colorbar(stats["tAWSS"],bar_range=bar_range_tawss,mode=cb_mode_tawss,color_by_array="WSS_time_average",title="TAWSS (Pa)",debug=debug)
    ResetCamera() 

    #=============================
    #       TAWSS - View 16_B
    #=============================
    renderView16_B, _= set_up_view(coords,"16_TAWSS_B",view_dic=view_dic[patient[0]])
    format_colorbar(stats["tAWSS"],bar_range=bar_range_tawss,mode=cb_mode_tawss,color_by_array="WSS_time_average",title="TAWSS (Pa)",debug=debug)
    ResetCamera() 
    GetActiveCamera().Azimuth(180)
    Render() 

    return renderView15, renderView15_B, renderView16, renderView16_B

def show_model(coords,all_results_00,truelumenvtp,cross_sec,view_dic,patient,thrombosedvtp=None): 
    #=============================
    #      True lumen model and Combined model - View 4
    #=============================
    renderView4, _= set_up_view(coords,"4_Model",view_dic=view_dic[patient[0]])
    #True Lumen
    tlDisplay = Show(truelumenvtp, renderView4, 'UnstructuredGridRepresentation')
    tlDisplay = set_solid_col(tlDisplay,opacity=0.3,color="purple")
    #Combined Model
    cbDisplay = Show(all_results_00, renderView4, 'UnstructuredGridRepresentation')
    cbDisplay = set_solid_col(cbDisplay,opacity=0.2,color="grey")
    if thrombosedvtp != None:
        thrombus, thDisplay = clip_box(thrombosedvtp, coords, myview=renderView4,box_type="thrombus")
        # thDisplay = Show(thrombus, renderView4, 'UnstructuredGridRepresentation')
        thDisplay = set_solid_col(thDisplay,opacity=0.1,color="grey")
        # thDisplay = set_solid_col(thDisplay,opacity=0.2,color="blue")

    #All cross-sections 
    for cross in cross_sec:
        src_cross = FindSource(cross)
        if src_cross is not None:
            print("Showing: {}".format(cross))
            crossDisplay = Show(src_cross, renderView4, 'GeometryRepresentation')
            # crossDisplay = set_solid_col(crossDisplay,opacity=1,color="green",rep="Wireframe")
            crossDisplay = set_solid_col(crossDisplay,opacity=1,color="green",disable_light=True)
    # GetActiveCamera().Azimuth(70)
    # GetActiveCamera().Zoom(3)
    Render()
    # import pdb; pdb.set_trace()

    return renderView4

#Return the bar_range and the name of the folder where the results will be saved 
#Set use_diff_max to False if you're using the config file 
def get_bar_range_vel(patient,bar_range_vel_max,use_diff_max):
    if use_diff_max: 
        y_maxs = {"1": 1, "2": 1.6, "3": 1.4}
        velmax_folder = "velmax_diff"
    else: #Set use_diff_max to False if you're using the config file 
        y_max = bar_range_vel_max
        y_maxs = {"1":y_max, "2": y_max, "3":y_max}
        velmax_folder = "velmax_{}".format(str(y_max).replace(".","_"))

    bar_range=[0,y_maxs[patient[0]]]  #Only used when mode is bar_range

    return bar_range, velmax_folder

#Add parameters to suffix 
def add_params_to_suffix(color_bar_mode_vel, bar_range, view_dic, patient):
    if color_bar_mode_vel == "bar_range":
        suffix = "_" + str(bar_range[1]).replace(".","_")
    else:
        suffix = ""
    if view_dic[patient[0]]["el"] != 45:
        suffix = suffix + "_el{}".format(view_dic[patient[0]]["el"])
    if view_dic[patient[0]]["az"] != -35:
        suffix = suffix + "_az{}".format(view_dic[patient[0]]["az"])

    return suffix 

# GetSources().keys(): all sources 
def pv_analyze(patient,Config1,use_diff_max=False,stream_nb_pts=None,out_dir=None,file_range=(0,None),nb_files=None,append_res=False,def_view=None,debug=False):
    #=============================
    #       Clean Pipeline and load flow control 
    #=============================
    print("Cleaning pipeline...")
    clean_pipe()

    #Controls what should be run
    print_header("Flow Control")
    fc = Config1.fc
    print(fc)

    #==========================================
    #       Parameters - Set using config file (.ini)
    #==========================================
    #PATHS, FILE READER + ANALYSIS, SLICE WORK
    centerline_path, vtu_folder, slices_tsv, slices_tsv_all, out_dir_conf, tl_path, etfn_path = Config1.get_paths()
    _, _, analysis_time_step, _, _ = Config1.get_analysis()
    slice_stride = Config1.get_slice() #If slice_centerline is used 
    if out_dir is None:
        out_dir = out_dir_conf

    #VIEW
    if def_view is None:
        az, el = Config1.get_view()
        def_view = {"az":az,"el":el}

    #VISUALIZATION: velocity/osi/tawss colorbar parameters
    palette, cb_mode_osi, color_bar_mode_vel, cb_mode_tawss, bar_range_tawss, bar_range_vel_max = Config1.get_vis()
    
    #SLICE + MODEL VISUALIZATION
    #cross_sec: list of cross section names
    wall_opacity, cross_sec = Config1.get_slice_vis() 
    show_thrombus = Config1.get_model_vis()

    #STREAMLINES parameters
    stream_radius, tube_radius, stream_max_length, stream_nb_pts_i, wall_opacity_vel = Config1.get_stream()
    if stream_nb_pts is None: 
        stream_nb_pts = stream_nb_pts_i

    #SAVING
    image_res = Config1.get_saving()

    #General 
    coords = None
    suffix = "" #Will be filled and added to files to keep track of parameters

    #Only used later if color_bar_mode_vel == "bar_range"
    bar_range, velmax_folder = get_bar_range_vel(patient,bar_range_vel_max,use_diff_max)
    print("color_bar_mode_vel: {} | bar_range velocity: {}| velmax_folder {}".format(color_bar_mode_vel,bar_range,velmax_folder))
    #=================END PARAMETERS=======================

    #Make Directory for all the results
    out_dir = os.path.join(out_dir,"patient{}".format(patient))
    os.makedirs(out_dir,exist_ok=True)
    print("Output Directory: {}".format(out_dir))
    
    #=============================
    #       Load Data and compute necessary info (converted pressure/velocity, centerline, OSI, tAWSS)
    #=============================
    all_results_00, stats = load_from_scratch(vtu_folder,centerline_path,file_range=file_range,nb_files=nb_files)
    if fc.save_model_bool:
        truelumenvtp = load_true_lumen(tl_path)  #Needed to show model 
        if "3F" in patient and show_thrombus:
            suffix += "_thrombus"
            print("Reading ET_FN model:{}".format(etfn_path))
            etfnvtp = load_true_lumen(etfn_path)
        else:
            etfnvtp = None

    all_times = all_results_00.TimestepValues
    if analysis_time_step > all_times[-1]:
        print("[WARNING] analysis_time_step = {} and last time step avaible = {}".format(analysis_time_step,all_times[-1]))
        analysis_time_step = all_times[-1]
        print("Set analysis_time_step = {}".format(analysis_time_step))

    # velocityms = FindSource("Velocity [m/s]")
    velocityms = stats["velocityms"]
    input_src_slices = stats["pressuremmHg"]
    SetActiveSource(input_src_slices)
    
    #=============================
    #       Slice work, see function
    #=============================
    slices_out_dic, slices_load_out_dic = do_slice_work(fc,input_src_slices,patient,stats,slices_tsv,out_dir,analysis_time_step,wall_opacity,bar_range,append_res,slice_stride)

    #=============================
    #       Parameters - Velocity
    #=============================
    
    #View 
    if def_view is None:
        def_view = {"az":-35,"el":45}
    view_dic = {"1": def_view, "2": def_view, "3": def_view}
    
    if True and view_dic[patient[0]]["el"] != 45:
        suffix = suffix + "_el{}".format(view_dic[patient[0]]["el"])
    if True and view_dic[patient[0]]["az"] != -35:
        suffix = suffix + "_az{}".format(view_dic[patient[0]]["az"])
        
    #Saving parameters
    backgrounds_slice = {"white":"WhiteBackground"}
    backgrounds_model = {"white":"WhiteBackground"}
    backgrounds = {"white":"WhiteBackground"}
    #Will iterate through all backgrounds 
    # backgrounds = {"white":"WhiteBackground", "black":"BlackBackground", "blueG":"BlueGrayBackground",
    #     "lightG":"LightGrayBackground", "trans":"transparent"} 
    
    # #=============================
    # #      Save Slices + Centerlines 
    # #=============================
    coords, _ = get_vtk(all_results_00) 
    do_slice_save(fc,view_dic,input_src_slices,patient,stats,slices_out_dic,slices_load_out_dic,
        coords,all_times,out_dir,image_res,analysis_time_step,wall_opacity,color_bar_mode_vel,bar_range,backgrounds_slice,suffix)

    #=============================
    #      True lumen model and Combined model - View 4
    #=============================
    if fc.save_model_bool: 
        print_header("Models")
        print(cross_sec)
        renderView4 = show_model(coords,all_results_00,truelumenvtp,cross_sec,view_dic,patient,etfnvtp)
    
    #=============================
    #       OSI and TAWSS - View 15, 15_B, 16, 16_B
    #=============================
    if fc.save_osi_tawss_bool:
        print_header("OSI and TAWSS")
        renderView15, renderView15_B, renderView16, renderView16_B = show_osi_tawss(coords,stats,cb_mode_osi,cb_mode_tawss,bar_range_tawss,view_dic,patient,debug=debug)

    #=============================
    #       Velocity Work (streamlines)
    #=============================
    if fc.do_streamlines_work_bool:
        print_header("Velocity Work")

        #Start streamlines from inlet
        np_centerline = vtk_to_numpy(sm.Fetch(stats["longest_centerline"]).GetPoints().GetData())
        stream_center = np_centerline[0]
        print("Center for streamlines: {}".format(stream_center))

        #=============================
        #      Change colorbar of velocity [m/s], range + blue to red rainbow
        #=============================    
        # import pdb; pdb.set_trace()
        renderView1, layout1 = set_up_view(coords,view_dic=view_dic[patient[0]],analysis_time_step=analysis_time_step)
        velocitymsDisplay,disp_LUT = format_colorbar(velocityms,bar_range=bar_range,mode=color_bar_mode_vel,debug=debug)
        velocitymsDisplay.Opacity = wall_opacity_vel
        
        #=============================
        #       Make streamlines
        #=============================
        
        # create a new 'Stream Tracer'
        streamlines_Np = StreamTracer(registrationName='Streamlines_Np', Input=velocityms,
            SeedType='Point Cloud')
        streamlines_Np.Vectors = ['POINTS', 'velocity']
        streamlines_Np.MaximumStreamlineLength = stream_max_length
        # init the 'Point Cloud' selected for 'SeedType'
        streamlines_Np.SeedType.Center = stream_center
        streamlines_Np.SeedType.NumberOfPoints = stream_nb_pts
        streamlines_Np.SeedType.Radius = stream_radius
        # Show(streamlines_Np)
        # Render()
        
        #=============================
        #       Add tube for thicker streamlines - View 1
        #=============================
        tube1 = Tube(registrationName='Tube1', Input=streamlines_Np)
        tube1.Scalars = ['POINTS', 'AngularVelocity']
        tube1.Vectors = ['POINTS', 'Normals']
        tube1.Radius = tube_radius
    
        tube1Display = Show(tube1, renderView1, 'UnstructuredGridRepresentation')
        tube1Display.SetScalarColoring("Velocity [m/s]", servermanager.GetAssociationFromString('POINTS'))
        tube1Display.SetScalarBarVisibility(renderView1, True)
        
        #Add to suffix parameters
        suffix = add_params_to_suffix(color_bar_mode_vel, bar_range, view_dic, patient)
        #====> Save main pictures (streamlines of the entire model)
        path_stream = os.path.join(out_dir,"main_streams",velmax_folder,"{}_{}p_step{}_stream{}.png".format(patient,str(stream_nb_pts),analysis_time_step,suffix))
        os.makedirs(os.path.dirname(path_stream),exist_ok=True)
        save_screen(path_stream,image_res,myview=renderView1)
        
        if fc.save_all_times_bool:
            save_all_times2(path_stream,input_src=tube1,image_res=image_res,all_times=all_times)
        
        # #=============================
        # #       Tube streamlines no wall - View 6
        # #=============================
        # renderView6, _= set_up_view(coords,"6_Tube_no_wall",view_dic=view_dic[patient[0]],analysis_time_step=analysis_time_step)
        # show_and_color(tube1,renderView6,show_bar=True)
        # ResetCamera()  

        #=============================
        #       Clip/crop box at the top - Entry Tear - View 2
        #=============================
        renderView2, _= set_up_view(coords,"2_Entry",view_dic=view_dic[patient[0]],analysis_time_step=analysis_time_step)
        entry_tube1, entry_tube1Display = clip_box(tube1, coords, myview=renderView2,box_type="entry")
        entry_velocityms, entry_velocitymsDisplay = clip_box(velocityms, coords, opacity=0.07,myview=renderView2,box_type="entry")
        ResetCamera()

        # #=============================
        # #       View with zoomed streamlines only (no walls) - View 3
        # #=============================
        # renderView3, _= set_up_view(coords,"3_No_wall",view_dic=view_dic[patient[0]],analysis_time_step=analysis_time_step)
        # show_and_color(entry_tube1,renderView3,show_bar=False)
        # ResetCamera() 

        #=============================
        #       Clip/crop box in the middle - Exit Tear - View 5
        #=============================
        if patient[0] == "3":
            #TODO: Change exit to fenetration eventually
            renderView5, _= set_up_view(coords,"5_Exit",view_dic=view_dic[patient[0]],analysis_time_step=analysis_time_step)
            exit_tube1, exit_tube1Display = clip_box(tube1, coords, myview=renderView5,box_type="fenestration")
            exit_velocityms, exit_velocitymsDisplay = clip_box(velocityms, coords, opacity=0.07,myview=renderView5,box_type="fenestration")
            ResetCamera()
        # elif patient[0] == "1":
        #     renderView5, _= set_up_view(coords,"5_Exit",view_dic=view_dic[patient[0]],analysis_time_step=analysis_time_step)
        #     exit_tube1, exit_tube1Display = clip_box(tube1, coords, myview=renderView5,box_type="exit")
        #     exit_velocityms, exit_velocitymsDisplay = clip_box(velocityms, coords, opacity=0.07,myview=renderView5,box_type="exit")
        #     ResetCamera()


    #=============================
    #       Save as picture
    #=============================
    print_header("Saving screenshots")
    if fc.do_streamlines_work_bool:
        path_stream = path_stream.replace("main_streams","all_streams")
        os.makedirs(os.path.dirname(path_stream),exist_ok=True)
        save_screen(path_stream,image_res,backgrounds,palette,myview=renderView1)
        save_screen(add_suffix(path_stream,"_entry"),image_res,backgrounds,palette,myview=renderView2)
        # save_screen(add_suffix(path_stream,"_entryTubes"),image_res,backgrounds,palette,myview=renderView3)
        # save_screen(add_suffix(path_stream,"_noWall"),image_res,backgrounds,palette,myview=renderView6)
        if fc.save_model_bool:
            save_screen(add_suffix(path_stream,"_model_tlcb"),image_res,backgrounds_model,palette,myview=renderView4)
    
    if fc.save_model_bool:
        model_path_2 = os.path.join(out_dir,"p{}.png".format(patient))
        save_screen(add_suffix(model_path_2,"_model_tlcb"+suffix),image_res,backgrounds_model,palette,myview=renderView4)
    
    if fc.save_osi_tawss_bool:
        backgrounds = {"trans":"transparent"}
        if cb_mode_tawss == "bar_range":
            suffix = "_" + str(bar_range_tawss[1]).replace(".","_")
        else:
            suffix = ""
        path_stream_osi_tawss = os.path.join(out_dir,"osi_tawss","p{}.png".format(patient))
        os.makedirs(os.path.dirname(path_stream_osi_tawss),exist_ok=True)
        save_screen(add_suffix(path_stream_osi_tawss,"_OSI"),image_res,backgrounds,palette,myview=renderView15)
        save_screen(add_suffix(path_stream_osi_tawss,"_OSI_180"),image_res,backgrounds,palette,myview=renderView15_B)
        save_screen(add_suffix(path_stream_osi_tawss,"{}_TAWSS".format(suffix)),image_res,backgrounds,palette,myview=renderView16)
        save_screen(add_suffix(path_stream_osi_tawss,"{}_TAWSS_180".format(suffix)),image_res,backgrounds,palette,myview=renderView16_B)
   
    # if fc.do_streamlines_work_bool and patient[0] in ["1"]:
    #     save_screen(add_suffix(path_stream,"_exitTubes"),image_res,backgrounds,palette,myview=renderView5)
    if fc.do_streamlines_work_bool and patient[0] in ["3"]:
        save_screen(add_suffix(path_stream,"_entryTubesFEN"),image_res,backgrounds,palette,myview=renderView5)
    
    #=============================
    #       Save state, so everything can be loaded into paraview
    #=============================
    print_header("Saving State File")
    path_state = os.path.join(out_dir,"p" + patient + "_FlowState.pvsm")
    # path_state = os.path.join(os.path.dirname(out_dir),"p" + patient + "_FlowState.pvsm")
    print(path_state)
    SaveState(path_state)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog="paraview_analyze.py",description='Processing + Visualization of Aortic Dissection CFD results')
    parser.add_argument('-c','--config', type=str,help="config file")
    args = parser.parse_args()
    print("Using config file: {}".format(args.config))

    process_paraview(args.config)

    
