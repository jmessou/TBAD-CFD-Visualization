#Functions that convert or compute some type of data

from paraview.simple import Calculator, ExtractTimeSteps, TemporalStatistics
import paraview 
from paraview.simple import *

from vtk.util.numpy_support import vtk_to_numpy
from vtk import vtkCellLocator
from vtk import mutable as vtk_mutable
import numpy as np

#Compute velocity in m/s and pressure in mmHg 
def convert_velocity_pressure(all_results_00,names={"velocity":"velocity","pressure":"pressure","Displacement":"displacement"}):
    # Velocity Calculator (m/s)
    velocityms = Calculator(registrationName='Velocity [m/s]', Input=all_results_00)
    velocityms.ResultArrayName = 'Velocity [m/s]'
    velocityms.Function = '{}/1000'.format(names["velocity"])

    # Pressure Calculator (mmHg)
    pressuremmHg = Calculator(registrationName='Pressure [mmHg]', Input=velocityms)
    pressuremmHg.ResultArrayName = 'Pressure [mmHg]'
    pressuremmHg.Function = '{}/133.3'.format(names["pressure"])

    return velocityms, pressuremmHg

#Compute WSS
def compute_wss_osi(all_results_00,max_step,first_step=0,names={"vWSS":"vWSS"}):
    # # Only select last results 
    # removeDudStep = ExtractTimeSteps(registrationName='RemoveDudStep', Input=all_results_00)
    # removeDudStep.SelectionMode = 'Select Time Range'
    # removeDudStep.TimeStepIndices = [first_step] #Not used based on the selection mode - leaving just in case 
    # removeDudStep.TimeStepRange = [first_step,max_step]
    src_in = all_results_00
    while hasattr(src_in,'Input'):
        src_in = src_in.Input
    print("TimestepValues: {}".format(src_in.TimestepValues))

    #Extracting first step no longer needed 
    removeDudStep = all_results_00 

    # Compute average of all fields
    avgofWSSvector = TemporalStatistics(registrationName='Avg of WSS (vector)', Input=removeDudStep)
    avgofWSSvector.ComputeMinimum = 0
    avgofWSSvector.ComputeMaximum = 0
    avgofWSSvector.ComputeStandardDeviation = 0

    # Get magnitude of WSS
    magnitudeofWSS = Calculator(registrationName='Magnitude of WSS', Input=removeDudStep)
    magnitudeofWSS.ResultArrayName = 'WSS_time'
    magnitudeofWSS.Function = 'mag({})'.format(names["vWSS"])

    # Take temporal average of WSS magnitude
    tAWSS = TemporalStatistics(registrationName='TAWSS (avg of mag)', Input=magnitudeofWSS)
    tAWSS.ComputeMinimum = 0
    tAWSS.ComputeMaximum = 0
    tAWSS.ComputeStandardDeviation = 0

    # Compute OSI: (1 - mag_of_avg / avg_of_mag)/2
    oSI = Calculator(registrationName='OSI ( 1 - mag_of_avg / avg_of_mag)/2', Input=tAWSS)
    oSI.ResultArrayName = 'OSI'
    oSI.Function = '(1-mag({}_average)/WSS_time_average)/2'.format(names["vWSS"])

    return tAWSS, oSI

#==========================================
#region          PRESSURE/FLOW FUNCTIONS
#==========================================

#Compute pressure and mag(mean(velocity)) using calculators
def compute_press_flow_python(slice_name,extractSurface, mean_press_name, mean_vel_name, mean_vel_mag_name,
    press_name='Pressure [mmHg]',vel_name='Velocity [m/s]'):

    #Compute Average Pressure at cross-section 
    pythonCalculator_press = PythonCalculator(registrationName='mean_press_' + slice_name, Input=extractSurface)
    pythonCalculator_press.Expression = "mean(inputs[0].PointData['{}'])".format(press_name)
    pythonCalculator_press.ArrayName = mean_press_name
    # pythonCalculator_press.CopyArrays = 0 
    
    #Compute Average Velocity at cross-section 
    pythonCalculator_flow = PythonCalculator(registrationName='mean_vel_' + slice_name, Input=pythonCalculator_press)
    pythonCalculator_flow.Expression = "mean(inputs[0].PointData['{}'],axis=0)".format(vel_name)
    pythonCalculator_flow.ArrayName = mean_vel_name
    
    #Take magnitude of average velocity
    #Take norm (not sure why Python Calculator does not accept axis for norm)
    Calculator_flow = Calculator(registrationName='mean_vel_mag_' + slice_name, Input=pythonCalculator_flow)
    Calculator_flow.ResultArrayName = mean_vel_mag_name
    Calculator_flow.Function = 'mag("{}")'.format(mean_vel_name)

    # print(sm.Fetch(Calculator_flow).GetPointData())

    return Calculator_flow

def compute_press_flow_integral(slice_name,integrateVariables, mean_press_name, mean_vel_mag_name,
    press_name='Pressure [mmHg]',vel_name='Velocity [m/s]'):

    #Take norm (not sure why Python Calculator does not accept axis for norm)
    Calculator_flow = Calculator(registrationName='mean_vel_mag_int_' + slice_name, Input=integrateVariables)
    Calculator_flow.ResultArrayName = mean_vel_mag_name
    Calculator_flow.Function = 'mag("{}")'.format(vel_name)

    #Just copy pressure 
    Calculator_press = Calculator(registrationName='mean_press_int_' + slice_name, Input=Calculator_flow)
    Calculator_press.ResultArrayName = mean_press_name
    Calculator_press.Function = '"{}"'.format(press_name)

    return Calculator_press

#endregion         END - PRESSURE/FLOW
#==========================================

#==========================================
#region         DISPLACEMENT FUNCTIONS
#==========================================
#Compute max(mag(displacement))  and mag(mean(displacement)) using calculators
def compute_disp_python(slice_name,extractSurface, mag_disp_name, mag_disp_max_name, mean_disp_name, mean_disp_mag_name,
    disp_name='Displacement'):

    #Take magnitude of Displacement
    #Take norm (not sure why Python Calculator does not accept axis for norm)
    Calculator_disp = Calculator(registrationName='mag_disp' + slice_name, Input=extractSurface)
    Calculator_disp.ResultArrayName = mag_disp_name
    Calculator_disp.Function = 'mag("{}")'.format(disp_name)

    #Compute Max of magnitude of Displacement at cross-section 
    pythonCalculator_disp = PythonCalculator(registrationName='mag_disp_max_' + slice_name, Input=Calculator_disp)
    pythonCalculator_disp.Expression = "max(inputs[0].PointData['{}'],axis=0)".format(mag_disp_name)
    pythonCalculator_disp.ArrayName = mag_disp_max_name
    
    #Compute Average Displacement at cross-section 
    pythonCalculator_disp2 = PythonCalculator(registrationName='mean_disp_' + slice_name, Input=pythonCalculator_disp)
    pythonCalculator_disp2.Expression = "mean(inputs[0].PointData['{}'],axis=0)".format(disp_name)
    pythonCalculator_disp2.ArrayName = mean_disp_name
    
    #Take magnitude of average Displacement
    #Take norm (not sure why Python Calculator does not accept axis for norm)
    Calculator_disp2 = Calculator(registrationName='mean_disp_mag_' + slice_name, Input=pythonCalculator_disp2)
    Calculator_disp2.ResultArrayName = mean_disp_mag_name
    Calculator_disp2.Function = 'mag("{}")'.format(mean_disp_name)

    # print(sm.Fetch(Calculator_flow).GetPointData())

    return Calculator_disp2

def compute_disp_integral(slice_name,integrateVariables, mean_disp_mag_name, disp_name='Displacement'):

    #Take norm (not sure why Python Calculator does not accept axis for norm)
    Calculator_disp = Calculator(registrationName='mean_disp_mag_int_' + slice_name, Input=integrateVariables)
    Calculator_disp.ResultArrayName = mean_disp_mag_name
    Calculator_disp.Function = 'mag("{}")'.format(disp_name)

    return Calculator_disp

#endregion         END - DISPLACEMENT
#==========================================

#==========================================
#region          HELPER FUNCTIONS
#==========================================

def threshold_id(Input,region_id,slice_name):
    #Use Threshold instead of query/extract
    threshold = Threshold(registrationName='Threshold_' + slice_name, Input=Input)
    threshold.Scalars = ['POINTS', 'RegionId']
    threshold.LowerThreshold = region_id
    threshold.UpperThreshold = region_id
    threshold.UpdatePipeline()
    return threshold

#Ideally, replace function above by this one
def threshold_face_id(Input,region_id,slice_name,array_name="ModelFaceID",cells=True,region_id_min=None,method=None,print_all=True,update_pipeline=True):
    if cells: 
        scalar_type = 'CELLS'
    else: 
        scalar_type = 'POINTS'

    if region_id_min is None and method is None: 
        region_id_min = region_id

    #Use Threshold instead of query/extract
    threshold = Threshold(registrationName='Threshold_' + slice_name, Input=Input)
    threshold.Scalars = [scalar_type, array_name]
    if region_id_min is not None:
        threshold.LowerThreshold = region_id_min
    if region_id is not None:
        threshold.UpperThreshold = region_id
    if method is not None:
        threshold.ThresholdMethod = method
    if print_all: 
        print(f"[Threshold] {slice_name} | {region_id_min} | {region_id} | {method} | {scalar_type} | {array_name} | Update {update_pipeline}")
    if update_pipeline:
        threshold.UpdatePipeline()
    return threshold

#Return area of surface
def get_area(Input):
    if Input is not None:
        # print(sm.Fetch(Input).GetPoints().GetNumberOfPoints())
        int1 = IntegrateVariables(Input=Input)
        area = sm.Fetch(int1).GetCellData().GetArray("Area").GetValue(0)
        # area = sm.Fetch(int1).GetPointData().GetArray("Area").GetValue(0)
        # print(f"Area Cell Data: {area1} | Area Point Data: {area}")
    else: 
        # print("Input is None")
        area = 0
    return area

#Return length of surface
def get_length(Input):
    if Input is not None:
        # print(sm.Fetch(Input).GetPoints().GetNumberOfPoints())
        int1 = IntegrateVariables(Input=Input)
        length = sm.Fetch(int1).GetCellData().GetArray("Length").GetValue(0)
        # length = sm.Fetch(int1).GetPointData().GetArray("Area").GetValue(0)
        # print(f"Length Cell Data: {length1} | Length Point Data: {length}")
    else: 
        # print("Input is None")
        length = 0
    return length

#Return the coordinates of the centroid of the points in the given source and the area
##the total number of points, and the first point in the array
def get_centroid(src):
    integral = IntegrateVariables(Input=src)
    integral_vtk = sm.Fetch(integral)
    try:
        area = integral_vtk.GetCellData().GetArray("Area").GetValue(0)
    except: 
        area = integral_vtk.GetCellData().GetArray("Volume").GetValue(0)
    coords_centroid = vtk_to_numpy(integral_vtk.GetPoints().GetData())

    # coords_i, pointCoords_i = get_vtk(src)
    # coords_centroid = (coords_i[0]['center'],coords_i[1]['center'],coords_i[2]['center'])
    # nb_points = pointCoords_i.shape[0]
    # return coords_centroid, nb_points, pointCoords_i[0]

    return coords_centroid[0], area  

#Add/Copy an array using the calculator 
def add_array(obj_src,reg_name,res_name,func,res_type='Int',set_active=True):
    # Create DOMAIN_ID using calculator
    calculator1 = Calculator(registrationName=reg_name, Input=obj_src)
    calculator1.AttributeType = 'Cell Data'
    calculator1.ResultArrayName = res_name
    calculator1.Function = str(func)
    calculator1.ResultArrayType = res_type

    if set_active:
        #Compute value
        SetActiveSource(calculator1)
    return calculator1    

#Check if elements of id_to_check are in the "GlobalElementID" array of slice_src
def is_id_in_src(slice_src, ids_to_check,print_warning=True): 
    id_in_slice = False 
    slice_ids = vtk_to_numpy(sm.Fetch(slice_src).GetCellData().GetArray("GlobalElementID"))
    id_in_slice = np.isin(ids_to_check,slice_ids)
    print(f"--{len(slice_ids)} elements in slice")

    if print_warning and np.all(id_in_slice): 
        print(f"[WARNING!] All IDs ({ids_to_check}) are in slice")
    return id_in_slice, slice_ids

#Find the closest cell to a point 
def find_closest_cell_to_point(src,ref_point,print_all=True):
    # Create a vtkCellLocator
    cell_locator = vtkCellLocator()
    cell_locator.SetDataSet(sm.Fetch(src))
    cell_locator.BuildLocator()

    # Need to store some info for vtk function 
    closest_cell_id = vtk_mutable(0)
    subId = vtk_mutable(0)
    dist2 = vtk_mutable(0.0)
    cp = [0.0, 0.0, 0.0]
    # Find the closest cell
    cell_locator.FindClosestPoint(ref_point, cp,closest_cell_id, subId,dist2)
    if print_all:
        print(f"Reference point: {ref_point}")
        print(f"Closest point: {cp}")
        print(f"Closest cell ID: {closest_cell_id}")
        print(f"Distance to closest cell: {dist2}")

    return cp, closest_cell_id, subId, dist2

#Find the closest point in points from the reference point
#Return minimum distance and the index of the closest point
def find_closest_point_np(ref_point,points,print_all=True): 
    if ref_point is None:
        print("No reference point, returning None.")
        return None, None, None

    ref_point = np.array(ref_point)
    # Replace None with [np.Inf,np.Inf,np.Inf]
    points = np.array([p if p is not None else [np.Inf] * len(ref_point) for p in points])
    
    if not np.any(points):
        print("[WARNING!] All points are None. Cannot compute distances.")
        return None, None, None

    # Calculate the distances from the reference point to each point
    distances = np.linalg.norm(points - ref_point, axis=1)

    # Find the smallest distance and its index
    min_index = np.argmin(distances)
    min_distance = distances[min_index]
    distances = [d if d is not np.Inf else None for d in distances]
    if print_all: 
        print(f"Distances: {distances} (min_index: {min_index})")
    return min_distance, min_index, distances
#endregion         END - HELPER
#==========================================