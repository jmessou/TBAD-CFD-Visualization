#Functions that convert or compute some type of data

from paraview.simple import Calculator, ExtractTimeSteps, TemporalStatistics
import paraview 
from paraview.simple import *


#Compute velocity in m/s and pressure in mmHg 
def convert_velocity_pressure(all_results_00):
    # Velocity Calculator (m/s)
    velocityms = Calculator(registrationName='Velocity [m/s]', Input=all_results_00)
    velocityms.ResultArrayName = 'Velocity [m/s]'
    velocityms.Function = 'velocity/1000'

    # Pressure Calculator (mmHg)
    pressuremmHg = Calculator(registrationName='Pressure [mmHg]', Input=velocityms)
    pressuremmHg.ResultArrayName = 'Pressure [mmHg]'
    pressuremmHg.Function = 'pressure/133.3'

    return velocityms, pressuremmHg

#Compute WSS
def compute_wss_osi(all_results_00,max_step,first_step=0):
    # # Only select last results 
    # removeDudStep = ExtractTimeSteps(registrationName='RemoveDudStep', Input=all_results_00)
    # removeDudStep.SelectionMode = 'Select Time Range'
    # removeDudStep.TimeStepIndices = [first_step] #Not used based on the selection mode - leaving just in case 
    # removeDudStep.TimeStepRange = [first_step,max_step]
    print("TimestepValues: {}".format(all_results_00.TimestepValues))

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
    magnitudeofWSS.Function = 'mag(vWSS)'

    # Take temporal average of WSS magnitude
    tAWSS = TemporalStatistics(registrationName='TAWSS (avg of mag)', Input=magnitudeofWSS)
    tAWSS.ComputeMinimum = 0
    tAWSS.ComputeMaximum = 0
    tAWSS.ComputeStandardDeviation = 0

    # Compute OSI: (1 - mag_of_avg / avg_of_mag)/2
    oSI = Calculator(registrationName='OSI ( 1 - mag_of_avg / avg_of_mag)/2', Input=tAWSS)
    oSI.ResultArrayName = 'OSI'
    oSI.Function = '(1-mag(vWSS_average)/WSS_time_average)/2'

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
#region          HELPER FUNCTIONS
#==========================================

def threshold_id(Input,region_id,slice_name):
    #Use Threshold instead of query/extract
    threshold = Threshold(registrationName='Threshold_' + slice_name, Input=Input)
    threshold.Scalars = ['POINTS', 'RegionId']
    threshold.LowerThreshold = region_id
    threshold.UpperThreshold = region_id

    return threshold

#Return area of surface
def get_area(Input):
    int1 = IntegrateVariables(Input=Input)
    return sm.Fetch(int1).GetCellData().GetArray("Area").GetValue(0)

    
#endregion         END - HELPER
#==========================================