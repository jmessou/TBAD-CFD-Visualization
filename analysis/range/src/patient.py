import warnings
import os
import pandas as pd
import numpy as np
from scipy.signal import argrelextrema

def format_true(bool_val):
    if bool_val is True:
        return  "<<<<====TRUE====>>>>"
    else: 
        return bool_val

#Mean Blood Pressure
def compute_meanbp(sys_bp,dias_bp):
    return (1/3)*(sys_bp+2*dias_bp)
    

#Patient class with clinical values and experiment results 
#Can analyze each region separately (region = column in result DF)
#Uses Region
class Patient:
    def __init__(self, patient_id, job_dir="",sys_bp=0, dias_bp=0, heart_rate=0):
        self.id=patient_id
        self.job_dir = job_dir

        #Can be stored using load_from_file
        self.sys_bp = sys_bp 
        self.dias_bp = dias_bp 
        self.mean_bp = compute_meanbp(self.sys_bp,self.dias_bp)
        self.amp = self.compute_bpamplitude()
        self.heart_rate = heart_rate 

        #Experiment info (add using add_exp_info)
        self.exp = ""
        self.inflow_s = -1
        self.R = -1 
        self.C = -1 
        self.exp_dic = None

        #Experiment variables
        #You can use the same format as long as region is present in columns of res_dir
        self.nb_cycles = 0 
        self.inlet = Region("inlet") 
        self.all = Region("All")
        self.inlet_res_dic = None 
        self.all_res_dic = None
        

        #Blood Pressure Amplitude
    def compute_bpamplitude(self):
        return (self.sys_bp-self.dias_bp)

    #Read data (Systolic/Diastolic BP and HR) from tsv file 
    def load_from_file(self,patient_info,patient_id=None):
        info = pd.read_csv(patient_info,sep="\t")
        info = info.set_index("Patient ID")
        #Replace ID if given
        if patient_id is None: 
            patient_id = self.id 
        else:
            self.id = patient_id 
        self.sys_bp = info["Systolic BP"][patient_id]
        self.dias_bp = info["Diastolic BP"][patient_id]
        self.heart_rate = info["Heart Rate"][patient_id] 
        self.mean_bp = compute_meanbp(self.sys_bp,self.dias_bp)
        self.amp = self.compute_bpamplitude()

        #Add them to the regions we will process
        self.all.set_reference(self.dias_bp,self.sys_bp,self.mean_bp,self.amp)
        self.inlet.set_reference(self.dias_bp,self.sys_bp,self.mean_bp,self.amp)

    #Add Experiment info 
    def add_exp_info(self,exp, in_flow_s, R, C):
        self.exp = exp
        self.inflow_s = in_flow_s
        self.R = R
        self.C = C
        self.exp_dic = {"Patient ID":self.id,"Within 10%":None,"Within 5%":None,"Exp":self.exp, "Job Directory":self.job_dir,
        "In":self.inflow_s, "R":self.R, "C":self.C} 

    #order: How many points on each side to use for the comparison to consider comparator(n, n+x) to be True.
    def compare_to_exp(self,res_dir,region="inlet",order=5):
        nb_min_max_same = True
        print("Using {} points on each side for comparison...".format(order))
        pressure = res_dir["pressure"] #View (will update original as well)
        #All regions?
        #faces = [x for x in pressure.columns if x != "step"] 
        faces = ["inlet", "All Min", "All Max"]
        #import pdb; pdb.set_trace()

        #Are we interested in regions that have a min or max at one point?
        faces += pressure["All Min Region"].unique().tolist()
        faces += pressure["All Max Region"].unique().tolist()

        #Get local min/max for each region (including all regions)
        for region in faces:
            min_idx = argrelextrema(pressure[region].values, np.less_equal,
                        order=order)[0]
            max_idx = argrelextrema(pressure[region].values, np.greater_equal,
                        order=order)[0]
            pressure[region+'_min'] = pressure.iloc[min_idx][region]
            pressure[region+'_max'] = pressure.iloc[max_idx][region]
            
            if len(min_idx)-1 != len(max_idx):
                print("[WARNING!] (nb_min ({}) - 1) != nb_max ({}), try increasing 'order' - i.e. the number of points used for comparison on each side".format(len(min_idx),len(max_idx)))
                print("Min idx: {}".format(min_idx))
                print("Max idx: {}".format(max_idx))
                nb_min_max_same = False
                print("Breaking loop! Will not set regions...")
                break

            #Set experiment variables and compare to reference values
            #Note that the first min is skipped since diastole if after systole
            if region == "inlet": #inlet 
                self.inlet.set_min(pressure,min_idx[1:]) #Skip 0
                self.inlet.set_max(pressure,max_idx)
            elif region == "All Min": #Minimum of all minima
                pressure[region+'_min Region'] = pressure.iloc[min_idx]["All Min Region"]
                self.all.set_min(pressure,min_idx[1:]) #Skip 0
            elif region == "All Max": #Maximum of all maxima
                pressure[region+'_max Region'] = pressure.iloc[max_idx]["All Max Region"]
                self.all.set_max(pressure,max_idx)

        
        self.nb_cycles = self.inlet.nb_cycles
        #import pdb; pdb.set_trace()
        #self.validate_exp()

        return res_dir, nb_min_max_same

    #Save everything in a dictonary to easily write to a tsv file
    def convert_to_dic(self):
        self.inlet_res_dic = self.merge_exp_res_dics(self.exp_dic,self.inlet.convert_to_dic())
        self.all_res_dic = self.merge_exp_res_dics(self.exp_dic,self.all.convert_to_dic())

    #Merge experiment dictionary with results
    def merge_exp_res_dics(self,exp_dic,res_dic_ref):
        res_dic = exp_dic.copy()
        res_dic.update(res_dic_ref)
        return res_dic

    def save(self,fpath):
        if os.path.exists(fpath):
            print("Overwriting {}...".format(fpath))
        with open(fpath,"w") as writer: 
            writer.write("{}\n".format(self.__str__())) 
        print("Done saving Patient Data to {}...".format(fpath))
   
    def save_tsv(self,fpath):
        if self.inlet_res_dic is None: 
            self.convert_to_dic()
        header = [*self.inlet_res_dic]
        res_df = pd.DataFrame(self.inlet_res_dic,index=[0],columns=header)
        all_res_df = pd.DataFrame(self.all_res_dic,index=[0],columns=header)
        res_df = pd.concat([res_df,all_res_df])
        if os.path.exists(fpath):
            res_df.to_csv(fpath,mode="a",float_format="%.2f",sep="\t",header=False,index=False)
        else: 
            res_df.to_csv(fpath,mode="w",float_format="%.2f",sep="\t",header=True,index=False)
        print("Done saving Patient Data to tsv {}...".format(fpath))
    

    def __str__(self):
        return "Patient ID: {} | Systolic BP: {} | Diastolic BP: {} | Heart Rate: {}\
        \nMean BP: {:.2g} | BP  Amplitude: {:.2g}\
        \nJob: {} | Nb Cycles: {}\
        \n{}\
        \n{}".format(
        self.id, self.sys_bp, self.dias_bp, self.heart_rate, 
        self.mean_bp, self.amp,
        self.job_dir,
        self.nb_cycles, 
        self.inlet,
        self.all)
 
#Region class 
#Compare to reference using the column that corresponds to region 
#Look at the min and the max of that region
#Uses Min_max_res
class Region:
    def __init__(self,region):
        self.region = region 

        self.nb_cycles = 0 
        self.min = Min_max_res(region,"Min")
        self.max = Min_max_res(region,"Max")
        self.mean = Min_max_res(region,"Mean")
        self.amp = Min_max_res(region,"Amp")
        self.w10 = None 
        self.w5 = None

        self.res_dic = None
    
    #Set systolic and diastolic BP
    def set_reference(self,min_ref,max_ref,mean_ref,amp_ref):
        self.min.set_reference(min_ref)
        self.max.set_reference(max_ref)
        self.mean.set_reference(mean_ref)
        self.amp.set_reference(amp_ref) #Amplitude 

    # Check if we are within 10% or 5% of the reference for both max and min
    def validate_exp(self):
        self.w10 = (self.min.w10 and self.max.w10)
        self.w5 = (self.min.w5 and self.max.w5)
        self.w10 = (self.min.w10 and self.max.w10 and self.mean.w10 and self.amp.w10)
        self.w5 = (self.min.w5 and self.max.w5 and self.mean.w5 and self.amp.w5)

    def set_min(self,pressure,idx):
        self.min.set_value(pressure,idx)
        self.update_values(self.min.nb_cycles)
        # self.nb_cycles = self.min.nb_cycles
        # self.validate_exp()

    def set_max(self,pressure,idx):
        self.max.set_value(pressure,idx)
        self.update_values(self.max.nb_cycles)
        # self.nb_cycles = self.max.nb_cycles
        # self.validate_exp()
        
    def convert_to_dic(self): 
        self.res_dic = {"Region": self.region, "Within 10%":self.w10, "Within 5%":self.w5}
        self.res_dic.update(self.min.convert_to_dic())
        self.res_dic.update(self.max.convert_to_dic())
        self.res_dic.update(self.mean.convert_to_dic())
        self.res_dic.update(self.amp.convert_to_dic())
        return self.res_dic

    #Update number of cycles, mean pressure, and amplitude 
    #Calling this after a min is set and after a max is set 
    # in case a max is set before a min
    def update_values(self,nb_cycles):
        self.nb_cycles = nb_cycles
        if self.max.value != [0] and self.min.value != [0]:
            max_arr = np.array(self.max.value)
            min_arr = np.array(self.min.value)
            mean_arr = compute_meanbp(max_arr,min_arr)
            self.mean.set_value(mean_arr,idx=None)
            self.amp.set_value(max_arr-min_arr,idx=None)
        self.validate_exp()

        
    def __str__(self):
        return "\n===Region: {} (Nb Cycles = {} | 10%: {} | 5%: {})\
    \n   Within 10%: {}\
    \n   Within 5%: {}\
    \n{}\
    \n{}\
    \n{}\
    \n{}".format(
    self.region, self.nb_cycles, self.w10, self.w5,
    format_true(self.w10), 
    format_true(self.w5),
    self.min.__str__(), 
    self.max.__str__(), 
    self.mean.__str__(), 
    self.amp.__str__())

#Compares to a reference and computes all the values related to it
# Can be used as min or max since it only deals with a column and a reference (numbers)
# The class itself only uses the type (min/max) to print
class Min_max_res:
    def __init__(self,region,type_val):
        self.region = region
        self.type = type_val # Min/Max 

        #Use set_reference
        self.ref = 0

        #Use set_value for these (might add to constructor)
        self.nb_cycles = 0 
        self.lvalue = 0 #Last min 
        self.lvalue_region = ""
        self.value = [0] #Last min 
        self.value_region = [""]

        #Use validate exp after setting other values
        self.lvalue_abs_diff = 0
        self.lvalue_rel_diff = 0 
        self.w10 = None 
        self.w5 = None

        self.res_dic = None 

    def set_reference(self,ref):
        self.ref = ref 
    
    # Compute difference from reference and 
    # check if we are within 10% or 5% of the reference
    def validate_exp(self):
        self.lvalue_abs_diff = self.lvalue-self.ref
        self.lvalue_rel_diff = 100*(self.lvalue_abs_diff)/self.ref
        self.w10 = bool(np.abs(self.lvalue_rel_diff) < 10) #Np bool acts differently
        self.w5 = bool(np.abs(self.lvalue_rel_diff) < 5)

    #pressure: panda DF or numpy array
    def set_value(self,pressure,idx):
        if "All" in self.region: #Columns are All Min OR All Max
            region = self.region + " " + self.type 
        else: 
            region = self.region
        #IF it's a numpy array...
        if type(pressure).__module__ == np.__name__:
            self.value = pressure.tolist()
        else:
            self.value = pressure.iloc[idx][region].tolist() #Value for all cycle
        self.lvalue = self.value[-1] #Value for the last cycle
        
        #Add the region where that value is  (only really applies for overall min/max)
        if type(pressure).__module__ != np.__name__ and region + " Region" in pressure.columns:
            self.value_region = pressure.iloc[idx][region + " Region"].tolist()
            self.lvalue_region = self.value_region[-1]
        else: 
            self.value_region = [self.region]*len(self.value)
            self.lvalue_region = self.region

        #Update number of cycles based of size of value array
        if self.nb_cycles != 0 and self.nb_cycles != len(self.value):
            warnings.warn("[ISSUE]: {} cycles stored but {} cycles for {}".format(
                    self.nb_cycles,len(self.value),len(self.type)),RuntimeWarning)
        else:
            self.nb_cycles = len(self.value_region)

        self.validate_exp()

    def convert_to_dic(self):
        self.res_dic = {self.type+" w10": self.w10, self.type+" w5": self.w5, 
        self.type+" Region": self.lvalue_region, self.type+" Ref":self.ref, self.type : self.lvalue,
        self.type+" Diff": self.lvalue_abs_diff, self.type+" % Diff": self.lvalue_rel_diff}
        return self.res_dic

    def __str__(self):
        return "-----Type: {}\
        \n       Last {}: {:.1f} (diff: {:.1f} | %-diff: {:.1f}%) [{}]\
        \n       {} Regions: {}\
        \n       {}        : {}\
        \n       Within 10%: {} | Within 5%: {}".format(
        self.type,
        self.type, self.lvalue,  self.lvalue_abs_diff, self.lvalue_rel_diff, self.lvalue_region,
        self.type, self.value_region, 
        self.type, ["{:.1f}".format(x) for x in self.value], 
        format_true(self.w10), 
        format_true(self.w5))