import argparse
import pandas as pd

import os
import re
import pdb

from patient import Patient 

# Parse arguments
def _parse_args():
    parser = argparse.ArgumentParser(prog="range_analyze.py",description='Read txt files from compflow output, plot pressure/flow at \
                    different faces, tell how far we are from the patient blood pressure. \
                    Regex to parse directory name: (.*)_p(\d)-sim-in(.*)-R(.*)-C(.*)')
    parser.add_argument('--converted_res_dir', type=str, default=".",
                    help="Folder with txt files from compflow and vtu/vtp result files")
    parser.add_argument('--out_dir', type=str, default="./analysis_output", help="Output directory \
                    (where plots and output files will be written)")
    parser.add_argument('--patient_info', type=str, default="./patient_info.tsv", help="tsv with \
                    patients BP and HR")
    parser.add_argument('--nb_points_comp', type=int, default=4, help="Number of points before and after \
                    the current one to look for a local min/max")
    parser.add_argument('--all_tsv_file', type=str, default="./aorta_all_results.tsv",help="tsv file where all \
                    the results are also added as new rows \
                    (good if you want to compare to other runs)")
    parser.add_argument('--exp', type=str, help="Experiment name \
                    (will try to extract from folder name if not provided)")
    parser.add_argument('--patient', type=int, help="Patient ID \
                    (will try to extract from folder name if not provided)")
    parser.add_argument('--in_flow', type=float, help="Inlet flow scaling factor \
                    (will try to extract from folder name if not provided)")
    parser.add_argument('-R', type=float, help="R scaling factor \
                    (will try to extract from folder name if not provided)")
    parser.add_argument('-C', type=float, help="C scaling factor \
                    (will try to extract from folder name if not provided)")

    return parser.parse_args()

#Read txt files as pandas DFs and store each DF as a dictionnary entry 
def read_txt_results(converted_res_dir,prefix="all_results"):
    res_dic_paths = {"flow": prefix+"-flows.txt", "pressure": prefix+"-pressures.txt",
            "avg": prefix+"-averages.txt", "avg_mm": prefix+"-averages-from_cm-to-mmHg-L_per_min.txt"}
    
    res_dic = {} #All files 
    for k,v in res_dic_paths.items():
        print("Reading: {} from {}".format(v,converted_res_dir))
        res_dic[k] = pd.read_csv(os.path.join(converted_res_dir,v),sep="\t")
        res_dic[k] = res_dic[k].dropna(axis=1, how='all')

    return res_dic

#Plot pressure or flow for a given region
def plot_res(res_dir,data_type="pressure",region="top", out_dir="./analysis_output",prefix="",suffix="",units="",plot_patient=False):
    #Plot values for the given region (T/M/B or a single face)
    faces = {}
    faces["top"] = ["inlet","b-trunk", "carotid", "subclavian"]
    faces["middle"] = ["inlet","celiac","superior-mesenteric","right-renal","left-renal"]
    faces["bottom"] = ["inlet","right-ext-illiac","right-int-illiac","left-int-illiac","left-ext-illiac"]
    if region in ["top", "middle", "bottom"]:
        faces_to_plot = faces[region]
    elif region == "all":
        faces_to_plot = [x for x in res_dir[data_type].columns if x != "step"]
    else:
        faces_to_plot = [region] 

    if plot_patient and data_type == "pressure":
        faces_to_plot += ["systolic BP", "diastolic BP"]

    title = "{}{}_{}{}".format(prefix,data_type,region,suffix)
    print("Plotting/Saving {} for faces: {}".format(data_type.upper(),faces_to_plot))
    axes = res_dir[data_type][faces_to_plot].plot.line(style='.-',title=title + " ({})".format(units),grid = 'on')
    #Plot parameters
    axes.minorticks_on()
    axes.grid(linestyle='--', linewidth=0.3)
    #axes.tick_params(axis='x',which='minor',bottom='off')

    output_path = os.path.join(out_dir,"{}.png".format(title))
    os.makedirs(out_dir,exist_ok=True)
    fig = axes.get_figure().savefig(output_path)
    print("Saved to {}\n".format(output_path))


#Make a plot for all the regions/types given 
def plot_all_regions(res_dir,out_dir=".",data_types=["pressure","flow"],regions="top",prefix="",suffix="",plot_patient=None):
    #Make plots
    for data_type in data_types:
        if data_type == "pressure":
            units = "mmHg"
        elif data_type == "flow":
            units = "L/min"
        for region in regions:
            plot_res(res_dir,data_type,region,out_dir=out_dir,prefix=prefix,suffix=suffix,units=units,plot_patient=plot_patient)

#10_p2-sim-in0.85-R2-C0.7
#Parse jobDir to retrieve experiment name, inflow, R, and C
def parse_job_dir(jobDir,args=None):
    exp=""
    patient="2" #Default
    in_flow=0
    R=0
    C=0

    #Parse 
    print("Parsing {}".format(jobDir))
    job_reg = re.compile(r'(.*)_p(\d).*-sim-in(.*)-R(.*)-C(.*)')
    parsed_reg = re.search(job_reg,jobDir)
    if parsed_reg is not None: 
        exp, patient, in_flow, R, C = parsed_reg.groups()
    
    #Overwrite values if they are provided from command-line
    if args is not None:
        if args.exp is not None: 
            exp = args.exp
        if args.patient is not None: 
            patient = args.patient
        if args.in_flow is not None: 
            in_flow = args.in_flow
        if args.R is not None: 
            R = args.R
        if args.C is not None: 
            C = args.C
    return exp, int(float(patient)), float(in_flow), float(R), float(C)

#Add patient Bloop Pressure to DF
def add_patient_bp(res_dir,myPatient):
    res_dir["pressure"]["systolic BP"] = myPatient.sys_bp
    res_dir["pressure"]["diastolic BP"] = myPatient.dias_bp
    return res_dir 

#Convert flow and pressure to L/min and mmHg respectively
def convert_flow_press(res_dir):
    #1 mmHg to dyn/mm^2	= 133.2
    #1 L to 1 mm^3 = 10^6
    res_dir["pressure_orig"] = res_dir["pressure"]
    res_dir["flow_orig"] = res_dir["flow"]

    res_dir["pressure"] = res_dir["pressure"]/133.224
    res_dir["flow"] = (res_dir["flow"])/10**6 
    # res_dir["flow"] = (res_dir["flow"]*60)/10**6 
    return res_dir  

#Add statistic per row as new columns (min/max/avg and region corresponding to mix/max)
def add_stats_col(res_dir):
    faces = [x for x in res_dir["pressure"].columns if x != "step"]
    res_dir["pressure"]["All Min"] = res_dir["pressure"][faces].min(axis="columns")
    res_dir["pressure"]["All Avg"] = res_dir["pressure"][faces].mean(axis="columns")
    res_dir["pressure"]["All Max"] = res_dir["pressure"][faces].max(axis="columns")
    #Store region associated to min/max
    res_dir["pressure"]["All Min Region"] = res_dir["pressure"][faces].idxmin(axis="columns")
    res_dir["pressure"]["All Max Region"] = res_dir["pressure"][faces].idxmax(axis="columns")
    return res_dir


if __name__ == '__main__':
    args = _parse_args()

    #Fixed 
    data_types = ["pressure", "flow"]
    regions = ["inlet","top", "middle", "bottom"]
    suffix=""
    plot_patient = True #Plot patient's systolic and diastolic pressure
    #Use the name of the converted-dir directory as prefix for the output files
    job_dir_prefix = False 

    #Arg dependent 
    order=args.nb_points_comp #Nub points before/after for local min/max
    all_tsv_file = args.all_tsv_file #File 
    jobDir = os.path.basename(os.path.realpath(args.converted_res_dir)).replace("-converted-results","")
    exp, patient, in_flow, R, C = parse_job_dir(jobDir,args)
    if job_dir_prefix:
        prefix=jobDir + "_"
    else: 
        prefix="{}_p{}-sim-in{:g}-R{:g}-C{:g}_".format(exp,patient,in_flow,R,C)
    print("\n-----------------------------\n")

    #Create, read patient values from tsv file, add experiment info
    px = Patient(patient,jobDir)
    px.load_from_file(args.patient_info,patient)
    px.add_exp_info(exp, in_flow, R, C)
    #print(px)
    print("\n-----------------------------\n")

    #Read results and add useful info
    res_dir = read_txt_results(args.converted_res_dir)
    res_dir = convert_flow_press(res_dir) #Call this before adding string columns
    res_dir = add_stats_col(res_dir)
    res_dir = add_patient_bp(res_dir,px)
    
    #Compare experiment results to patient BP
    res_dir, nb_min_max_same = px.compare_to_exp(res_dir,order=order)
    #Could change it into a while loop, will try it twice for now 
    if not nb_min_max_same:
        print("--Doubling order")
        order = order*2
        res_dir, nb_min_max_same = px.compare_to_exp(res_dir,order=order)

    # suffix += "_d{:.1f}_s{:.1f}".format(px.all.min.lvalue_rel_diff,px.all.max.lvalue_rel_diff)
    suffix += "_d{:.1f}_s{:.1f}".format(px.inlet.min.lvalue_rel_diff,px.inlet.max.lvalue_rel_diff)

    print("Columns/Available faces: {}".format(res_dir[data_types[0]].columns))
    print("\n-----------------------------\n")

    #Plot flow/pressure for multiple regions
    plot_all_regions(res_dir,args.out_dir,data_types,regions,prefix,suffix,plot_patient)
    print("\n-----------------------------\n")

    print(px)
    print("\n-----------------------------\n")

    px.save(os.path.join(args.out_dir,prefix+suffix+"_summary.txt"))
    px.save_tsv(all_tsv_file)
    
