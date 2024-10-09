import argparse
import pandas as pd

import os
import re
import pdb

from patient import Patient 
import matplotlib.pyplot as plt
import numpy as np 

# Parse arguments
def _parse_args():
    parser = argparse.ArgumentParser(prog="range_analyze.py",description='Read txt files from compflow or svFSI output, plot pressure/flow at \
                    different faces, tell how far we are from the patient blood pressure. Optionally perform a mesh convergence \
                    analysis (see --convergence_analysis argument for details). \
                    Regex to parse directory name: (.*)_p(\d)-sim-in(.*)-R(.*)-C(.*)')
    parser.add_argument('--converted_res_dir', type=str, default=".",
                    help="(1) Folder with txt files from avFSI/compflow and vtu/vtp result files\n" \
                        "(OR 2) A path to a list of such folders (see example)")
    parser.add_argument('--out_dir', type=str, default="./analysis_output", help="Output directory \
                    (where plots and output files will be written)")
    parser.add_argument('--patient_info', type=str, default="./patient_info.tsv", help="tsv with \
                    patients BP and HR")
    parser.add_argument('--nb_points_comp', type=int, default=400, help="Number of points before and after \
                    the current one to look for a local min/max (good estimate: nb points per cycle/2 )")
    parser.add_argument('--all_tsv_file', type=str, default="./aorta_all_results.tsv",help="tsv file where all \
                    the results are also added as new rows \
                    (good if you want to compare to other runs)")
    parser.add_argument('-p', '--column_prefix', default="", type=str, help="Prefix to remove from column names  \
                     to match face names")
    parser.add_argument('--show_cycle_diff', action='store_true', help='Will show cycle-to-cycle convergence')
    parser.add_argument('--convergence_analysis', action='store_true', help='Will perform a mesh convergence analysis if ' \
        'a list of folders is provided or other folders with the right naming convention are found (see read_exp_dir())')               
    parser.add_argument('--keep_first_n_points', type=int, help='Number of points to use ' \
        '(useful if a simulation was not finished or ran for longer than necessary)') 
    parser.add_argument('--keep_first_n_cycles', type=int, help='Number of cycles to show ' \
        '(useful if a simulation was not finished or ran for longer than necessary)') 
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

#==========================================
#region       READING FUNCTIONS
#==========================================

#Either accepts a list of directories to process 
#or a single directory (string put in a list with 1 item)
#deduce_list: Will find other meshes for the same experiment 
def read_exp_dir(converted_res_dir,deduce_list=False): 
    mesh_sizes = None
    exp_dir = ""
    if os.path.isfile(converted_res_dir):
        df = pd.read_csv(converted_res_dir,sep="\t")

        # Check number of columns
        if df.shape[1] == 2:
            all_dirs =  df[0].tolist()
            mesh_sizes = df[1].tolist()
        else:
            all_dirs = df[0].tolist()
    else: 
        #TODO: Remove if for other people (only works when results are stored as ours)
        if deduce_list: 
            converted_res_dir = os.path.abspath(converted_res_dir)
            all_dirs = [converted_res_dir]
            #Specific to our jobs 
            # root/exp_dir/exp_dir-converted-results
            root_dir = os.path.dirname(converted_res_dir)
            exp_dir = os.path.basename(root_dir)
            exp_dirs = os.path.join(exp_dir,os.path.basename(converted_res_dir))
            root_dir = os.path.dirname(root_dir)

            mesh_nb = converted_res_dir.split("_p")[0][-2:]
            mesh_sizes_str = ["20","21","22","24"]
            if mesh_nb in mesh_sizes_str:
                mesh_sizes_dic = {"20":2,"21":1.5,"22":1,"24":0.8}
                mesh_sizes = [mesh_sizes_dic[mesh_nb]]
                for mesh_nb_i in mesh_sizes_str: 
                    if mesh_nb_i != mesh_nb:
                        exp_dirs_i = exp_dirs.replace(f"_{mesh_nb}_",f"_{mesh_nb_i}_")
                        mesh_dir = os.path.join(root_dir,exp_dirs_i)
                        if os.path.exists(mesh_dir):
                            all_dirs.append(mesh_dir)
                            mesh_sizes.append(mesh_sizes_dic[mesh_nb_i])
        else:
            all_dirs = [converted_res_dir]
    return all_dirs, mesh_sizes, exp_dir

#Read txt files as pandas DFs and store each DF as a dictionnary entry 
def read_txt_results(converted_res_dir,prefix="all_results",column_prefix="",keep_first_n_points=None):
    res_dic_paths = {"flow": prefix+"-flows.txt", "pressure": prefix+"-pressures.txt",
            "avg": prefix+"-averages.txt", "avg_mm": prefix+"-averages-from_cm-to-mmHg-L_per_min.txt"}
    res_dic_paths_svfsi = {"flow": "B_FS_Velocity_flux.txt", "pressure": "B_FS_Pressure_average.txt",
            "avg": None, "avg_mm": None}
    if "fl" in converted_res_dir:
        res_dic_paths_svfsi = {"flow": "B_NS_Velocity_flux.txt", "pressure": "B_NS_Pressure_average.txt",
                "avg": None, "avg_mm": None}

    if not os.path.exists(os.path.join(converted_res_dir,res_dic_paths["pressure"])):
        res_dic_paths = res_dic_paths_svfsi
        svfsi = True 
        sep = " "
    else:
        svfsi = False
        sep = "\t"

    res_dic = {} #All files 
    for k,v in res_dic_paths.items():
        if v is not None:
            print("Reading: {} from {}".format(v,converted_res_dir))
            clean_csv(os.path.join(converted_res_dir,v))
            res_dic[k] = pd.read_csv(os.path.join(converted_res_dir,v),sep=sep, skipinitialspace=True)
            res_dic[k] = res_dic[k].dropna(axis=1, how='all')

    #Get flow rate (first row is area)
    if svfsi: 
        res_dic["pressure"] = restrict_columns(res_dic["pressure"],column_prefix=column_prefix)
        res_dic["flow"] = restrict_columns(res_dic["flow"],column_prefix=column_prefix)
        res_dic["flow"] = res_dic["flow"].drop(res_dic["flow"].index[0]).reset_index(drop=True)
        res_dic["pressure"].reset_index(drop=True, inplace=True)

    print(f"{len(res_dic['pressure'])} points found")
    #Keep first n points
    if keep_first_n_points is not None:
        res_dic["pressure"] = res_dic["pressure"].iloc[:keep_first_n_points]
        res_dic["flow"] = res_dic["flow"].iloc[:keep_first_n_points]
    return res_dic

def restrict_columns(df,column_prefix=""):
    if column_prefix != "":
        df_filt = df.filter(like=column_prefix)
        if not df_filt.empty:   
            df = df_filt.rename(columns=lambda x: x.replace(column_prefix, ''))
        else:
            print("Prefix {} not found".format(column_prefix))
    return df 

def clean_csv(input_file):
    # Read the file and remove trailing spaces
    with open(input_file, 'r') as file:
        lines = file.readlines()
        lines = [line.rstrip() for line in lines]

    # Write the cleaned lines back to the file
    with open(input_file, 'w') as file:
        for line in lines:
            file.write(line + '\n')

#endregion    END - READING 
#==========================================

#==========================================
#region       PROCESSING HELPER FUNCTIONS
#==========================================

#10_p2-sim-in0.85-R2-C0.7
#Parse jobDir to retrieve experiment name, inflow, R, and C
def parse_job_dir(jobDir,args_in=None):
    exp=""
    patient="2" #Default
    in_flow=0
    R=0
    C=0
    viscosity = "constant"

    #Parse 
    print("Parsing {}".format(jobDir))
    job_reg = re.compile(r'(.*)_p(\d).*-sim-in(.*)-R(.*)-C(.*)')
    parsed_reg = re.search(job_reg,jobDir)
    if parsed_reg is not None: 
        exp, patient, in_flow, R, C = parsed_reg.groups()
        if "_" in C:
            C, viscosity = C.split("_")

    #Overwrite values if they are provided from command-line
    if args_in is not None:
        if args_in.exp is not None: 
            exp = args_in.exp
        if args_in.patient is not None: 
            patient = args_in.patient
        if args_in.in_flow is not None: 
            in_flow = args_in.in_flow
        if args_in.R is not None: 
            R = args_in.R
        if args_in.C is not None: 
            C = args_in.C
    return exp, int(float(patient)), float(in_flow), float(R), float(C), viscosity

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

#endregion    END - PROCESSING HELPER 
#==========================================

#==========================================
#region       PLOTTING FUNCTIONS
#==========================================
#Plot pressure or flow for a given region
def plot_res(res_dir,data_type="pressure",region="top", out_dir="./analysis_output",prefix="",suffix="",units="",plot_patient=False):
    #Plot values for the given region (T/M/B or a single face)
    faces = {}
    faces["top"] = ["inlet","b-trunk", "carotid", "subclavian"]
    faces["middle"] = ["inlet","celiac","superior-mesenteric","right-renal","left-renal"]
    faces["bottom"] = ["inlet","right-ext-illiac","right-int-illiac","left-int-illiac","left-ext-illiac"]
    # faces["top"] = ["inlet","outlet","wall"]
    # faces["middle"] = ["inlet","outlet","wall"]
    # faces["bottom"] = ["inlet","outlet","wall"]
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
    # axes = res_dir[data_type][faces_to_plot].plot.line(style='.-',title=title + " ({})".format(units),grid = 'on')
    axes = res_dir[data_type][faces_to_plot].plot.line(style='.-',title=title + "\n({})".format(units),grid = 'on')
    #Plot parameters
    axes.minorticks_on()
    axes.grid(linestyle='--', linewidth=0.3)
    #axes.tick_params(axis='x',which='minor',bottom='off')
    # axes.set_ylim(-5, 130)
    output_path = os.path.join(out_dir,"{}.png".format(title))
    os.makedirs(out_dir,exist_ok=True)
    # axes.get_figure().tight_layout(rect=[0.1, 0.1, 0.9, 0.9])
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

#endregion    END - PLOTTING
#==========================================

#==========================================
#region       HIGH-LEVEL FUNCTION/CLASSES
#==========================================

#Attributes needed to parse dir 
class Args_in:
    def __init__(self,exp=None,patient=None,in_flow=None,R=None,C=None):
        self.exp = exp
        self.patient = patient 
        self.in_flow = in_flow
        self.R = R
        self.C = C

#show_diff = True #Print difference between current and previous value for min/max/amp/mean...
#plot_patient = True #Plot patient's systolic and diastolic pressure
#This function processes a single experiment directory
def range_analyze(args,converted_res_dir,out_dir,all_tsv_file,show_diff=True,plot_patient=True,data_types=["pressure", "flow"],regions=["inlet","top", "middle", "bottom"],save=True,keep_first_n_points=None,keep_first_n_cycles=None):
    #Fixed 
    suffix=""
    #Use the name of the converted-dir directory as prefix for the output files
    job_dir_prefix = False 
    job_dir_prefix = True #REMMMMOOOOOVVVEEEE for directories with right naming conventions
    
    #Arg dependent 
    column_prefix =  args.column_prefix 
    patient_info = args.patient_info
    order=args.nb_points_comp #Nb points before/after for local min/max
    args_in = Args_in(args.exp,args.patient,args.in_flow,args.R,args.C)

    jobDir = os.path.basename(os.path.realpath(converted_res_dir)).replace("-converted-results","")
    exp, patient, in_flow, R, C, viscosity = parse_job_dir(jobDir,args_in)

    if job_dir_prefix:
        prefix=jobDir + "_"
    else: 
        if viscosity != "constant":
            prefix="{}_p{}-sim-in{:g}-R{:g}-C{:g}_{}_".format(exp,patient,in_flow,R,C,viscosity)
        else: 
            prefix="{}_p{}-sim-in{:g}-R{:g}-C{:g}".format(exp,patient,in_flow,R,C)

    prefix += column_prefix 

    print("\n-----------------------------\n")

    #Create, read patient values from tsv file, add experiment info
    px = Patient(patient,jobDir)
    px.load_from_file(patient_info,patient)
    px.add_exp_info(exp, in_flow, R, C, viscosity)
    #print(px)
    print("\n-----------------------------\n")

    #Read results and add useful info
    res_dir = read_txt_results(converted_res_dir,column_prefix=column_prefix,keep_first_n_points=keep_first_n_points)
    res_dir = convert_flow_press(res_dir) #Call this before adding string columns
    res_dir = add_stats_col(res_dir)
    res_dir = add_patient_bp(res_dir,px)
    
    #May fail if the right order is not used for comparison
    #Putting in a try/except so the user can use the graphs to debug    
    try: 
        #Compare experiment results to patient BP
        res_dir, nb_min_max_same = px.compare_to_exp(res_dir,order=order,keep_first_n_cycles=keep_first_n_cycles)
        #Could change it into a while loop, will try it twice for now 
        if not nb_min_max_same:
            print("--Doubling order")
            order = order*2
            res_dir, nb_min_max_same = px.compare_to_exp(res_dir,order=order,keep_first_n_cycles=keep_first_n_cycles)

        # suffix += "_d{:.1f}_s{:.1f}".format(px.all.min.lvalue_rel_diff,px.all.max.lvalue_rel_diff)
        suffix += "_d{:.1f}_s{:.1f}".format(px.inlet.min.lvalue_rel_diff,px.inlet.max.lvalue_rel_diff)

        print("Columns/Available faces: {}".format(res_dir[data_types[0]].columns))
        print("\n-----------------------------\n")

        if save:
            #Plot flow/pressure for multiple regions
            plot_all_regions(res_dir,out_dir,data_types,regions,prefix,suffix,plot_patient)
        print("\n-----------------------------\n")

        print(px.__str__(show_diff=show_diff))
        print("\n-----------------------------\n")

        if save:
            px.save(os.path.join(out_dir,prefix+suffix+"_summary.txt"))
        px.save_tsv(all_tsv_file)
    except:
        if save:
            #Plot flow/pressure for multiple regions
            plot_all_regions(res_dir,out_dir,data_types,regions,prefix,suffix,plot_patient)
        print("\n-----------------------------\n") 
        raise

#endregion    END - HIGH-LEVEL FUNCTION/CLASSES
#==========================================

#==========================================
#region       MESH CONVERGENCE ANALYSIS FUNCTIONS
#==========================================

#Example: 
# Experiment	Max edge size (mm)
# .../2_115_22_p1-sim-...-converted-results	1
# .../2_115_20_p1-sim-...-converted-results	2
def write_convergence_list(fname, all_dirs, mesh_sizes,out_dir="."):
    if mesh_sizes is None: 
        print("[WRITE_CONVERGENCE_LIST] No mesh size provided...")
        return 
    file_path = os.path.join(out_dir,fname)
    if len(all_dirs) != len(mesh_sizes):
        raise ValueError("Both lists must have the same length.")
    
    os.makedirs(out_dir,exist_ok=True)
    print(f"Saving convergence list to {os.path.abspath(file_path)}")
    with open(file_path,'w', newline='') as file:
        file.writelines(f"Experiment\tMax edge size (mm)\n")
        file.writelines(f"{item1}\t{item2}\n" for item1, item2 in zip(all_dirs, mesh_sizes))

#Receives list or numpy array and compute difference between current and previous value
def compute_diff(value):
    value_arr = np.array(value)
    _value_abs_diff = value_arr[1:] - value_arr[0:-1]
    value_abs_diff = [None] + _value_abs_diff.tolist()
    value_rel_diff = [None] + (100*_value_abs_diff/value_arr[0:-1]).tolist()

    return value_abs_diff, value_rel_diff

#Do convergence analysis and plot it 
def plot_convergence(all_tsv_file,mesh_sizes,region="inlet",x_axis="Max Edge size (mm)",columns_to_plot  = ['Min', 'Max', 'Mean', 'Amp']): 
    title_str = os.path.basename(all_tsv_file).split("_conve")[0]
    fname = all_tsv_file.replace("_convergence_results.tsv",f"_convergence_graph_{region.lower()}.png")
    print(f"Reading: {all_tsv_file}")
    df = pd.read_csv(all_tsv_file,sep="\t")
    df_inlet = df[df['Region'] == region]
    df_inlet[x_axis] = mesh_sizes
    
    if x_axis == "Max Edge size (mm)":
        ascending = False 
    else: 
        ascending = True 

    # import pdb; pdb.set_trace()
    df_inlet = df_inlet.sort_values(by=x_axis,ascending=ascending)

    plt.figure(figsize=(8, 5)) #figsize=(12, 8)
    # Plot each column against x_axis
    for column in columns_to_plot:
        abs_diff, rel_diff = compute_diff(df_inlet[column].values.tolist())
        print(f"==={column}===")
        print(f"[Diff {column}] {' '.join(f'{v:.3g}' if v is not None else 'None' for v in abs_diff)}")
        print(f"[% - Diff {column}] {','.join(f'{v:.3g}' if v is not None else 'None' for v in rel_diff)}")
        plt.plot(df_inlet[x_axis], df_inlet[column], marker='o', label=column)

        # Add text annotations for each point
        for rel_diff_i, exp, value in zip(rel_diff,df_inlet[x_axis], df_inlet[column]):
            # plt.text(exp, value, f'{value}', fontsize=9, ha='right')
            if rel_diff_i is not None:
                plt.text(exp, value + 0.5, f'{value:.3g} ({rel_diff_i:.3g}%)', fontsize=9, ha='center', va='bottom', 
                    bbox=dict(facecolor='white', alpha=0.6, edgecolor='none'))
            else: 
                plt.text(exp, value + 0.5, f'{value:.3g} ()', fontsize=9, ha='center', va='bottom', 
                    bbox=dict(facecolor='white', alpha=0.6, edgecolor='none'))
    plt.xlabel(x_axis)
    plt.ylabel('Pressure (mmHg)')
    plt.title(f"{title_str}\nMesh convergence analysis, rel. diff. with previous value shown in ()", pad=30)
    if x_axis == "Max Edge size (mm)":
        plt.gca().invert_xaxis()
    # plt.legend()
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=4)
    plt.grid(True)

    # Save the plot to a file
    plt.tight_layout()
    plt.savefig(fname, format='png')

#endregion    END - MESH CONVERGENCE ANALYSIS
#==========================================

if __name__ == '__main__':
    args = _parse_args()
    data_types = ["pressure", "flow"]
    regions=["inlet","top", "middle", "bottom"]
    deduce_list = True #Will find the other meshes using the current one (only valid for naming conventions like ours)
    save = True #Save results (will be turned off for other meshes for convergence analysis)
    plot_patient = True
    show_diff = args.show_cycle_diff #Show cycle to cycle difference 
    convergence_analysis_bool = args.convergence_analysis

    out_dir = args.out_dir
    all_tsv_file = args.all_tsv_file
    converted_res_dirs,mesh_sizes,exp_dir = read_exp_dir(args.converted_res_dir,deduce_list)

    #Find similar experiment with different mesh sizes
    if deduce_list and mesh_sizes is not None:
        convergence_list_name = f"{exp_dir}_convergence_list.txt"
        write_convergence_list(convergence_list_name, converted_res_dirs, mesh_sizes,out_dir=out_dir)
    else: 
        convergence_analysis_bool = False

    #Only do convergence analysis when final mesh is provided 
    if convergence_analysis_bool and ("_22_" not in converted_res_dirs[0] or "2_" != exp_dir[:2]):
        converted_res_dirs = converted_res_dirs[:1]
        convergence_analysis_bool = False
    
    #Will store all the results here 
    if convergence_analysis_bool: 
        all_tsv_file = os.path.join(out_dir,f"{exp_dir}_convergence_results.tsv")

    for i,converted_res_dir in enumerate(converted_res_dirs): 
        range_analyze(args,converted_res_dir,out_dir,all_tsv_file,show_diff=show_diff,plot_patient=plot_patient,
            data_types=data_types,regions=regions,save=save,
            keep_first_n_points=args.keep_first_n_points,keep_first_n_cycles=args.keep_first_n_cycles)
        save=False
        
    #Start the convergence analysis
    if convergence_analysis_bool: 
        plot_convergence(all_tsv_file,mesh_sizes)
