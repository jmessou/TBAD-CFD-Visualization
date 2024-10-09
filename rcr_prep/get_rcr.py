import argparse
import pandas as pd
import os

import pdb

#Read file with outlet flow fractions and store it in a tsv file
def read_flow_fractions():
    outlet_order = pd.read_csv(args.outlet_order)

# Add an output order for the output file
def add_output_order(flow_fractions,outlet_order_p):
    if outlet_order_p is not None and os.path.isfile(outlet_order_p): 
        outlet_order = pd.read_csv(outlet_order_p)
        outlet_order["Output order"] = outlet_order.index.values.astype(int)
        outlet_order = outlet_order.set_index("Outlet")

        flow_fractions = flow_fractions.set_index("Outlet")
        flow_fractions["Output order"] = outlet_order["Output order"]
        flow_fractions = flow_fractions.reset_index()
    else:
        flow_fractions["Output order"] = flow_fractions.index.values.astype(int)
    
    return flow_fractions

def comp_rcr(flow_fractions, R=1, C=1):
    ctotal_orig = 7.1 #Baumler
    rtotal_orig = 0.084 #Baumler
    k = 0.9
    ctotal_scaled = C*ctotal_orig
    rtotal_scaled =  R*rtotal_orig
    flow_fractions["Rp"] = (1-k)*(rtotal_scaled/flow_fractions["Fraction"])
    flow_fractions["Rd"] = (k)*(rtotal_scaled/flow_fractions["Fraction"])
    flow_fractions["C"] = ctotal_scaled*flow_fractions["Fraction"]

    return flow_fractions 

#Add RCR Boundary Conditions to svFSI file
#Needs file to already have matching Neumann BCs as RCR
#Written in a simple way, without full parsing the file
def add_rcr_to_svfsi(svfsi_file,flow_fractions):
    if type(flow_fractions) == str: 
        flow_fractions = pd.read_csv(flow_fractions)

    print("Replacing BCs in {}".format(svfsi_file))
    output_svfsi_file_path = add_suffix(svfsi_file,"_rcr")
    lines, empty = read_file(svfsi_file)

    if not empty:
        new_lines = replace_rcr(lines,flow_fractions)
        write_file(new_lines, output_svfsi_file_path)   
    # import pdb; pdb.set_trace()

#Go through the lines and replace resistance by RCR if necessary and 
#add RCR values 
def replace_rcr(lines,flow_fractions):
    new_lines = []
    skip = False 
    distal_pressure = 0 
    init_pressure = 0
    for i, line in enumerate(lines): 
        if not skip:
            new_lines.append(line)
        skip = False 

        if "Resistance" in line or ("RCR" in line and "dependence" in line): 
            skip = True 
            face = lines[i-2].split(":")[1].split(" ")[1].replace("lumen_","")
            face_row = flow_fractions[flow_fractions['Outlet names'].str.contains(face, regex=False)]
            if len(face_row) > 1:
                raise Exception("Too many rows ({}) matching face ({}): {}".format(len(face_row),face,face_row))
            elif len(face_row) == 0:
                raise Exception("No match for face ({}).".format(face))

            #If it was a resistance, make it a RCR
            new_lines[-1] = lines[i].replace("Resistance","RCR") 
            
            #Get sample line with right spaces and no content
            start_index = lines[i].find("Time")
            empty_line = lines[i][:start_index]
            rp,c,rd = (face_row["Rp"].values[0],face_row["C"].values[0],face_row["Rd"].values[0])
            
            new_lines.append("{}RCR values: ({:.6g},{:.6g},{:.6g})".format(empty_line,rp,c,rd))
            if "Resistance" in line:
                new_lines.append("{}Distal pressure: {}".format(empty_line,distal_pressure))
                new_lines.append("{}Initial pressure: {}".format(empty_line,init_pressure))

    return new_lines 


def read_file(file_path):
    lines = []
    empty = False
    with open(file_path) as file:
        lines = [line.rstrip() for line in file]
    if len(lines) == 0:
        empty = True
        print("File {} not found or empty.".format(file_path))
    return lines, empty

def write_file(lines, output_file_path):
    print("Writing: {}".format(output_file_path))
    with open(output_file_path, 'w') as file:
        for line in lines:
            file.write(line + '\n')

#Add suffix to a path
def add_suffix(path,suffix):
    path_split = path.rsplit(".",1)
    path = path_split[0] + suffix + "." + path_split[1]
    
    return path

# Parse arguments
def _parse_args():
    parser = argparse.ArgumentParser(prog="get_rcr.py",description='Write RCR boundaries based on outlet flow fractions \
                    and R/C scaling factors. Replace RCR BCs in svFSI.inp file if provided (see add_rcr_to_svfsi() function for expected format).')
    parser.add_argument('--flow_fractions_p', type=str, default="./outlet_flow_fractions_bm.txt",
                    help="Path to outlet flow fractions")
    parser.add_argument('-s','--svfsi_file', type=str,
                    help="Path to svFSI file where RCR Boundary conditions will be written.")
    parser.add_argument('--outlet_names', type=str, default="outlet_names.txt",
                    help="File with potential face names (needed if svfsi_file provided).")
    parser.add_argument('--out_dir', type=str, default=".", help="Output directory \
                    (where RCR boundary conditions will be written)")
    parser.add_argument('--outlet_order_p', type=str, default=None, help="Path to file \
                    with the order in which the outlets will be written to the output file")
    parser.add_argument('-R', type=float, default="1", help="R scaling factor")
    parser.add_argument('-C', type=float, default="1", help="C scaling factor")

    return parser.parse_args()

if __name__ == '__main__':
    args = _parse_args()

    columns_to_write = ["Outlet","Rp","C","Rd"]
    flow_fractions = pd.read_csv(args.flow_fractions_p,sep="\t")
    flow_fractions = add_output_order(flow_fractions,args.outlet_order_p)
    flow_fractions = comp_rcr(flow_fractions, args.R, args.C) #Scale R and C
    flow_fractions = flow_fractions.set_index("Output order")
    flow_fractions = flow_fractions.sort_index() #Sort in the output order we want
    print(flow_fractions)
    if os.path.isdir(args.out_dir):
        output_path = os.path.join(args.out_dir,"R{}-C{}.txt".format(args.R,args.C))
        print("Writing to: {}".format(output_path))
        flow_fractions[columns_to_write].to_csv(output_path,index=False,float_format='%.6g',sep=" ")

    if args.svfsi_file: 
        outlet_names_df = pd.read_csv(args.outlet_names, comment='#',header=None)
        outlet_names_df['Outlet names'] = outlet_names_df.apply(lambda row: ','.join(map(str, row)), axis=1)
        outlet_names_df = outlet_names_df.rename(columns={0: 'Outlet'})
        columns_to_write += ['Outlet names']
        flow_fractions = pd.merge(flow_fractions,outlet_names_df[['Outlet',"Outlet names"]],on="Outlet",how="inner")
        add_rcr_to_svfsi(args.svfsi_file,flow_fractions[columns_to_write])