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

# Parse arguments
def _parse_args():
    parser = argparse.ArgumentParser(prog="get_rcr.py",description='Write RCR boundaries based on outlet flow fractions \
                    and R/C scaling factors')
    parser.add_argument('--flow_fractions_p', type=str, default="./outlet_flow_fractions_bm.txt",
                    help="Path to outlet flow fractions")
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
    output_path = os.path.join(args.out_dir,"R{}-C{}.txt".format(args.R,args.C))
    print("Writing to: {}".format(output_path))
    flow_fractions[columns_to_write].to_csv(output_path,index=False,float_format='%.6g',sep=" ")
