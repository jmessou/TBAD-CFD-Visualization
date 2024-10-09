import argparse
from argparse import RawTextHelpFormatter
import warnings
import os

import pandas as pd
import xml.etree.ElementTree as ET

# Parse arguments
def _parse_args():
    #~/bin/pvpython convert_obj_to_vtp.py
    parser = argparse.ArgumentParser(prog="convert_txt_to_mdl.py",
        description="Convert txt file with ID/face names to mdl file "
            "that is supported by SimVascular", formatter_class=RawTextHelpFormatter)
    parser.add_argument("-i",'--input_file_face_id', required=True,type=str, help="Input file (.txt) "
            "- output of get_face_names.py")
    parser.add_argument('--save_path', type=str, default=None, help="Output file (.mdl) path"
            "Defaults to the input file name without the 'face_id.txt'")

    return parser.parse_args()

#==========================================
#region          HIGH-LEVEL FUNCTIONS
#==========================================
def convert_df_to_mdl(face_df,save_path,template_path="template.mdl"):
    print("Converting face ids/face names to Simvascular .mdl file using: \n"
        "---template: {}\n"
        "---output_path: {}".format(template_path,save_path))

    # Get face info from face_df
    faces_list = face_df.apply(lambda row: get_mdl_face_dic(row.name,row["face_names"],row["face_type"]), axis=1).tolist()
    
    # Read Simvascular mdl template 
    tree = ET.parse(template_path)
    root = tree.getroot()

    # Remove faces from template
    faces_section = root.find('.//faces')
    faces_section.clear() 

    #Add faces found in face_df 
    for face_info in faces_list:
        ET.SubElement(faces_section,'face', face_info)

    #Fix indentation for added elements and save
    ET.indent(tree, space = '\t', level=0)
    tree.write(save_path, encoding="UTF-8", xml_declaration=True)


#See simvacular template 
def get_mdl_face_dic(face_id,face_name,face_type="cap"):
    mdl_face_dic = {'id': str(face_id), 'name': face_name, 'type': face_type, 'visible': 'true', 'opacity': '1', 'color1': '1', 'color2': '1', 'color3': '1'}
    return mdl_face_dic 
#def 
#endregion         END - HIGH-LEVEL FUNCTIONS
#==========================================

if __name__ == '__main__':
    args = _parse_args()

    if args.save_path is None: 
        args.save_path = args.input_file_face_id.replace("_face_id.txt",".mdl")
    
    face_df = pd.read_csv(args.input_file_face_id, index_col=0, sep=" ")
    template_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),"template.mdl")

    convert_df_to_mdl(face_df,save_path=args.save_path,template_path=template_path)

