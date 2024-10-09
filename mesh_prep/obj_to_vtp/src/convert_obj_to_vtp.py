import argparse
from argparse import RawTextHelpFormatter
import paraview
import warnings
import os
import glob
paraview.compatibility.major = 5
paraview.compatibility.minor = 11

from paraview.simple import *
from paraview_utils import get_array_names, load_vtp_and_get_arrays, print_header

from get_face_names import get_face_names
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()


# Parse arguments
def _parse_args():
    #~/bin/pvpython convert_obj_to_vtp.py
    parser = argparse.ArgumentParser(prog="convert_obj_to_vtp.py",
        description="(1) Convert .obj file to .vtp file using Paraview.\n"
            "(2) Add 1 to 'GroupIds' and replace it with 'ModelFaceID' to support face identication in Simvascular.\n"
            "(3) Automatically attribute face names and generate separate .mdl file (See SimVascular).\n"
            "Run as PATH_TO_PARAVIEW/bin/pvpython convert_obj_to_vtp.py", formatter_class=RawTextHelpFormatter)
    parser.add_argument("-i",'--input_file_obj', required=True,type=str, help="Input .obj file to be converted\n"
            "Convert all available .obj files if directory is provided.")
    parser.add_argument('--model_type', type=str, default="solid",
            help="Model type (fluid or solid), 13 faces for fluid, 14 faces for solid.")
    parser.add_argument("-o",'--out_dir', type=str, default=None, help="Output directory "
            "where the converted model will be written.\n" 
            "Defaults to the input obj file folder if not provided.")
    parser.add_argument('--dont_get_face', action='store_false', help='Associate face id to respective face name (see get_face_names.py).')

    return parser.parse_args()

#==========================================
#region          HIGH-LEVEL FUNCTIONS
#==========================================

#(1) Convert .obj file to .vtp file using Paraview.
#(2) Add 1 to 'GroupIds' and replace it with 'ModelFaceID' to support face identication in Simvascular.
def convert_obj_to_vtp(obj_path,out_dir,model_type="ignore"):
    obj_fname = os.path.basename(obj_path)
    vtp_path = os.path.join(out_dir,obj_fname.replace(".obj",".vtp"))

    
    obj_src = WavefrontOBJReader(registrationName=obj_fname, FileName=obj_path)

    print_header("Make ModelFaceID array and save vtp file")

    # create ModelFaceID using Calculator (GroupIds + 1)
    calculator1 = Calculator(registrationName='Calculator1', Input=obj_src)
    calculator1.AttributeType = 'Cell Data'
    calculator1.ResultArrayName = 'ModelFaceID'
    calculator1.Function = 'GroupIds+1'
    calculator1.ResultArrayType = 'Int'

    #Compute value
    SetActiveSource(calculator1)

    #Get Array Names 
    calculator1_vtk = paraview.servermanager.Fetch(calculator1)
    point_data_array_names, _ = get_array_names(calculator1_vtk,data_type="PointData")
    cell_data_array_names, _ = get_array_names(calculator1_vtk,data_type="CellData")

    if model_type == "fluid" or model_type== "solid":
        faces_range = calculator1_vtk.GetCellData().GetArray("ModelFaceID").GetRange(0)
        print("ModelFaceID Range: {}".format(faces_range))
        compare_nb_faces(model_type,faces_range[1],return_values=False)

    #Remove GroupIds from arrays to save 
    cell_data_array_names = [x for x in cell_data_array_names if x != "GroupIds"]

    #Save vtp file, write other arrays by default
    print("Saving {}".format(vtp_path))
    SaveData(vtp_path, calculator1,ChooseArraysToWrite=1,CellDataArrays=cell_data_array_names)
    print_header("[TEST]")
    print("Loading {}".format(vtp_path))
    load_vtp_and_get_arrays(vtp_path)
    
    return calculator1, calculator1_vtk, faces_range[1]


#def 
#endregion         END - HIGH-LEVEL FUNCTIONS
#==========================================

#==========================================
#region          UTILS FUNCTIONS
#==========================================

#Check if the fluid/solid has the correct number of faces
def compare_nb_faces(model_type,nb_faces,return_values=False):
    correct_nb_faces = False
    if model_type == "fluid":
        expected_nb = 13 
    elif model_type == "solid":
        expected_nb = 14 
    else:
        warnings.warn("Model type '{}' not supported. Please choose between 'fluid' and 'solid'".format(model_type))

    if nb_faces == expected_nb:
        correct_nb_faces = True 
    else: 
        warnings.warn("Expected {} faces for {}, got {}.".format(expected_nb,model_type,nb_faces))

    
    #Could split into 2 functions, leaving as is for now 
    if return_values: 
        return expected_nb, nb_faces
    else: 
        return correct_nb_faces
#endregion         END - UTILS FUNCTIONS
#==========================================


if __name__ == '__main__':
    args = _parse_args()

    if args.out_dir is None: 
        out_dir = args.input_file_obj
    else:
        out_dir = args.out_dir

    if os.path.isdir(args.input_file_obj):
        input_file_objs = glob.glob(os.path.join(args.input_file_obj, '*.obj'))
    else:
        # out_dir = os.path.dirname(out_dir)
        input_file_objs = [args.input_file_obj]

    #Process all obj files in directory or the single provided 
    for input_file_obj in input_file_objs:
        print("--------------------------------------------------")
        print("--------------------------------------------------")
        calculator1, calculator1_vtk,nb_faces = convert_obj_to_vtp(input_file_obj,out_dir,args.model_type)
        
        #Optional, get be done as separete script
        if args.dont_get_face:
            input_fname = os.path.basename(input_file_obj)
            face_name_template_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),"face_names_template.txt")
            get_face_names(calculator1,model_type=args.model_type,out_dir=out_dir,
                fname_prefix=input_fname.replace(".obj",""), face_name_template_path=face_name_template_path)
        print()
