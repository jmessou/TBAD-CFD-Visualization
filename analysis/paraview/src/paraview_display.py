import paraview 
from paraview.simple import *
from paraview_calc import threshold_face_id

#==========================================
#region          DISPLAY CLASSES
#==========================================

class ColorBar(): 
    def __init__(self,name,min_val=0,max_val=1.6,mode="visible",palette="rainbow"):
        self.name = name
        self.min = min_val
        self.max = max_val
        self.palette = palette
        self.mode = mode
        self.update_range()

    def update_range(self): 
        self.range = [self.min,self.max] 

    def __str__(self):
        return f"[{self.name}] | Mode: {self.mode} | Min: {self.min} | Max: {self.max} | Palette: {self.palette}"

#Reference to a specific region that needs to be shown 
class RegionToShow(): 
    def __init__(self,name,src,opacity=1,color=None):
        self.name = name
        self.src = src
        self.opacity = opacity
        self.color = color

    def __str__(self):
        return f"[Region] {self.name} | Opacity: {self.opacity} | Color: {self.color}"

#endregion         END - DISPLAY CLASSES
#==========================================

#==========================================
#region          DISPLAY FUNCTIONS
#==========================================

#TODO Use format_colorbar or ColorBy directly
def show_and_color(src,myview=None,opacity=None,show_bar=True,array_name="Velocity [m/s]"): 
    if myview is None:
        myview = GetActiveView()

    srcDisplay = Show(src, myview, 'UnstructuredGridRepresentation')
    srcDisplay.SetScalarColoring(array_name, servermanager.GetAssociationFromString('POINTS'))
    # ColorBy(srcDisplay, value=('POINTS', array_name))
    if opacity is not None:
        srcDisplay.Opacity = opacity
    # srcLUT = GetColorTransferFunction("velocity [m/s]")
    # srcLUTColorBar = GetScalarBar(srcLUT, myview)
    # srcLUTColorBar.RangeLabelFormat = '%-#6.1f'
    srcDisplay.SetScalarBarVisibility(myview, show_bar)

    return srcDisplay 

def disable_light(disp):
    print(disp.Diffuse)
    print(disp.Ambient)
    disp.Diffuse = 0.8
    disp.Ambient = 0.2

#Set a solid color
#For surfaces, use disable_light to disable it darker surfaces when the object if seen from below
def set_solid_col(disp,opacity=1,color="",rep="Surface",disable_light=False):
    disp.Representation = rep
    
    disp.ColorArrayName = ['POINTS', '']
    disp.DiffuseColor = get_color(color)
    disp.AmbientColor = get_color(color)
    if disable_light:
        disp.Diffuse = 0
        disp.Ambient = 1
    disp.Opacity = opacity   
    return disp

#TODO, paraview function instead
#Might be a paraview function for it 
def get_color(color): 
    colors = {
    'red': [1.0, 0.0, 0.0],
    'green_real': [0.0, 1.0, 0.0],
    'blue': [0.0, 0.0, 1.0],
    'yellow': [1.0, 1.0, 0.0],
    'orange': [1.0, 0.5, 0.0],
    'purple_real': [0.5, 0.0, 0.5],
    'cyan': [0.0, 1.0, 1.0],
    'magenta': [1.0, 0.0, 1.0],
    'white': [1.0, 1.0, 1.0],
    'black': [0.0, 0.0, 0.0],
    'grey': [0.65, 0.65, 0.65],
    'green': [0.0, 1.0, 0.5],
    'purple': [0.6666666666666666, 0.3333333333333333, 1.0]
    }
    return colors[color]  
#Will set up the color bar blue to red rainbow for the source "src_name"
#and will impose the range bar_range 
#Example: format_colorbar("velocityms", [0.0,0.8,0.0, 1.0])
#mode: data, visible (default), overtime, bar_range
def format_colorbar(src_name, bar_range=[0.0,1] ,preset="rainbow",
    color_by_array='Velocity [m/s]',title='Velocity Magnitude [m/s]',
    separate=False,opacity=1,mode="visible",show=True,show_bar=True,debug=False,
    input_src_color=None,threshold=False):
    
    if debug:
        print("{} | {} | Opacity: {}".format(color_by_array,mode,opacity))
    myview = GetActiveView()
    if threshold: 
        # disp_background  = Show(src_name) 
        # set_solid_col(disp_background,opacity=0.1,color="white")
        # Render() 
        # mode = "bar_range"
        # if color_by_array == "WSS_time_average":
        #     import pdb; pdb.set_trace()
        # src_threshold_solid = threshold_face_id(src_name,bar_range[1],color_by_array,array_name=color_by_array,cells=True,region_id_min=bar_range[0],method='Above Upper Threshold')
        # disp_background  = Show(src_threshold_solid) 
        # set_solid_col(disp_background,opacity=1,color="white")
        # Render()
        print(bar_range)
        # if color_by_array == "WSS_time_average":
        #     import pdb; pdb.set_trace()
        src_threshold = threshold_face_id(src_name,bar_range[1],color_by_array,array_name=color_by_array,cells=True,region_id_min=bar_range[0])
        src_name = src_threshold
    disp  = Show(src_name) 
    # if color_by_array == "WSS_time_average":
    #     import pdb; pdb.set_trace()
    if input_src_color is not None:
        # ColorBy(disp, None)
        set_solid_col(disp,opacity=opacity,color=input_src_color)
        disp_LUT = None
    else:
        ColorBy(disp, value=('POINTS', color_by_array),separate=separate)
        disp_LUT = GetColorTransferFunction(color_by_array, disp,separate=separate)

        #Apply preset
        if preset == "rainbow":
            preset = 'Blue to Red Rainbow'
        elif preset == "cool": 
            preset = "Cool to Warm"
        else: 
            preset = "Cool to Warm"
        disp_LUT.ApplyPreset(preset, True)
        # if color_by_array == "WSS_time_average":
            # import pdb; pdb.set_trace()
        set_scale_range(disp,mode=mode,bar_range=bar_range,disp_LUT=disp_LUT,debug=debug)
        
        disp.Opacity = opacity

        if show_bar:
            disp.SetScalarBarVisibility(myview, show_bar)
            format_bar_labels(GetScalarBar(disp_LUT, myview),title)

    if show is False:
        Hide(src_name)

    return disp, disp_LUT 

def set_scale_range(disp,mode="visible",bar_range=[0,1],disp_LUT=None,myview=None,debug=False):
    if myview is None:
        myview = GetActiveView()
        Render()
    #Scale Range 
    if mode == "visible": 
        disp.RescaleTransferFunctionToVisibleRange(myview)
    elif mode == "data":
        disp.RescaleTransferFunctionToDataRange(True)
    elif mode == "overtime":
        disp.RescaleTransferFunctionToDataRangeOverTime() 
    elif mode == "bar_range":
        if disp_LUT is not None:
            disp_LUT.RescaleTransferFunction(bar_range[0],bar_range[1])
        else:
            print("Please provide bar disp_LUT")
    else: 
        if debug:
            print("No Range Mode")


def format_bar_labels(srcLUTColorBar,title): 
    srcLUTColorBar.Title = title
    srcLUTColorBar.ComponentTitle = ''
    #srcLUTColorBar.TitleBold = 1
    #srcLUTColorBar.LabelBold = 1
    srcLUTColorBar.AutomaticLabelFormat = 0
    srcLUTColorBar.AddRangeLabels = 1
    if "osi" in title.lower() or "wss" in title.lower():
        label_format = '%-#6.1f' 
        srcLUTColorBar.HorizontalTitle = 1
        #label_format = '%.2e' #TODO Remove
    else:
        label_format = '%-#6.1f' 
        # label_format = '%.2e' #TODO Remove
    srcLUTColorBar.LabelFormat = label_format
    srcLUTColorBar.RangeLabelFormat = label_format
    # srcLUTColorBar.Orientation = 'Horizontal'
    # srcLUTColorBar.ScalarBarLength = 0.8

#Will set up the color bar blue to red rainbow for the source "src_name"
#and will impose the range bar_range 
#Example: format_colorbar("velocityms", [0.0,0.8,0.0, 1.0])
def format_colorbar2(src_name, bar_range=[0.0,1, 0.0, 1.0] ,palette="cool",title='Velocity Magnitude [m/s]'):
    # get 2D transfer function for 'velocityms'
    srcTF2D = GetTransferFunction2D(src_name)
    srcTF2D.ScalarRangeInitialized = 1
    srcTF2D.Range = bar_range 

    # get color transfer function/color map for 'src'
    srcLUT = GetColorTransferFunction(src_name)
    srcLUT.AutomaticRescaleRangeMode = 'Never'
    srcLUT.TransferFunction2D = srcTF2D
   
    if palette == "rainbow": 
        srcLUT.RGBPoints = [bar_range[0], 0.0, 0.0, 1.0, bar_range[1], 1.0, 0.0, 0.0]
        srcLUT.ColorSpace = 'HSV'
        srcLUT.NanColor = [0.498039215686, 0.498039215686, 0.498039215686]
    elif palette == "cool":
        srcLUT.RGBPoints = [bar_range[0], 0.23137254902, 0.298039215686, 0.752941176471, 0.4*(bar_range[1]-bar_range[0])+bar_range[0], 0.865, 0.865, 0.865, bar_range[1], 0.705882352941, 0.0156862745098, 0.149019607843]
        srcLUT.ColorSpace = 'RGB'
        srcLUT.NanColor = [1, 1, 0]
    else: #cool 
        srcLUT.RGBPoints = [bar_range[0], 0.23137254902, 0.298039215686, 0.752941176471, 0.4*(bar_range[1]-bar_range[0])+bar_range[0], 0.865, 0.865, 0.865, bar_range[1], 0.705882352941, 0.0156862745098, 0.149019607843]
        srcLUT.ColorSpace = 'RGB' 
        srcLUT.NanColor = [1, 1, 0] 

    srcLUT.ScalarRangeInitialized = 1.0

    # get opacity transfer function/opacity map for 'src'
    srcPWF = GetOpacityTransferFunction(src_name)
    srcPWF.Points = [bar_range[0], 0.0, 0.5, 0.0, bar_range[1], 1.0, 0.5, 0.0]
    srcPWF.ScalarRangeInitialized = 1
    
    renderView1 = GetActiveView()
    # get color legend/bar for srcLUT in view renderView1
    srcLUTColorBar = GetScalarBar(srcLUT, renderView1)
    srcLUTColorBar.WindowLocation = 'Any Location'
    srcLUTColorBar.Position = [0.8012956669498725, 0.3109090909090909]
    srcLUTColorBar.Title = title
    srcLUTColorBar.ComponentTitle = ''
    #srcLUTColorBar.TitleBold = 1
    #srcLUTColorBar.LabelBold = 1
    srcLUTColorBar.ScalarBarLength = 0.33000000000000007
    srcLUTColorBar.AutomaticLabelFormat = 0
    srcLUTColorBar.LabelFormat = '%-#6.1f'
    srcLUTColorBar.AddRangeLabels = 1
    srcLUTColorBar.RangeLabelFormat = '%-#6.1f'
    #srcLUTColorBar.Visibility = 1
    
    # show color legend
    #velocitymsDisplay.SetScalarBarVisibility(renderView1, True)
    

#endregion         END - DISPLAY FUNCTIONS
#==========================================