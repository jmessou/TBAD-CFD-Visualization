import configparser 

#==========================================
#region          config CLASS
#==========================================

#Reads config file and sends content to main file
#Main file less bulky this way 
class ConfigAll:
    def __init__(self,config_file='example.ini'):
        self.config = configparser.ConfigParser()
        self.config.read(config_file)
        self.fc = FlowControl() 
        self.set_flow()

    def set_flow(self):
        self.fc.slice_centerline_bool = self.config['FLOW'].getboolean("slice_centerline_bool") 
        self.fc.load_and_slice_centerline_bool = self.config['FLOW'].getboolean("load_and_slice_centerline_bool") 
        self.fc.process_slice_bool = self.config['FLOW'].getboolean("process_slice_bool") 
        self.fc.convert_results_bool = self.config['FLOW'].getboolean("convert_results_bool")
        self.fc.save_centerline_bool = self.config['FLOW'].getboolean("save_centerline_bool") 
        self.fc.save_slices_bool =self.config['FLOW'].getboolean("save_slices_bool") 
        self.fc.save_all_times_bool = self.config['FLOW'].getboolean("save_all_times_bool") 
        self.fc.save_model_bool = self.config['FLOW'].getboolean("save_model_bool")
        self.fc.save_osi_tawss_bool = self.config['FLOW'].getboolean("save_osi_tawss_bool") 
        self.fc.do_streamlines_work_bool = self.config['FLOW'].getboolean("do_streamlines_work_bool") 

        if self.config['FLOW'].getboolean("no_show"):
            self.fc.no_show()
        if self.fc.process_slice_bool and self.config['ANALYSIS'].getboolean("loop_timestep"):
            print("[WARNING] loop_timestep = True, disabling all visualization.")
            self.fc.no_show()

    def get_FlowControl(self):
        return self.fc

    def get_general(self):
       return self.config['GENERAL'].get("patient")

    def get_paths(self):
        return self.config['PATHS']["Centerline"], \
        self.config['PATHS']["VTU_folder"], \
        self.config['PATHS']["Slices_tsv"], \
        self.config['PATHS']["Slices_tsv_all"], \
        self.config['PATHS'].get("output_dir"), \
        self.config['PATHS'].get("true_lumen_path"), \
        self.config['PATHS'].get("etfn_path")
    
    def get_analysis(self):
        file_range = self.config['ANALYSIS'].get("file_range").split(",")
        file_range = [int(float(x)) for x in file_range]
        analysis_time_step = self.config['ANALYSIS'].getint("analysis_time_step")
        nb_files = self.config['ANALYSIS'].getint("nb_files")
        loop_timestep = self.config['ANALYSIS'].getboolean("loop_timestep")
        file_stride = self.config['ANALYSIS'].getint("file_stride")
        return file_range, loop_timestep, analysis_time_step, nb_files, file_stride
    
    def get_slice(self):
        return self.config['SLICE PROCESSING'].getint("slice_stride",fallback=10) 

    def get_view(self):
        return self.config['VIEW'].getfloat("azimuth"),  \
        self.config['VIEW'].getfloat("elevation")
    
    #velocity/osi/tawss colorbar parameters
    def get_vis(self):
        cb_mode_osi = "visible"
        palette = self.config['VISUALIZATION'].get("palette",fallback="rainbow")
        color_bar_mode_vel = self.config['VISUALIZATION'].get("color_bar_mode_vel",fallback="visible")
        cb_mode_tawss = self.config['VISUALIZATION'].get("bar_range",fallback="bar_range")
        #Only used when mode is bar_range
        bar_range_tawss_min = self.config['VISUALIZATION'].getfloat("bar_range_tawss_min",fallback=0)
        bar_range_tawss_max = self.config['VISUALIZATION'].getfloat("bar_range_tawss_max",fallback=2)
        bar_range_tawss = [bar_range_tawss_min, bar_range_tawss_max]
        bar_range_vel_max = self.config['VISUALIZATION'].getfloat("bar_range_vel_max",fallback=1.6)
        return palette, cb_mode_osi, color_bar_mode_vel, cb_mode_tawss, bar_range_tawss, bar_range_vel_max

    def get_slice_vis(self):
        return self.config['SLICE VISUALIZATION'].getfloat("wall_opacity"), \
        self.config['SLICE VISUALIZATION'].get("cross_sec").split(",")

    def get_model_vis(self):
        return self.config['MODEL VISUALIZATION'].getboolean("show_thrombus")

    def get_stream(self):
        return self.config['STREAMLINES'].getfloat("stream_radius"), \
        self.config['STREAMLINES'].getfloat("tube_radius"), \
        self.config['STREAMLINES'].getfloat("stream_max_length"), \
        self.config['STREAMLINES'].getint("stream_nb_pts"), \
        self.config['STREAMLINES'].getfloat("wall_opacity_vel")

    def get_saving(self): 
        image_res = self.config['SAVING'].get("image_res").split(",")
        image_res = [int(float(x)) for x in image_res]
        return image_res
        
    def save_config(self,output_file='example2.ini',config=None):
        if config is None:
            config = self.config
        with open(output_file, 'w') as configfile:
            config.write(configfile)
            
#endregion         config CLASS
#==========================================

#==========================================
#region          Flow control CLASS
#==========================================
class FlowControl:
    def __init__(self):
        #Running Options
        self.save_centerline_bool = True #Show and save centerline 
        self.slice_centerline_bool = True #Will slice center line and analyze cross sections
        self.load_and_slice_centerline_bool = True #Will load slices info and analyze cross sections
        self.process_slice_bool = True #Compute flow/pressure at slices
        self.save_slices_bool =True #Show and save slices 
        self.save_all_times_bool = True #Save slices at all timesteps if save_slices_bool is on
        self.convert_results_bool = True
        self.save_model_bool = True #Show and save model with cross-sections and TL/FL
        self.save_osi_tawss_bool = True #Show and save OSI/TAWSS results
        self.do_streamlines_work_bool = True #Streamlines and old work

    def disable_all(self):
        #Running Options
        self.slice_centerline_bool = False #Will slice center line and analyze cross sections
        self.load_and_slice_centerline_bool = False #Will load slices info and analyze cross sections
        self.process_slice_bool = False #Compute flow/pressure at slices
        self.convert_results_bool = False
        self.save_slices_bool =False #Show and save slices 
        self.save_all_times_bool = False #Save slices at all timesteps if save_slices_bool is on
        self.save_centerline_bool = False #Show and save centerline 
        self.save_model_bool = False #Show and save model with cross-sections and TL/FL
        self.save_osi_tawss_bool = False #Show and save OSI/TAWSS results
        self.do_streamlines_work_bool = False #Streamlines and old work

    def enable_all(self): 
        #Running Options
        self.slice_centerline_bool = True #Will slice center line and analyze cross sections
        self.convert_results_bool = True
        self.load_and_slice_centerline_bool = True #Will load slices info and analyze cross sections
        self.process_slice_bool = True #Compute flow/pressure at slices
        self.save_slices_bool =True #Show and save slices 
        self.save_all_times_bool = True #Save slices at all timesteps if save_slices_bool is on
        self.save_centerline_bool = True #Show and save centerline 
        self.save_model_bool = True #Show and save model with cross-sections and TL/FL
        self.save_osi_tawss_bool = True #Show and save OSI/TAWSS results
        self.do_streamlines_work_bool = True #Streamlines and old work

    def no_show(self): #Do work that doesn't need any view 
        #Running Options
        self.save_slices_bool =False #Show and save slices 
        self.save_all_times_bool = False #Save slices at all timesteps if save_slices_bool is on
        self.save_centerline_bool = False #Show and save centerline 
        self.save_model_bool = False #Show and save model with cross-sections and TL/FL
        self.save_osi_tawss_bool = False #Show and save OSI/TAWSS results
        self.do_streamlines_work_bool = False #Streamlines and old work 

    def __str__(self):
        return '''====Flow Control==== \
        \nslice_centerline: {} | load_and_slice_centerline: {} \
        \nprocess_slice_bool: {} | convert_results: {} \
        \nsave_slices: {} | save_all_times: {} | save_centerline: {} \
        \nsave_model_bool: {} | save_osi_tawss: {} | do_velocity_work: {} \
        \n-----------------'''.format(self.slice_centerline_bool, self.load_and_slice_centerline_bool, 
        self.process_slice_bool, self.convert_results_bool,
        self.save_slices_bool, self.save_all_times_bool, self.save_centerline_bool,
        self.save_model_bool, self.save_osi_tawss_bool, self.do_streamlines_work_bool)

#endregion         Flow Control CLASS
#==========================================