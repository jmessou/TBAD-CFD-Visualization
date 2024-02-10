//Only tested in Linux 
// Wrapper executable that takes 
// (1) The output folder of svpost (exportDir, Default: ., content: all_results*vtp and all_results*vtu)
// (2) The original job directory (jobDir, Default: .., content: mesh-complete),
//  then calls sv4guiSimulationUtils::CreateFlowFiles()
// that is in the Simvascular software
#include <string>
#include <iostream>
#include <filesystem>
#include <regex>

#include "sv4gui_SimulationUtils.h" //Simvascular library with flow calculations

namespace fs = std::filesystem;
using namespace std;

// bool regex_search_str(std::string const & pattern)
// {
//     std::regex self_regex(pattern,
//             std::regex_constants::ECMAScript | std::regex_constants::icase);
    
// }

//std::vector<std::string> meshFaceFileNames;
//Get all the filenames in a directory that followa pattern 
//Note: Default is ECMAScript regex 
std::vector<std::string> GetPatternFileNames(std::string const & dataDir, std::string const & fpattern, 
                                        bool filename_only = false)
{
    std::vector<std::string> patternFileNames;
    fs::path entry_path; 
    std::regex fpattern_regex(fpattern);

    cout << "Exploring Directory: " << dataDir << endl;
    cout << "Pattern: " << dataDir << endl;
    cout << endl;

    if (fs::is_directory(dataDir))
    {
        for (const auto & entry : fs::directory_iterator(dataDir))
        {
            entry_path = entry.path();
            //std::cout << entry_path.generic_string() << std::endl;
            if (std::regex_search(entry_path.generic_string(), fpattern_regex)) //Only add vtp files
            {
                std::cout << entry_path << std::endl;
                if (filename_only)
                    entry_path = entry_path.filename();
                patternFileNames.push_back(entry_path.generic_string());
            }
        }
    }

    cout << "--------------------------------------------" << endl;

    return patternFileNames;
}

void compute_flow(std::string const & exportDir, std::string const & jobDir)
{
    // Input Dirs/Files
    string meshFaceDir=jobDir+"/mesh-complete/mesh-surfaces";
    //Escape onces for C, other \ has to stay in pattern string
    std::vector<std::string> meshFaceFileNames = GetPatternFileNames(meshFaceDir,".*\\.vtp",true);
    std::vector<std::string> vtxFilePaths = GetPatternFileNames(exportDir,"all_results_.*\\.vtp");
    // Output Files
    string outPressureFilePath=exportDir+"/all_results-pressures.txt";
    string outFlowFilePath=exportDir+"/all_results-flows.txt";
    string outAverageFilePath=exportDir+"/all_results-averages.txt";
    string outAverageUnitsFilePath=exportDir+"/all_results-averages-from_cm-to-mmHg-L_per_min.txt";

    // Bools/others
    bool useComboFile = false; 
    bool skipWalls = false; 
    string unit = "mm";
    bool calculateFlows = false;

    if (meshFaceFileNames.size() > 0 && vtxFilePaths.size()  > 0)
    {
        string unit2 = "mm";
        calculateFlows = sv4guiSimulationUtils::CreateFlowFiles(outFlowFilePath, outPressureFilePath
                                                        , outAverageFilePath, outAverageUnitsFilePath
                                                        , vtxFilePaths,useComboFile
                                                        , meshFaceDir, meshFaceFileNames
                                                        , unit, skipWalls);
    }

    if (calculateFlows)
        cout << "Flow calculations were successful!" << endl;
    else
        cout << "Flow calculations FAILED!" << endl;
}

// ============
// MAIN PROGRAM
// ============
int main(int argc, char* argv[])
{
    string tmpstr; 
    bool RequestedHelp; 
    int iarg, arglength;

    // default directory is current directory
    char exportDir[255];
    exportDir[0]='.';
    exportDir[1]='/';
    exportDir[2]='\0';

    // default directory is upper directory
    char jobDir[255];
    jobDir[0]='.';
    jobDir[1]='.';
    jobDir[2]='/';
    jobDir[3]='\0';

    /* argc is the number of strings on the command-line */
    /*  starting with the program name */
    for(iarg=1; iarg<argc; iarg++){
        arglength = strlen(argv[iarg]);
        /* replace 0..arglength-1 with argv[iarg] */
        tmpstr.replace(0,arglength,argv[iarg],0,arglength);
        if(tmpstr=="-h"){
            RequestedHelp = true;
            cout << endl;
            cout << "usage:" <<endl;
            cout << "  compflow -e exportDir -j jobDir" << endl;
            cout << endl;
            cout << "COMMAND-LINE ARGUMENT SUMMARY" << endl;
            cout << "  -h             : Display usage and command-line argument summary" << endl;
            cout << "  -e exportDir   : Directory containing result files (default .)" << endl;
            cout << "  -j jobDir      : Directory containing the mesh-complete folder (default ..)"<< endl;
            cout << "END COMMAND-LINE ARGUMENT SUMMARY" << endl;
            cout << endl;
        }
        else if(tmpstr=="-e"){
            iarg++;
            exportDir[0]='\0';
            sprintf(exportDir,"%s/",argv[iarg]);
        }
        else if(tmpstr=="-j"){
            iarg++;
            jobDir[0]='\0';
            sprintf(jobDir,"%s/",argv[iarg]);
        }
        /* reset tmpstr for next argument */
        tmpstr.erase(0,arglength);
    }
    if(RequestedHelp){
        cout << endl;
        cout << "Exiting before due to -h flag";
        cout << endl;
        return(0);
    }
    cout << endl;
    cout << "Running compute_flow using: " << endl;
    cout << "    exportDir: " << exportDir << endl;
    cout << "    jobDir: " << jobDir << endl;
    cout << endl;
    cout << "--------------------------------------------" << endl;
    compute_flow(exportDir, jobDir);
    
    return 0;
}

