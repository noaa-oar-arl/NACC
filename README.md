# The NOAA-EPA Atmosphere-Chemistry Coupler (NACC)

The NOAA-EPA Atmosphere-Chemistry Coupler (NACC) is adapted from the Meteorology-Chemistry Interface Processor (MCIP), and can ingest output from the [Finite Volume Cubed Sphere (FV3)](https://www.gfdl.noaa.gov/fv3/) version of the [Global Forecast System (GFS)](https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/global-forcast-system-gfs), Regional (i.e., Limited Area Model; LAM) [FV3-based Short Range Weather (SRW)-Application](https://ufscommunity.org/srwa/), and the [Weather Research and Forecasting (WRF) Model](http://www.wrf-model.org) to prepare the meteorology files that are used within the CMAQ Modeling System. Where possible, MCIP/NACC uses data directly from the meteorological model to maximize consistency with the CMAQ Modeling System. When specific atmospheric fields are not explicitly output by WRF or FV3-based systems, MCIP/NACC uses scientific algorithms to create those fields for CMAQ.  MCIP/NACC output is used by the emissions model (for example, to provide time-varying temperatures for mobile emissions) and by the CCTM to define the atmospheric conditions. A scientific overview of MCIP is in [Otte and Pleim (2010)](https://www.geosci-model-dev.net/3/243/2010/).  

NACC was developed by Dr. Patrick C. Campbell and Dr. Youhua Tang, with support from the entire NOAA-ARL research team (Contact:  Patrick.C.Campbell@noaa.gov).  Unless otherwise necessary, hereafter the documentation of the MCIP-based NACC code will only be referred to as the NACC.

When using NACC, please appropriately cite the following:

 - Campbell, P. C., Tang, Y., Lee, P., Baker, B., Tong, D., Saylor, R., Stein, A., Huang, J., Huang, H.-C., Strobach, E., McQueen, J., Pan, L., Stajner, I., Sims, J., Tirado-Delgado, J., Jung, Y., Yang, F., Spero, T. L., and Gilliam, R. C.: Development and evaluation of an advanced National Air Quality Forecasting Capability using the NOAA Global Forecast System version 16, Geosci. Model Dev., 15, 3281–3313, https://doi.org/10.5194/gmd-15-3281-2022, 2022.

NACC performs the following functions using the output (history) file from WRF or FV3-based GFS and SRW-App (LAM) systems:

-   Defines the computational domain for the CCTM. The CCTM typically uses a smaller computational domain than the meteorological model, and the lateral boundary cells from the meteorological model generally are not used by CCTM.

-   Extracts meteorological model output on the computational domain that is prescribed for the CCTM.

-   Processes all required meteorological fields for the emissions model and the CCTM. Meteorological fields such as atmospheric temperature, pressure, humidity, and winds are acquired directly from the meteorological model (i.e., "passed through").

-   Uses the available meteorological fields to compute additional fields that are required by the CCTM but are not part of the meteorological model's output stream, such as the Jacobian which is used for coordinate transformations.

-   Outputs files that contain meteorological and geospatial information used by the emissions model and the CCTM.  The output can be either in I/O API or netCDF.

NACC is written in FORTRAN, and this version runs on a single (serial; WRF or FV3GFS) or multiple (parallel; only for FV3GFS-based) processors in a Unix/Linux environment. NACC is driven by a C-shell or K-shell script with several run-time options that are defined through a FORTRAN namelist. It is typical to use NACC to process hourly output fields from the meteorological model for each one-day period.

NACC is often updated concurrently with the CCTM.  The changes to NACC are documented with each update to the software, and a "Frequently Asked Questions" (FAQ) file exists that is specific to NACC.

As of NACCv2.0.0 (based on MCIPv5.0), WRF, FV3-GFS, and FV3-SRW Appi (LAM) based models are the only meteorological models that can be processed with NACC, but NACC could be expanded to process data from other meteorological models.

NACC can be used to determine the spatial region that is processed by CMAQ. NACC can process the full meteorological modeling domain, uniformly trim cells from that domain, or "window" a rectilinear subset of that domain. Configuration options for NACC include the time periods over which to extract data from the meteorological model output files, horizontal and vertical grid definitions, native or collapsed vertical layer definition and selections for integrating satellite cloud observations into NACC output.

## Files, configuration, and environment variables

All NACC configurations are established at run-time (rather than at compile time) via Fortran namelist variables rather than environment variables, which is a distinction from the rest of the CMAQ programs. The user does not need to directly edit the NACC namelist file. All configuration settings are contained in the NACC run script (run_nacc*.csh or run_nacc*.ksh), which automatically creates a new namelist file each time the script is executed.  The NACC output files are listed in Table 1, and the NACC output files are listed in Table 2.


## Compilation Configuration

All model configuration options for NACC are set during execution. System compiler options must be set in the provided Linux Makefile to build the program for different operating system/compiler combinations. Example compiler paths, flags, and library locations are provided in the default Makefile.


## Execution Configuration Variables

The variables listed here are set by the user in the NACC run script, and they are used during execution of the program.

-   `APPL [default: None]`  
    Application name; scenario ID for file naming
-   `CoordName [default: None]`  
    Coordinate system name of the NACC output grid that is written to the GRIDDESC file. Additional information about the parameters in the GRIDDESC file can be found in the [I/O API Documentation](https://www.cmascenter.org/ioapi/documentation/all_versions/html/GRIDDESC.html).
-   `GridName [default: None]`  
    Model grid name of the NACC output grid that is written to the GRIDDESC file. 
-   `DataPath [default: $CMAQ_DATA]`  
    Input/output data directory path
-   `InMetDir [default: None]`  
    Path of the input data directory containing the WRF‑ARW, FV3-GFS, or FV3-SRW App (LAM) output data files
-   `InGeoDir [default: None]`  
    Path of the input data directory containing the WRF Geogrid file, or similar pre-processed "geofile" (e.g., LAI, LANDUSEF) used for FV3-GFS and FV3-SRW-App (LAM)
-   `InVIIRSDir [default: None]`  
    Path of the input data directory containing VIIRS input file   
-   `OutDir [default: None]`  
    Path of the NACC output data directory
-   `ProgDir [default: None]`  
    Working directory containing the NACC executable
-   `WorkDir [default: $OutDir]`  
    Temporary working directory for Fortran links and the namelist file
-   `InMetFiles [WRF;  default: None]`  
    List of input meteorology files for WRF, including the directory path for each file; without modifying NACC, up to 300 meteorological model output files are allowed as input to a single NACC execution
-   `file_mm [WRF or FV3;  default: None]`  
    Input 2d/3d meteorology files for WRF or FV3, including the directory path for each file; without modifying NACC, up to 300 meteorological model output files are allowed as input to a single NACC execution
-   `file_sfc [FV3;  default: None]`  
    Input 2d-only surface meteorology files for FV3, including the directory path for each file; without modifying NACC, up to 300 meteorological model output files are allowed as input to a single NACC execution
-   `IfGeo [default: F]`  
    Binary flag indicating the availability of an input WRF or FV3-GFS Geogrid file; options include T (true) or F (false)
-   `InGeoFile [WRF; default: None]`  
    Name and location of input WRF Geogrid file
-   `file_geo [FV3;  default: None]`  
    Input geographic file for WRF or FV3, including the directory path for each file; without modifying NACC, up to 300 meteorological model output files are allowed as input to a single NACC execution
-   `InVIIRSFile [FV3; default: None]`  
    Name and location of input VIIRS File
-   `file_viirs_gvf [FV3;  default: None]`  
    Input VIIRS File
-   `InMetModel [default: 2]`
 Choose input meteorological model.
    -   `2`: WRF-ARW
    -   `3`: FV3-GFS
    -   `4`: FV3-SRW App (LAM)
-   `dx_out [FV3;  default: 12000 m]`
    Output dx grid met resolution in meters (must match DX in `domains` setting below)
-   `dy_out [FV3;  default: 12000 m]`
    Output dy grid met resolution in meters (must match DY in `domains` setting below)
-   `LPV: [default: 0]`  
    Compute and output potential vorticity. This must be activated to support the [CCTM O3 potential vorticity scaling](../../CCTM/docs/ReleaseNotes/Potential_Vorticity_Scaling.md).
    -   `0`: Do not compute and output potential vorticity
    -   `1`: Compute and output potential vorticity
-   `LWOUT [default: 0]`  
    Output vertical velocities.
    -   `0`: Do not output vertical velocity
    -   `1`: Output vertical velocity
-   `LUVBOUT [default: 0]`  
    Output u- and v-component winds on B staggered grid.
    -   `0`: Do not output u- and v-component winds on B-grid
    -   `1`: Output u- and v-component winds on B-grid (in addition to the C-grid)
-   `IFDIAG_PBL [default: False]`  
    Logical to choose if PBLH is taken from FV3 or rediagnosed by NACC.
    -   `False`: Use FV3 PBLH
    -   `True`: Recalculate PBLH in NACC
-   `IFVIIRS_GVF [default: False]`  
    Logical to choose if the vegetation fraction is taken from FV3 or VIIRS GVF.
    -   `False`: Use FV3 vegetation fraction
    -   `True`: Use VIIRS GVF 
-   `IFVIIRS_LAI [default: False]`
    Logical to choose if the LAI is taken from FV3 or VIIRS.
    -   `False`: Use FV3 LAI
    -   `True`: Use VIIRS LAI
-   `IFFENGSHA_DUST [default: False]`
    Logical to choose if the NOAA-ARL Fengsha Windblown dust is used in CMAQ:  based on input soil parameters in "GeoFile" and output in METCRO2D file.
    -   `False`: Do not input/output Fengsha WB Dust soil parameters
    -   `True`: Input/output of Fengsha WB Dust soil parameters
-   `IFBIOSEASON [default: False]`
    Logical to choose if the NOAA-ARL time-dependend bioseasons withc is used in CMAQ:  based on summer/winter seasonal information and T2,  and is output in METCRO2D file.
    -   `False`: Do not calculate and output SEASON variable in METCRO2D
    -   `True`:  Calculate and output SEASON variable in METCRO2D
-   `IFCANOPY [default: False]`
    Logical to choose if the ECCC in-canopy shading effects are used in CMAQ:  based on input canopy parameters in "GeoFile" and output in METCRO2D file.
    -   `False`: Do not input/output ECCC canopy parameters
    -   `True`: Input/output of ECCC canopy parameters
-   `MCIP_START [format: YYYY-MM-DD-HH:MM:SS.SSSS]`  
    Beginning date and time (UTC) of data to output from NACC. The start date and time must be contained within the input data from WRF, FV3-GFS, or FV3-SRW App (LAM).
-   `MCIP_END [format: YYYY-MM-DD-HH:MM:SS.SSSS]`  
    End date and time (UTC) of data to output from NACC. The end date and time must be contained within the input data from WRF, FV3-GFS, or FV3-SRW App (LAM).
-   `INTVL [default: 60]`  
    Output interval in minutes. This setting determines the amount of model time contained in each output time step. The output interval for NACC can be less frequent than the incoming meteorological model output (e.g., process 30-minute data for CCTM from 15-minute WRF output).
-   `CTMLAYS [default: -1.0]`  
    Set CTM layers.  Should be in descending order starting at 1 and ending with 0.  There is currently a maximum of 200 layers allowed.
    To use all of the layers from the input meteorology without collapsing (or explicitly specifying), set CTMLAYS = -1.0.
-   `CUTLAY_COLLAPX [default: 0]`
    Set number of layers to cut off model top before collapsing. This should be used with caution, and could be used when model top is extremely high leading to unrealistic mass-weighted linear interpolation of topmost collapsed layers.
-   `MKGRID [default: T]`  
    Determines whether to output static (GRID) meteorology files
-   `IOFORM [default: 1]`  
    Choose output format.
    -   `1`: Models-3 I/O API
    -   `2`: netCDF
-   `BTRIM [default: 5]`
    The number of boundary points to remove on each of the four horizontal sides of the meteorology output to define the NACC output domain. Setting BTRIM = 0 will specify the maximum extent of the input meteorology domain. To remove the WRF lateral boundaries, set BTRIM = 5 (recommended).
    This setting affects the output NACC horizontal domain by reducing the input meteorology domain by 2*BTRIM + 2*NTHIK + 1, where NTHIK is the lateral boundary thickness (from the BDY files). The extra point reflects the conversion from the grid points (dot points) to grid cells (cross points).
    To crop a subset of the input meteorology ("window"), set BTRIM = -1; this setting causes BTRIM to be replaced by the information provided by X0, Y0, NCOLS, and NROWS (see below).
-   `X0 [used only if BTRIM = -1]`  
    The *x*-coordinate of the lower-left corner of the full NACC cross-point domain (including the NACC lateral boundary) based on the input WRF‑ARW domain. X0 refers to the east-west direction.
-   `Y0 [used only if BTRIM = -1]`  
    The *y*-coordinate of the lower-left corner of the full NACC cross-point domain (including the NACC lateral boundary) based on the input WRF‑ARW domain. Y0 refers to the north-south direction.
-   `NCOLSIN [used only if BTRIM = -1]`  
    Number of columns in the output NACC domain (excluding NACC lateral boundaries)
-   `NROWSIN [used only if BTRIM = -1]`  
    Number of rows in the output NACC domain (excluding NACC lateral boundaries)
-   `LPRT_COL [default: 0]`  
    Column cell coordinate for diagnostic outputs on the NACC modeling domain
-   `LPRT_ROW [default: 0]`  
    Row cell coordinate for diagnostic outputs on the NACC modeling domain
-   `WRF_LC_REF_LAT [optional; used only for Lambert conformal projections; default: -999.0]`  
    WRF Lambert Conformal reference latitude. Use this setting to force the reference latitude in the output NACC data. If not set, NACC will use the average of the two true latitudes.
-   `projparm [FV3GFS or FV3-SRW App (LAM) Inputs  Only; used  to define projection parameters for subset of global grid; example set: GDTYP=2, P_ALP=33., P_BET=45., P_GAM=-97., XCENT=-97., YCENT=40.]`  
    Required definition of grid projection parameters for FV3GFS or FV3-SRW App (LAM). Use this setting to subset the FV3-GFS global or FV3-SRW App (LAM) domain that is regridded to use in CMAQ.
    For FV3-GFS global inputs, GRDTYP options 2 (Lambert Conformal, e.g., Regional CMAQ ) and 6 (Polar Stereographic, e.g., Hemispheric CMAQ) are currently supported.  For FV3-SRW App (LAM) Inputs, only 
    GRDTYP option 2 is supported.
-   `domains [FV3GFS or FV3-SRW App (LAM) Inputs Only; used  to define domain grid information for subset of global grid; example set XORIG=-2508000., YORIG=-1716000., DX=12000., DY=12000., NROWS=442, NCOLS=265]`  
    Required definition of grid domain parameters for FV3GFS or FV3-SRW App (LAM). Use this setting to subset the FV3-GFS or FV3-SRW App (LAM) global domain that is regridded to use in CMAQ.
-   `ntimes [WRF, FV3GFS or FV3-SRW App (LAM) Inputs Only; default = 0]`  

    Number of times to process for the model
    
  
    
    
    

## Compiling and Running

**Compile NACC**

NACC is compiled with a Makefile. The configuration options in the Makefile include the compiler and compiler flags to use for building the executable. Note that this version of NACC is either serial (WRF, FV3GFS, or FV3-SRW App (LAM)) or parallelized code (FV3GFS or FV3-SRW App (LAM)), so MPI libraries are required for ONLY the parallel version.  The parallel version is meant to be used with the global FV3GFSv16, as it significantly speeds up the total IO throughput wall clock time. 
Note:  To read the full vertical resolution of the global FV3GFSv16, NACC must be compiled with a modified version of the IOAPI library that has an adjusted maximum possible input vertical levels of MXLAYS3 =  200 in PARMS3.EXT.  An example of the modified IOAPI library to use is included in the NACC/lib/ioapi-3.2 repository. See the [CMAS IOAPI Tutorial](https://www.cmascenter.org/ioapi/documentation/all_versions/html/TUTORIAL.html) to find the location of where this change to MXLAYS3 takes place.  The only other library needed to run NACC on your system is the NETCDF library. 

The Makefile is located in the directory with the NACC source code (e.g., `NACC/serial/src` or `NACC/parallel/src`). To compile NACC, simply invoke the Makefile at the command line:

```
./make |& tee make.nacc.log
```

To port NACC to different compilers, change the compiler names, locations, and flags in the Makefile.

**Run NACC**

Set the run script settings according to the execution configuration variables described above. Run NACC to produce meteorology input data for the CCTM:

```
cd $NACC_home/serial/scripts
./run-nacc-fv3.csh |& tee run_nacc.log
```


**Table 1. NACC input files**

|**File Name**|**Format**|**Description**|**Required**|
|------------|------------------------------|-----------------------------------------------------|---------------------|
|InMetFiles|netCDF (WRF, FV3-GFS, or FV3-SRW App (LAM))|List of WRF, FV3-GFS, or FV3-SRW App (LAM) output files for input to NACC|required|
|InSfcFiles|netCDF (FV3-GFS or FV3-SRW App (LAM))|List of FV3-GFS or FV3-SRW App (LAM) output files for input to NACC|required (only FV3-GFS or FV3-SRW App (LAM))|
|InGeoFile|netCDF (WRF, FV3-GFS, or FV3-SRW App (LAM))|Output from WRF Geogrid processor | optional; only required if fractional land use, LAI, etc are not part of the WRF, FV3-GFS, or FV3-SRW App (LAM) output.  Offline Pre-processed NOAA-ARL "geofiles" with LAI (VIIRS 2018-2020 climatology) and LANDUSEF (based on 12-month climatological IGBP-MODIS) for the global GFSv16 Gaussian NetCDF Grid are available via FTP by request (Contact:  Patrick C. Campbell; Patrick.C.Campbell@noaa.gov)|
|InVIIRSFile|netCDF (FV3-GFS or FV3-SRW App (LAM))|Input from VIIRS data |optional; only if global NetCDF VIIRS Input COARDS file is provided. Global ~ 4km VIIRS NetCDF Grid are available via FTP by request. (Contact:  Patrick C. Campbell; Patrick.C.Campbell@noaa.gov)|

**Table 2. NACC output files**

|**File Name**|**Format**|**Description**|**Required**|
|--------------------|-----------------|------------------------------------------------------------------|---------------------------|
|GRIDDESC|ASCII|Grid description file with coordinate and grid definition information|required|
|GRID_BDY_2D|I/O API|Time-independent 2-D boundary meteorology file|required|
|GRID_CRO_2D|I/O API|Time-independent 2-D cross-point meteorology file|required|
|GRID_CRO_3D|I/O API|Time-independent 3-D cross-point meteorology file|required|
|GRID_DOT_2D|I/O API|Time-independent 2-D dot-point meteorology file|required|
|LUFRAC_CRO|I/O API|Time-independent fractional land use by category|created if fractional land use was provided in WRF's, FV3-GFS's, or FV3-SRW App's (LAM) output or in Geogrid output|
|MET_BDY_3D|I/O API|Time-varying 3-D boundary meteorology file|required|
|MET_CRO_2D|I/O API|Time-varying 2-D cross-point meteorology file|required|
|MET_CRO_3D|I/O API|Time-varying 3-D cross-point meteorology file|required|
|MET_DOT_3D|I/O API|Time-varying 3-D dot-point meteorology file|required|
|MOSAIC_CRO|I/O API|Time-varying 3-D output from mosaic land use|created if the Noah Mosaic land-surface model was run in WRF|
|SOI_CRO|I/O API|Time-varying soil properties in each soil layer|created if a land-surface model was run in WRF, FV3-GFS, or FV3-SRW App (LAM)|

The location of the NACC output files is set in the NACC script using the `$OutDir` variable. The names of the NACC output files are generic and do not have any information about the model grid that they are simulating or the time period that is covered. These attributes can be controlled by the NACC script. For example, the name of the grid can be used in the output directory path. In addition, the default naming convention for all NACC output files appends the `APPL` environment variable to the file name to identify files by the time period that is represented by the file. All of the file naming variables for the NACC outputs are set in the run script, and they can be easily tailored to fit each user's application or style.

***Current Version of NACC***

NACCv2.1.2 (Frozen Jun 01, 2023)

**Previous Versions of NACC**
v2.1.1 (Apr 11, 2023)
v2.1.0 (Dec 14, 2022)
v2.0.0 (Jul 12, 2022)
v1.3.2 (Aug 04, 2021)
v1.3.1 (Dec 07, 2021)
v1.3.0 (Nov 24, 2021)
v1.2.1 (Nov 10, 2021)
v1.2.0 (Sep 16, 2021)
v1.1.0 (Jul 06, 2020)
v1.0.0 (May 06, 2020)
