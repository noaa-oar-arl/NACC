# The NOAA-ARL Atmosphere-Chemistry Coupler (NACC)

The NOAA-ARL Atmosphere-Chemistry Coupler (NACC) is adapted from the Meteorology-Chemistry Interface Processor (MCIP), and can ingest output from the [Finite Volume Cubed Sphere (FV3)](https://www.gfdl.noaa.gov/fv3/) version of the [Global Forecast System (GFS)](https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/global-forcast-system-gfs) and the [Weather Research and Forecasting (WRF) Model](http://www.wrf-model.org) to prepare the meteorology files that are used within the CMAQ Modeling System. Where possible, MCIP uses data directly from the meteorological model to maximize consistency with the CMAQ Modeling System. When specific atmospheric fields are not explicitly output by WRF, MCIP uses scientific algorithms to create those fields for CMAQ.  MCIP output is used by the emissions model (for example, to provide time-varying temperatures for mobile emissions) and by the CCTM to define the atmospheric conditions. A scientific overview of MCIP is in [Otte and Pleim (2010)](https://www.geosci-model-dev.net/3/243/2010/).  

NACC was developed by Dr. Patrick C. Campbell and Dr. Youhua Tang, with support from the entire NOAA-ARL research team (Contact:  Patrick.C.Campbell@noaa.gov).  Unless otherwise necessary, hereafter the documentation of the MCIP-based NACC code will only be referred to as the NACC.

NACC performs the following functions using the output (history) file from WRF or FV3-GFS:

-   Defines the computational domain for the CCTM. The CCTM typically uses a smaller computational domain than the meteorological model, and the lateral boundary cells from the meteorological model generally are not used by CCTM.

-   Extracts meteorological model output on the computational domain that is prescribed for the CCTM.

-   Processes all required meteorological fields for the emissions model and the CCTM. Meteorological fields such as atmospheric temperature, pressure, humidity, and winds are acquired directly from the meteorological model (i.e., "passed through").

-   Uses the available meteorological fields to compute additional fields that are required by the CCTM but are not part of the meteorological model's output stream, such as the Jacobian which is used for coordinate transformations.

-   Outputs files that contain meteorological and geospatial information used by the emissions model and the CCTM.  The output can be either in I/O API or netCDF.

NACC is written in FORTRAN, and this version runs on a single (serial) or multiple (parallel) processors in a Unix/Linux environment. NACC is driven by a C-shell or K-shell script with several run-time options that are defined through a FORTRAN namelist. It is typical to use NACC to process hourly output fields from the meteorological model for each one-day period.

NACC is often updated concurrently with the CCTM.  The changes to NACC are documented with each update to the software, and a "Frequently Asked Questions" (FAQ) file exists that is specific to NACC.

As of NACCv1.0 (based on MCIPv5.0), WRF and FV3-GFS are the only meteorological models that can be processed with NACC, but NACC could be expanded to process data from other meteorological models.

NACC can be used to determine the spatial region that is processed by CMAQ. NACC can process the full meteorological modeling domain, uniformly trim cells from that domain, or "window" a rectilinear subset of that domain. Configuration options for NACC include the time periods over which to extract data from the meteorological model output files, horizontal and vertical grid definitions, native or collapsed vertical layer definition and selections for integrating satellite cloud observations into NACC output.

## Files, configuration, and environment variables

All NACC configurations are established at run-time (rather than at compile time) via Fortran namelist variables rather than environment variables, which is a distinction from the rest of the CMAQ programs. The user does not need to directly edit the NACC namelist file. All configuration settings are contained in the NACC run script (run_nacc*.csh or run_nacc*.ksh), which automatically creates a new namelist file each time the script is executed.  The NACC output files are listed in Table 1, and the NACC output files are listed in Table 2.


## Compilation Configuration

All model configuration options for NACC are set during execution. System compiler options must be set in the provided Linux Makefile to build the program for different operating system/compiler combinations. Example compiler paths, flags, and library locations are provided in the default Makefile.


## Execution Configuration Variables

The variables listed here are set by the user in the NACC script (run_mcip.csh), and they are used during execution of the program.

-   `APPL [default: None]`  
    Application name; scenario ID for file naming
-   `CoordName [default: None]`  
    Coordinate system name of the NACC output grid that is written to the GRIDDESC file. Additional information about the parameters in the GRIDDESC file can be found in the [I/O API Documentation](https://www.cmascenter.org/ioapi/documentation/all_versions/html/GRIDDESC.html).
-   `GridName [default: None]`  
    Model grid name of the NACC output grid that is written to the GRIDDESC file. 
-   `DataPath [default: $CMAQ_DATA]`  
    Input/output data directory path
-   `InMetDir [default: None]`  
    Path of the input data directory containing the WRF‑ARW or FV3-GFS output data files
-   `InGeoDir [default: None]`  
    Path of the input data directory containing the WRF Geogrid file, or similar pre-processed "geofile" (e.g., LAI, LANDUSEF) used for FV3-GFS
-   `OutDir [default: $CMAQ_HOME/data/mcip]`  
    Path of the NACC output data directory
-   `ProgDir [default: $CMAQ_HOME/PREP/mcip/src]`  
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
-   `InMetModel [default: 2]`
 Choose input meteorological model.
    -   `2`: WRF-ARW
    -   `3`: FV3-GFS
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
-   `MCIP_START [format: YYYY-MM-DD-HH:MM:SS.SSSS]`  
    Beginning date and time (UTC) of data to output from NACC. The start date and time must be contained within the input data from WRF or FV3-GFS.
-   `MCIP_END [format: YYYY-MM-DD-HH:MM:SS.SSSS]`  
    End date and time (UTC) of data to output from NACC. The end date and time must be contained within the input data from WRF or FV3-GFS.
-   `INTVL [default: 60]`  
    Output interval in minutes. This setting determines the amount of model time contained in each output time step. The output interval for NACC can be less frequent than the incoming meteorological model output (e.g., process 30-minute data for CCTM from 15-minute WRF output).
-   `CTMLAYS [default: -1.0]`  
    Set CTM layers.  Should be in descending order starting at 1 and ending with 0.  There is currently a maximum of 200 layers allowed.
    To use all of the layers from the input meteorology without collapsing (or explicitly specifying), set CTMLAYS = -1.0.
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
-   `projparm [FV3GFS-Only; used  to define projection parameters for subset of global grid; example set: GDTYP=2, P_ALP=33., P_BET=45., P_GAM=-97., XCENT=-97., YCENT=40.]`  
    Required definition of grid projection parameters for FV3GFS. Use this setting to subset the FV3-GFS global domain that is regridded to use in CMAQ.
-   `domains [FV3GFS-Only; used  to define domain grid information for subset of global grid; example set XORIG=-2508000., YORIG=-1716000., DX=12000., DY=12000., NROWS=442, NCOLS=265]`  
    Required definition of grid domain parameters for FV3GFS. Use this setting to subset the FV3-GFS global domain that is regridded to use in CMAQ.
-   `ntimes [WRF and FV3GFS; default = 0]`  
    Number of times to process for the model
    
  
    
    
    

## Compiling and Running

**Compile NACC**

NACC is compiled with a Makefile. The configuration options in the Makefile include the compiler and compiler flags to use for building the executable. Note that this version of NACC is either serial or parallelized code (time splitting), so MPI libraries are required for ONLY the parallel version.  The Makefile is located in the directory with the NACC source code (e.g., `NACC/serial/src`). To compile NACC, simply invoke the Makefile at the command line:

```
./make |& tee make.nacc.log
```

To port NACC to different compilers, change the compiler names, locations, and flags in the Makefile.

**Run NACC**

Set the run script settings according to the execution configuration variables described above. Run NACC to produce meteorology input data for the CCTM:

```
cd $CMAQ_HOME/PREP/mcip/scripts
./run-nacc-fv3.csh |& tee run_nacc.log
```


**Table 1. NACC input files**

|**File Name**|**Format**|**Description**|**Required**|
|------------|------------------------------|-----------------------------------------------------|---------------------|
|InMetFiles|netCDF (WRF or FV3-GFS)|List of WRF or FV3-GFS output files for input to NACC|required|
|InSfcFiles|netCDF (FV3-GFS)|List of FV3-GFS output files for input to NACC|required (only FV3-GFS)|
|InGeoFile|netCDF (WRFor FV3-GFS)|Output from WRF Geogrid processor|optional; only required if fractional land use, LAI, etc are not part of the WRF or FV3-GFS output|


**Table 2. NACC output files**

|**File Name**|**Format**|**Description**|**Required**|
|--------------------|-----------------|------------------------------------------------------------------|---------------------------|
|GRIDDESC|ASCII|Grid description file with coordinate and grid definition information|required|
|GRID_BDY_2D|I/O API|Time-independent 2-D boundary meteorology file|required|
|GRID_CRO_2D|I/O API|Time-independent 2-D cross-point meteorology file|required|
|GRID_CRO_3D|I/O API|Time-independent 3-D cross-point meteorology file|required|
|GRID_DOT_2D|I/O API|Time-independent 2-D dot-point meteorology file|required|
|LUFRAC_CRO|I/O API|Time-independent fractional land use by category|created if fractional land use was provided in WRF's or FV3-GFS's output or in Geogrid output|
|MET_BDY_3D|I/O API|Time-varying 3-D boundary meteorology file|required|
|MET_CRO_2D|I/O API|Time-varying 2-D cross-point meteorology file|required|
|MET_CRO_3D|I/O API|Time-varying 3-D cross-point meteorology file|required|
|MET_DOT_3D|I/O API|Time-varying 3-D dot-point meteorology file|required|
|MOSAIC_CRO|I/O API|Time-varying 3-D output from mosaic land use|created if the Noah Mosaic land-surface model was run in WRF|
|SOI_CRO|I/O API|Time-varying soil properties in each soil layer|created if a land-surface model was run in WRF or FV3-GFS|
|mcip.nc|netCDF|contains both time-independent and time-varying output variables that contain 2-D layers (either only in 2-D or in 3-D, where the third dimension could be atmospheric layers, soil layers, land use categories, mosaic categories, etc.)|required, if IOFORM=2|
|mcip_bdy.nc|netCDF|contains time-independent and time-varying output along the domain perimeter|required, if IOFORM=2| 

The default location of the NACC output files is the `$CMAQ_HOME/data/mcip/$GridName` directory, but it can be changed in the NACC script using the `$OutDir` variable. The names of the NACC output files are generic and do not have any information about the model grid that they are simulating or the time period that is covered. These attributes can be controlled by the NACC script. For example, the name of the grid can be used in the output directory path. In addition, the default naming convention for all NACC output files appends the `APPL` environment variable to the file name to identify files by the time period that is represented by the file. All of the file naming variables for the NACC outputs are set in the run script, and they can be easily tailored to fit each user's application or style.

**Previous Versions of NACC**
