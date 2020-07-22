# ROMSForcingCreator
Scripts for making ROMS forcing files in Matlab. The repository is split into sections for making the lateral boundary conditions (`lbc`), surface boundary conditions (`sbc`), tides (`tds`),  the grid file (`grid`) and river sources (`rvr`).

## Install and use:
1. Download repository.
2. Edit the `make_*.m` file for each forcing type, e.g. `make_sbc.m` to make surface boudnary conditions.
3. Run each `make_*.m` file in Matlab.

More thorough details are provided below:
### Lateral boundary conditions
Edit the file `make_lbc.m`. First set the gridname with `grdname`; the output filename with `bryname`; the timings of the forcing product (e.g. `MinYear=1992; MaxYear=2015;` for interannual forcing; `MinYear=2006; MaxYear=2006;` for climatology forcing); the boundaries of the box to extract data from ECCO2 (e.g. `ECCObounds = [410 525 80 125]`); and the `RunName` (for metadata).

You will also need to match the vertical stretching settings which you chose in the `*.in` file.

The `DataProduct` selects the source data for making the boundary conditions; match the example strings shown in `make_lbc.m`. The `ForcingType` chooses the type of forcing, e.g. climatology, interannual forcing, single year, and is chosen with a number matching the examples in `make_lbc.m`.
### Surface boundary conditions

### Grid

### River sources
Requirements: Antarctic Mapping Tools, available at [Matlab File Exchange](https://au.mathworks.com/matlabcentral/fileexchange/47638-antarctic-mapping-tools)

# Authors:
David E. Gwyther
Eva Cougnon
Ben Galton-Fenzi
