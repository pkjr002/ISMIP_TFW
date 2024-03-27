This archive provides the scripts and routines used as part of the publication "ISMIP6 Antarctica: a multi-model ensemble of the Antarctic ice sheet evolution over the 21st century", published in The Cryosphere, https://tc.copernicus.org/articles/14/3033/2020/

Seroussi, H., Nowicki, S., Payne, A. J., Goelzer, H., Lipscomb, W. H., Abe-Ouchi, A., Agosta, C., Albrecht, T., Asay-Davis, X., Barthel, A., Calov, R., Cullather, R., Dumas, C., Galton-Fenzi, B. K., Gladstone, R., Golledge, N. R., Gregory, J. M., Greve, R., Hattermann, T., Hoffman, M. J., Humbert, A., Huybrechts, P., Jourdain, N. C., Kleiner, T., Larour, E., Leguy, G. R., Lowry, D. P., Little, C. M., Morlighem, M., Pattyn, F., Pelle, T., Price, S. F., Quiquet, A., Reese, R., Schlegel, N.-J., Shepherd, A., Simon, E., Smith, R. S., Straneo, F., Sun, S., Trusel, L. D., Van Breedam, J., van de Wal, R. S. W., Winkelmann, R., Zhao, C., Zhang, T., and Zwinger, T.: ISMIP6 Antarctica: a multi-model ensemble of the Antarctic ice sheet evolution over the 21st century, The Cryosphere, 14, 3033–3070, https://doi.org/10.5194/tc-14-3033-2020, 2020.

Contact: Helene Seroussi, Helene.seroussi@jpl.nasa.gov

Further information on ISMIP6 and ISMIP6 Antarctica Projections can be found here:
http://www.climate-cryosphere.org/activities/targeted/ismip6
http://www.climate-cryosphere.org/wiki/index.php?title=ISMIP6-Projections-Antarctica

Users should cite the original publication when using all or part of the data. 
In order to document CMIP6's scientific impact and enable ongoing support of CMIP, users are also obligated to acknowledge CMIP6, ISMIP6 and the participating modeling groups.


Archive overview
-----------------------------------------------
README.txt - this information
figures_paper.m - Matlab script to reproduce figures in Seroussi et al. (2020)
scalars_paper.m - Matlab script to compute scalars from two-dimensional fields submitted by ice sheet models
specifics.m - Matlab script including parameters from ice flow models necessary but not included in the output files
WriteNetCDFComputedOutputs.m - Matlab script to write NetCDF files of scalar outputs

Data/ - Directory with datasets needed to compute scalars and reproduce figures
Data/af2_el_ismip6_ant_01.nc - 1 km grid containing area distorsion
Data/sectors_4km.nc - Antarctic regions and sectors at 4 km resolution
Data/sectors_8km.nc - Antarctic regions and sectors at 8 km resolution
Data/sectors_16km.nc - Antarctic regions and sectors at 16 km resolution
Data/sectors_32km.nc - Antarctic regions and sectors at 32 km resolution
Data/sectors_8km_iceonly.nc - Antarctic regions and sectors at 8 km resolution, present-day ice covered areas only

Some additional files might be needed to reproduce the figures:
distinguishable_colors - Matlab file to create colors for ice sheet models: https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors
hatchfill2 - Matlab file to add hatches on collapsed ice shelves: https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/53593/versions/10/previews/hatchfill2.m/index.html
thickness_grid - ice thickness on 8 km standard ISMIP6 grid used to compare with modeled thickness,
for example based on BedMachineAntarctica https://nsidc.org/data/NSIDC-0756
velocity_grid - ice velocity on 8 km standard ISMIP6 grid used to compare with modeled velocity, for
example based on MEaSURE dataset: https://nsidc.org/data/NSIDC-0754/versions/1 
ISMIP6 datasets from original submission and regridded at 8 km to create scalars and some 2d figures
- these can be found on the ISMIP6 archive, please contact ismip6@gmail.com for instructions to get
access

------------------------------------------------

Data usage notice:
If you use any of these results, please acknowledge the work of the people involved in producing them. Acknowledgements should have language similar to the below.

"We thank the Climate and Cryosphere (CliC) effort, which provided support for ISMIP6 through sponsoring of workshops, hosting the ISMIP6 website and wiki, and promoted ISMIP6. We acknowledge the World Climate Research Programme, which, through it's Working Group on Coupled Modelling, coordinated and promoted CMIP5 and CMIP6. We thank the climate modeling groups for producing and making available their model output, the Earth System Grid Federation (ESGF) for archiving the CMIP data and providing access, the University at Buffalo for ISMIP6 data distribution and upload, and the multiple funding agencies who support CMIP5 and CMIP6 and ESGF. We thank the ISMIP6 steering committee, the ISMIP6 model selection group and ISMIP6 dataset preparation group for their continuous engagement in defining ISMIP6."

You should also refer to and cite the following papers:

Seroussi, H., Nowicki, S., Payne, A. J., Goelzer, H., Lipscomb, W. H., Abe-Ouchi, A., Agosta, C., Albrecht, T., Asay-Davis, X., Barthel, A., Calov, R., Cullather, R., Dumas, C., Galton-Fenzi, B. K., Gladstone, R., Golledge, N. R., Gregory, J. M., Greve, R., Hattermann, T., Hoffman, M. J., Humbert, A., Huybrechts, P., Jourdain, N. C., Kleiner, T., Larour, E., Leguy, G. R., Lowry, D. P., Little, C. M., Morlighem, M., Pattyn, F., Pelle, T., Price, S. F., Quiquet, A., Reese, R., Schlegel, N.-J., Shepherd, A., Simon, E., Smith, R. S., Straneo, F., Sun, S., Trusel, L. D., Van Breedam, J., van de Wal, R. S. W., Winkelmann, R., Zhao, C., Zhang, T., and Zwinger, T.: ISMIP6 Antarctica: a multi-model ensemble of the Antarctic ice sheet evolution over the 21st century, The Cryosphere, 14, 3033–3070, https://doi.org/10.5194/tc-14-3033-2020, 2020.

Nowicki, S., Goelzer, H., Seroussi, H., Payne, A. J., Lipscomb, W. H., Abe-Ouchi, A., Agosta, C., Alexander, P., Asay-Davis, X. S., Barthel, A., Bracegirdle, T. J., Cullather, R., Felikson, D., Fettweis, X., Gregory, J. M., Hattermann, T., Jourdain, N. C., Kuipers Munneke, P., Larour, E., Little, C. M., Morlighem, M., Nias, I., Shepherd, A., Simon, E., Slater, D., Smith, R. S., Straneo, F., Trusel, L. D., van den Broeke, M. R., and van de Wal, R.: Experimental protocol for sea level projections from ISMIP6 stand-alone ice sheet models, The Cryosphere, 14, 2331–2368, https://doi.org/10.5194/tc-14-2331-2020, 2020.

