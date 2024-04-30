function status=WriteNetCDFComputedOutputs(file_name,variablename,variablelongname,variableunit,time_data,lim_vector,lim_vector_sector1,lim_vector_sector2,lim_vector_sector3,lim_vector_sector4,lim_vector_sector5,lim_vector_sector6,lim_vector_sector7,lim_vector_sector8,lim_vector_sector9,lim_vector_sector10,lim_vector_sector11,lim_vector_sector12,lim_vector_sector13,lim_vector_sector14,lim_vector_sector15,lim_vector_sector16,lim_vector_sector17,lim_vector_sector18,lim_vector_region1,lim_vector_region2,lim_vector_region3,ice_density,ocean_density);
%Routine to write NetCDF scalar values for manuscript https://tc.copernicus.org/articles/14/3033/2020/
%Contact: Helene Seroussi helene.seroussi@jpl.nasa.gov

%First a few checks
if min(time_data)<1849, error(['time should be at least 1849, not ' int2str(min(time_data)) '']); end
if max(time_data)>2116, error(['time should be no more that 2016, not ' int2str(max(time_data)) '']); end
if any(isnan(lim_vector)), error('NaN in total vector'); end
if any(isnan(lim_vector_region1)), error('NaN in vector region 1'); end
if any(isnan(lim_vector_region2)), error('NaN in vector region 2'); end
if any(isnan(lim_vector_region3)), error('NaN in vector region 3'); end
if any(isnan(lim_vector_sector1)), error('NaN in vector sector 1'); end
if any(isnan(lim_vector_sector2)), error('NaN in vector sector 2'); end
if any(isnan(lim_vector_sector3)), error('NaN in vector sector 3'); end
if any(isnan(lim_vector_sector4)), error('NaN in vector sector 4'); end
if any(isnan(lim_vector_sector5)), error('NaN in vector sector 5'); end
if any(isnan(lim_vector_sector6)), error('NaN in vector sector 6'); end
if any(isnan(lim_vector_sector7)), error('NaN in vector sector 7'); end
if any(isnan(lim_vector_sector8)), error('NaN in vector sector 8'); end
if any(isnan(lim_vector_sector9)), error('NaN in vector sector 9'); end
if any(isnan(lim_vector_sector10)), error('NaN in vector sector 10'); end
if any(isnan(lim_vector_sector11)), error('NaN in vector sector 11'); end
if any(isnan(lim_vector_sector12)), error('NaN in vector sector 12'); end
if any(isnan(lim_vector_sector13)), error('NaN in vector sector 13'); end
if any(isnan(lim_vector_sector14)), error('NaN in vector sector 14'); end
if any(isnan(lim_vector_sector15)), error('NaN in vector sector 15'); end
if any(isnan(lim_vector_sector16)), error('NaN in vector sector 16'); end
if any(isnan(lim_vector_sector17)), error('NaN in vector sector 17'); end
if any(isnan(lim_vector_sector18)), error('NaN in vector sector 18'); end

%Write netcdf
mode = netcdf.getConstant('NETCDF4');
mode = bitor(mode,netcdf.getConstant('CLASSIC_MODEL'));
ncid=netcdf.create([file_name],mode);

%General attributes
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Project','ISMIP6 Antarctica Projections');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Description','Rescaled scalars from 2d fields');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Grid','Datum = WGS84, earth_radius = 6378137., earth_eccentricity = 0.081819190842621, falseeasting = -3040000., falsenorthing = -3040000., standard_parallel = -71., central_meridien = 0, EPSG=3031');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'proj','+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'proj4','+init=epsg:3031');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Date',date());

%Define dimensions
time_id  = netcdf.defDim(ncid,'time',length(time_data));
unit_id  = netcdf.defDim(ncid,'unit',1);

%Define variables
time_var_id= netcdf.defVar(ncid,'time','NC_FLOAT',[time_id]);
netcdf.putAtt(ncid,time_var_id,'standard_name','time');
netcdf.putAtt(ncid,time_var_id,'long_name','time in years');
netcdf.putAtt(ncid,time_var_id,'units', 'years since 0000-00-00 00:00:00');
netcdf.putAtt(ncid,time_var_id,'calendar', '365 days');
netcdf.putAtt(ncid,time_var_id,'axis', 'T');

totalnumber=[variablename '_antarctica'];
variablenumber=[variablename ' Antarctica '];
variablelongnumber=[variablelongname ' Antarctica '];
eval(['total_var_id= netcdf.defVar(ncid,variablename,''NC_FLOAT'',[time_id]);'])
eval(['netcdf.putAtt(ncid,total_var_id,''standard_name'',variablename);'])
eval(['netcdf.putAtt(ncid,total_var_id,''long_name'',variablelongnumber);'])
eval(['netcdf.putAtt(ncid,total_var_id,''units'',variableunit);'])

for indsector=1:18,
sectornumber=[variablename '_sector_' int2str(indsector) ''];
variablenumber=[variablename ' sector ' int2str(indsector) ''];
variablelongnumber=[variablelongname ' sector ' int2str(indsector) ''];
eval(['sector' int2str(indsector) '_var_id= netcdf.defVar(ncid,sectornumber,''NC_FLOAT'',[time_id]);'])
eval(['netcdf.putAtt(ncid,sector' int2str(indsector) '_var_id,''standard_name'',variablenumber);'])
eval(['netcdf.putAtt(ncid,sector' int2str(indsector) '_var_id,''long_name'',variablelongnumber);'])
eval(['netcdf.putAtt(ncid,sector' int2str(indsector) '_var_id,''units'',variableunit);'])
end

for indregion=1:3,
regionnumber=[variablename '_region_' int2str(indregion) ''];
variablenumber=[variablename ' region ' int2str(indregion) ''];
variablelongnumber=[variablelongname ' region ' int2str(indregion) ''];
eval(['region' int2str(indregion) '_var_id= netcdf.defVar(ncid,regionnumber,''NC_FLOAT'',[time_id]);'])
eval(['netcdf.putAtt(ncid,region' int2str(indregion) '_var_id,''standard_name'',variablenumber);'])
if indregion==1,
variablelongnumber=[variablelongname ' region ' int2str(indregion) ' (West Antarctica)'];
elseif indregion==2,
variablelongnumber=[variablelongname ' region ' int2str(indregion) ' (East Antarctica)'];
elseif indregion==3,
variablelongnumber=[variablelongname ' region ' int2str(indregion) ' (Peninsula)'];
else error('');
end
eval(['netcdf.putAtt(ncid,region' int2str(indregion) '_var_id,''long_name'',variablelongnumber);'])
eval(['netcdf.putAtt(ncid,region' int2str(indregion) '_var_id,''units'',variableunit);'])
end

rhoi_var_id= netcdf.defVar(ncid,'rhoi','NC_FLOAT',[unit_id]);
netcdf.putAtt(ncid,rhoi_var_id,'standard_name','rhoi');
netcdf.putAtt(ncid,rhoi_var_id,'long_name','model specific ice density');
netcdf.putAtt(ncid,rhoi_var_id,'units', 'kg/m^3');

rhow_var_id= netcdf.defVar(ncid,'rhow','NC_FLOAT',[unit_id]);
netcdf.putAtt(ncid,rhow_var_id,'standard_name','rhow');
netcdf.putAtt(ncid,rhow_var_id,'long_name','model specific ocean density');
netcdf.putAtt(ncid,rhow_var_id,'units', 'kg/m^3');

netcdf.endDef(ncid);

%Write variables
netcdf.putVar(ncid,time_var_id,time_data');
netcdf.putVar(ncid,total_var_id,lim_vector');
for indsector=1:18,
eval(['netcdf.putVar(ncid,sector' int2str(indsector) '_var_id,transpose(lim_vector_sector' int2str(indsector) '));'])
end
for indregion=1:3,
eval(['netcdf.putVar(ncid,region' int2str(indregion) '_var_id,transpose(lim_vector_region' int2str(indregion) '));'])
end
netcdf.putVar(ncid,rhoi_var_id,ice_density);
netcdf.putVar(ncid,rhow_var_id,ocean_density);

%Close netcdf
netcdf.close(ncid)
status=1;
