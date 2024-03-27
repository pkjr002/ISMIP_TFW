clear all; close all; clc

%Routine to recompute scalar values of manuscript https://tc.copernicus.org/articles/14/3033/2020/
%Contact: Helene Seroussi helene.seroussi@jpl.nasa.gov

step=13;
% error('Add the correct path where the model outputs can be found below on the line below')
datapath='/Volumes/Geek_Boi-Ed/ISMIP6'; %Change path for ISMIP6 model outputs

yts=365.25*24*3600;
% model_list={'AWI_PISM1','DOE_MALI','ILTS_PIK_SICOPOLIS','IMAU_IMAUICE1','IMAU_IMAUICE2','JPL1_ISSM','LSCE_GRISLI','NCAR_CISM','PIK_PISM1','PIK_PISM2','UCIJPL_ISSM','ULB_FETISH32','ULB_FETISH16','UTAS_ElmerIce','VUB_AISMPALEO','VUW_PISM'}; 
model_list={'AWI_PISM1','DOE_MALI','ILTS_PIK_SICOPOLIS1','IMAU_IMAUICE1','IMAU_IMAUICE2','JPL1_ISSM','LSCE_GRISLI','NCAR_CISM','PIK_PISM1','PIK_PISM2','UCIJPL_ISSM','ULB_FETISH32','ULB_FETISH16','UTAS_ElmerIce','VUB_AISMPALEO','VUW_PISM'}; 

model_list2=model_list;
for i=1:numel(model_list2)
	   model_list2{i} = strrep(model_list2{i},'_','\_');
end



% experiments_list={...
% 	'hist_open','hist_std','ctrl_proj_open','ctrl_proj_std',...
% 	'exp01','exp02','exp03','exp04','exp05','exp06','exp07','exp08','exp09','exp10','exp11','exp12','exp13',...
% 	'expA1','expA2','expA3','expA4','expA5','expA6','expA7','expA8',...
% };

experiments_list={...
	'hist_std','ctrl_proj_std',...
	'exp01','exp02','exp03','exp04','exp05','exp06','exp07','exp08','exp09','exp10','exp11','exp12','exp13',...
	'expA1','expA2','expA3','expA4','expA5','expA6','expA7','expA8',...
};


colors = distinguishable_colors(length(model_list));

%Load some datasets needed for all steps
scale_file =['./Data/af2_el_ismip6_ant_01.nc']; %File to rescale distorsion due to polar stereographic projection
scalefac   = double(ncread(scale_file,'af2'));
sectors_32km=ncread(['./Data/sectors_32km.nc'],'sectors'); %18 sectors at 32 km resolution
regions_32km=ncread(['./Data/sectors_32km.nc'],'regions'); %3 regions (West, East and Peninsula)
sectors_16km=ncread(['./Data/sectors_16km.nc'],'sectors'); %18 sectors at 16 km resolution
regions_16km=ncread(['./Data/sectors_16km.nc'],'regions'); %3 regions (West, East and Peninsula)
sectors_8km=ncread(['./Data/sectors_8km.nc'],'sectors'); %18 sectors at 8 km resolution 
regions_8km=ncread(['./Data/sectors_8km.nc'],'regions'); %3 regions (West, East and Peninsula)
sectors_4km=ncread(['./Data/sectors_4km.nc'],'sectors'); %18 sectors at 4 km resolution
regions_4km=ncread(['./Data/sectors_4km.nc'],'regions'); %3 regions (West, East and Peninsula)
numsectors=max(sectors_4km(:));
numregions=max(regions_4km(:));


if step==13, % {{{Compute_ivaf_exp_basins

	for iexp=1:length(experiments_list),
			expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			expename=experiments_list{iexp};
			time_vector=[];
			ivaf_vector=[];
			for isector=1:numsectors,
				eval(['ivaf_vector_sector' num2str(isector) '=[];'])
			end
			for iregion=1:numregions,
				eval(['ivaf_vector_region' num2str(iregion) '=[];'])
            end
            
            disp('==>this1')

			modelname=model_list{imodel};
			specificsv1;
			scalefac_model=scalefac(1:resolution:end,1:resolution:end);
			path=['' datapath '/' group '/' simul ''];
			eval(['exp_name=' expename '_name;'])
% 			exp_directory=[path '/' exp_name ];
            exp_directory=[path '/' expename]; 
            
            disp('==>this2')
            %disp(path)
            
			eval(['isexp=is' expename ';'])
			if isexp,
				if exist(exp_directory,'dir')==0,
					error(['directory ' exp_directory ' not found']);
				end
                disp('==>here')
				field='lithk';
				if is_lithk==1,
					expthickness_file=[exp_directory '/lithk_AIS_' group '_' simul '_' expename '.nc'];
					expbed_file=[exp_directory '/topg_AIS_' group '_' simul '_' expename '.nc'];
					expmask_file=[exp_directory '/sftgif_AIS_' group '_' simul '_' expename '.nc'];
					expgroundmask_file=[exp_directory '/sftgrf_AIS_' group '_' simul '_' expename '.nc'];
					thicknessdata = double(ncread(expthickness_file,'lithk')); 
					beddata = double(ncread(expbed_file,'topg')); 
					mask = double(ncread(expmask_file,'sftgif')); 
					groundmask = double(ncread(expgroundmask_file,'sftgrf')); 
					if(max(mask(:))>1 | min(mask(:))<0), disp(['mask should be between 0 and 1 in model ' modelname '']); end
					pos=find(mask>1);
					mask(pos)=1;
					pos=find(groundmask>1);
					groundmask(pos)=1;
					time = double(ncread(expthickness_file,'time'))/yearday_model;  
					if strcmpi(modelname,'VUW_PISM'),
						disp(['time should be in days in model ' modelname '']);
						time=time/(3600*24);
					end
					time=round(time,1);
					if strcmpi(expename,'ctrl') | strcmpi(expename,'asmb') |strcmpi(expename,'abmb') | strcmpi(expename,'hist_std') | strcmpi(expename,'ctrl_open') | strcmpi(expename,'hist_open'),
						start_model_year=initial_model_year;
					else
						start_model_year=expe_model_year;
					end
					time_vector(1:length(time))=start_model_year+time;
					if strcmpi(expename,'hist_std') | strcmpi(expename,'hist_open'),
						pos_time=find(time_vector<=2015);
					else
						pos_time=find(time_vector>=2016 & time_vector<=2101);
					end
					time_vector=time_vector(pos_time);
					for itime=1:length(time_vector),
						thickness_i=thicknessdata(:,:,itime);
						bed_i=beddata(:,:,itime);
						mask_i=mask(:,:,itime);
						groundmask_i=groundmask(:,:,itime);
						pos=find(mask_i==0);
						thickness_i(pos)=0;
						bed_i(pos)=0;
						groundmask_i(pos)=0;
						posnan=find(isnan(mask_i));
						thickness_i(posnan)=0;
						bed_i(posnan)=0;
						mask_i(posnan)=0;
						groundmask_i(posnan)=0;
						posnan=find(isnan(thickness_i));
						thickness_i(posnan)=0;
						bed_i(posnan)=0;
						mask_i(posnan)=0;
						shelfmask_i(posnan)=0;
						posnan=find(isnan(shelfmask_i));
						thickness_i(posnan)=0;
						bed_i(posnan)=0;
						mask_i(posnan)=0;
						groundmask_i(posnan)=0;
						posnan=find(isnan(bed_i));
						thickness_i(posnan)=0;
						bed_i(posnan)=0;
						mask_i(posnan)=0;
						groundmask_i(posnan)=0;
						hf_i=thickness_i+ocean_density/ice_density*min(bed_i,0);
						volaf=sum(hf_i(:).*mask_i(:).*groundmask_i(:).*scalefac_model(:))*(resolution*1000)^2;
						ivaf_total=volaf; %in m^3
						ivaf_vector(itime)=ivaf_total;
						for isector=1:numsectors,
							eval(['sectors=sectors_' num2str(resolution) 'km;'])
							pos_sector=find(sectors==isector);
							ivaf_total_sector=sum((thickness_i(pos_sector)+ocean_density/ice_density*min(bed_i(pos_sector),0)).*groundmask_i(pos_sector).*mask_i(pos_sector).*scalefac_model(pos_sector))*(resolution*1000)^2; %in m^3
							eval(['ivaf_vector_sector' num2str(isector) '(itime)=ivaf_total_sector;'])
						end
						for iregion=1:numregions,
							eval(['regions=regions_' num2str(resolution) 'km;'])
							pos_region=find(regions==iregion);
							ivaf_total_region=sum((thickness_i(pos_region)+ocean_density/ice_density*min(bed_i(pos_region),0)).*groundmask_i(pos_region).*mask_i(pos_region).*scalefac_model(pos_region))*(resolution*1000)^2; %in m^3
							eval(['ivaf_vector_region' num2str(iregion) '(itime)=ivaf_total_region;'])
						end
					end %end of time

					if ~strcmpi(expename,'hist_std') & ~strcmpi(expename,'hist_open'),
						if length(time)<85, error(['run is too short: ' int2str(length(time)) ' years']); end
						if (size(thicknessdata,3)<85 | size(thicknessdata,3)>122), error('field has the wrong size'); end
					end
					if length(time)~=size(thicknessdata,3), error(['lenght or time and data are not consistent in model ' modelname  ' for experiment ' exp_directory]); end
				end %end of lithk

				modelname
				size(start_model_year+time)
				size(time_vector)
				min(time_vector)
				max(time_vector)

				if exist(['ComputedScalarsPaper/' group '/' simul '/' expename ''],'dir')==0,
					eval(['mkdir ComputedScalarsPaper/' group '/' simul '/' expename '']);
				end
				expivaf_file=['ComputedScalarsPaper/' group '/' simul '/' expename '/computed_ivaf_AIS_' group '_' simul '_' expename '.nc'];
				expivaf_file
				status=WriteNetCDFComputedOutputs(expivaf_file,'ivaf','ice volume above floatation','m^3',time_vector,ivaf_vector,...
					ivaf_vector_sector1,ivaf_vector_sector2,ivaf_vector_sector3,ivaf_vector_sector4,ivaf_vector_sector5,...
					ivaf_vector_sector6,ivaf_vector_sector7,ivaf_vector_sector8,ivaf_vector_sector9,ivaf_vector_sector10,...
					ivaf_vector_sector11,ivaf_vector_sector12,ivaf_vector_sector13,ivaf_vector_sector14,ivaf_vector_sector15,...
					ivaf_vector_sector16,ivaf_vector_sector17,ivaf_vector_sector18,...
					ivaf_vector_region1,ivaf_vector_region2,ivaf_vector_region3,...
					ice_density,ocean_density);

			end %end of isexp
		end %end of model
	end %end of experiment

end %}}}