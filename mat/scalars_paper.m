%Routine to recompute scalar values of manuscript https://tc.copernicus.org/articles/14/3033/2020/
%Contact: Helene Seroussi helene.seroussi@jpl.nasa.gov

step=1;
% error('Add the correct path where the model outputs can be found below on the line below')
datapath='/Volumes/Geek_Boi-Ed/ISMIP6/'; %Change path for ISMIP6 model outputs

yts=365.25*24*3600;
model_list={'AWI_PISM1','DOE_MALI','ILTS_PIK_SICOPOLIS','IMAU_IMAUICE1','IMAU_IMAUICE2','JPL1_ISSM','LSCE_GRISLI','NCAR_CISM','PIK_PISM1','PIK_PISM2','UCIJPL_ISSM','ULB_FETISH32','ULB_FETISH16','UTAS_ElmerIce','VUB_AISMPALEO','VUW_PISM'}; 

model_list2=model_list;
for i=1:numel(model_list2)
	   model_list2{i} = strrep(model_list2{i},'_','\_');
end

experiments_list={...
	'hist_open','hist_std','ctrl_proj_open','ctrl_proj_std',...
	'exp01','exp02','exp03','exp04','exp05','exp06','exp07','exp08','exp09','exp10','exp11','exp12','exp13',...
	'expA1','expA2','expA3','expA4','expA5','expA6','expA7','expA8',...
};

colors = distinguishable_colors(length(model_list));

%Load some datasets needed for all steps
scale_file =['' datapath '/Data/af2_el_ismip6_ant_01.nc']; %File to rescale distorsion due to polar stereographic projection
scalefac   = double(ncread(scale_file,'af2'));
sectors_32km=ncread(['' datapath '/Data/sectors_32km.nc'],'sectors'); %18 sectors at 32 km resolution
regions_32km=ncread(['' datapath '/Data/sectors_32km.nc'],'regions'); %3 regions (West, East and Peninsula)
sectors_16km=ncread(['' datapath '/Data/sectors_16km.nc'],'sectors'); %18 sectors at 16 km resolution
regions_16km=ncread(['' datapath '/Data/sectors_16km.nc'],'regions'); %3 regions (West, East and Peninsula)
sectors_8km=ncread(['' datapath '/Data/sectors_8km.nc'],'sectors'); %18 sectors at 8 km resolution 
regions_8km=ncread(['' datapath '/Data/sectors_8km.nc'],'regions'); %3 regions (West, East and Peninsula)
sectors_4km=ncread(['' datapath '/Data/sectors_4km.nc'],'sectors'); %18 sectors at 4 km resolution
regions_4km=ncread(['' datapath '/Data/sectors_4km.nc'],'regions'); %3 regions (West, East and Peninsula)
numsectors=max(sectors_4km(:));
numregions=max(regions_4km(:));

if step==1, % {{{Compute_icearea_exp_basins

	for iexp=1:length(experiments_list),
			expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			expename=experiments_list{iexp};
			time_vector=[];
			icearea_vector=[];
			for isector=1:numsectors,
				eval(['icearea_vector_sector' num2str(isector) '=[];'])
			end
			for iregion=1:numregions,
				eval(['icearea_vector_region' num2str(iregion) '=[];'])
			end
			modelname=model_list{imodel};
			specifics;
			scalefac_model=scalefac(1:resolution:end,1:resolution:end);
			path=['' datapath '/' group '/' simul ''];
			eval(['exp_name=' expename '_name;'])
			exp_directory=[path '/' exp_name '/'];

			eval(['isexp=is' expename ';'])
			if isexp,
				if exist(exp_directory,'dir')==0,
					error(['directory ' exp_directory ' not found']);
				end

				field='sftgif';
				if is_sftgif==1,
					exp_file=[exp_directory '/' field '_AIS_' group '_' simul '_' expename '.nc'];
					data = double(ncread(exp_file,field)); 
					if(max(data(:))>1 | min(data(:))<0), disp(['data should be between 0 and 1 in model ' modelname '']); end
					pos=find(data>1);
					data(pos)=1;
					time = double(ncread(exp_file,'time'))/yearday_model;  
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
						sftgif_i=data(:,:,itime);
						posnan=find(isnan(sftgif_i)); 
						sftgif_i(posnan)=0; 
						icearea_total=sum(sftgif_i(:).*scalefac_model(:))*(resolution*1000)^2; %in m^2
						icearea_vector(itime)=icearea_total; 
						for isector=1:numsectors,
							eval(['sectors=sectors_' num2str(resolution) 'km;'])
							pos_sector=find(sectors==isector);
							icearea_total_sector=sum(sftgif_i(pos_sector).*scalefac_model(pos_sector))*(resolution*1000)^2/(1); %in m^2 
							eval(['icearea_vector_sector' num2str(isector) '(itime)=icearea_total_sector;'])
						end
						for iregion=1:numregions,
							eval(['regions=regions_' num2str(resolution) 'km;'])
							pos_region=find(regions==iregion);
							icearea_total_region=sum(sftgif_i(pos_region).*scalefac_model(pos_region))*(resolution*1000)^2/(1); %in m^2 
							eval(['icearea_vector_region' num2str(iregion) '(itime)=icearea_total_region;'])
						end
					end %end of time

					if ~strcmpi(expename,'hist_std') & ~strcmpi(expename,'hist_open'),
						if length(time)<85, error(['run is too short: ' int2str(length(time)) ' years']); end
						if (size(data,3)<85 | size(data,3)>122), error('field has the wrong size'); end
					end
					if length(time)~=size(data,3), error(['lenght or time and data are not consistent in model ' modelname  ' for experiment ' exp_directory]); end
				end %end

				modelname
				size(start_model_year+time)
				size(time_vector)
				min(time_vector)
				max(time_vector)

				if exist(['ComputedScalarsPaper/' group '/' simul '/' expename ''],'dir')==0,
					eval(['mkdir ComputedScalarsPaper/' group '/' simul '/' expename '']);
				end
				expicearea_file=['ComputedScalarsPaper/' group '/' simul '/' expename '/computed_icearea_AIS_' group '_' simul '_' expename '.nc'];
				expicearea_file
				status=WriteNetCDFComputedOutputs(expicearea_file,'icearea','ice extent','m^2',time_vector,icearea_vector,...
					icearea_vector_sector1,icearea_vector_sector2,icearea_vector_sector3,icearea_vector_sector4,icearea_vector_sector5,...
					icearea_vector_sector6,icearea_vector_sector7,icearea_vector_sector8,icearea_vector_sector9,icearea_vector_sector10,...
					icearea_vector_sector11,icearea_vector_sector12,icearea_vector_sector13,icearea_vector_sector14,icearea_vector_sector15,...
					icearea_vector_sector16,icearea_vector_sector17,icearea_vector_sector18,...
					icearea_vector_region1,icearea_vector_region2,icearea_vector_region3,ice_density,ocean_density);

			end %end of isexp
		end %end of model
	end %end of experiment
end %}}}
if step==2, % {{{Compute_icearea_exp_basins_init_noinit

	%For the models that do not have a historical run, take the first step of the control to be the historical
	experiments_list={'hist'}; %Just create one hist depending on what models use for the control
	model_list={'DOE_MALI','PIK_PISM2','UTAS_ElmerIce'}; 

	for iexp=1:length(experiments_list),
			expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			expename=experiments_list{iexp};
			time_vector=[];
			icearea_vector=[];
			for isector=1:numsectors,
				eval(['icearea_vector_sector' num2str(isector) '=[];'])
			end
			for iregion=1:numregions,
				eval(['icearea_vector_region' num2str(iregion) '=[];'])
			end
			modelname=model_list{imodel};
			specifics;
			scalefac_model=scalefac(1:resolution:end,1:resolution:end);
			path=['' datapath '/' group '/' simul ''];

			if (strcmpi(modelname,'DOE_MALI') | strcmpi(modelname,'UTAS_ElmerIce')),
				expename='ctrl_proj_std'
			elseif strcmpi(modelname,'PIK_PISM2'),
				expename='ctrl_proj_open'
			end
			eval(['exp_name=' expename '_name;'])
			exp_directory=[path '/' exp_name '/'];
			if exist(exp_directory,'dir')==0,
				error(['directory ' exp_directory ' not found']);
			end

			field='sftgif';
			if is_sftgif==1,
				exp_file=[exp_directory '/' field '_AIS_' group '_' simul '_' expename '.nc'];
				data = double(ncread(exp_file,field)); 
				if(max(data(:))>1 | min(data(:))<0), disp(['data should be between 0 and 1 in model ' modelname '']); end
				pos=find(data>1);
				data(pos)=1;
				time = double(ncread(exp_file,'time'))/yearday_model;  
				time=round(time,1);
				start_model_year=expe_model_year;
				time_vector(1:length(time))=start_model_year+time;
				pos_time=find(time_vector<=2015);
				if size(pos_time)>1, error('should not have several values before 2015'); end
				time_vector=time_vector(pos_time);
				for itime=1:length(time_vector),
					sftgif_i=data(:,:,itime);
					posnan=find(isnan(sftgif_i)); 
					sftgif_i(posnan)=0; 
					icearea_total=sum(sftgif_i(:).*scalefac_model(:))*(resolution*1000)^2; %in m^2
					icearea_vector(itime)=icearea_total; 
					for isector=1:numsectors,
						eval(['sectors=sectors_' num2str(resolution) 'km;'])
						pos_sector=find(sectors==isector);
						icearea_total_sector=sum(sftgif_i(pos_sector).*scalefac_model(pos_sector))*(resolution*1000)^2/(1); %in m^2 
						eval(['icearea_vector_sector' num2str(isector) '(itime)=icearea_total_sector;'])
					end
					for iregion=1:numregions,
						eval(['regions=regions_' num2str(resolution) 'km;'])
						pos_region=find(regions==iregion);
						icearea_total_region=sum(sftgif_i(pos_region).*scalefac_model(pos_region))*(resolution*1000)^2/(1); %in m^2
						eval(['icearea_vector_region' num2str(iregion) '(itime)=icearea_total_region;'])
					end
				end %end of time

				if length(time)~=size(data,3), error(['lenght or time and data are not consistent in model ' modelname  ' for experiment ' exp_directory]); end
			end %end of lithk

			modelname
			size(start_model_year+time)
			size(time_vector)
			min(time_vector)
			max(time_vector)


			if (strcmpi(modelname,'DOE_MALI') | strcmpi(modelname,'UTAS_ElmerIce')),
				if exist(['ComputedScalarsPaper/' group '/' simul '/hist_std'],'dir')==0,
					eval(['mkdir ComputedScalarsPaper/' group '/' simul '/hist_std']);
				end
				expicearea_file=['ComputedScalarsPaper/' group '/' simul '/hist_std/computed_icearea_AIS_' group '_' simul '_hist_std.nc'];
			elseif (strcmpi(modelname,'PIK_PISM2')),
				if exist(['ComputedScalarsPaper/' group '/' simul '/hist_open'],'dir')==0,
					eval(['mkdir ComputedScalarsPaper/' group '/' simul '/hist_open']);
				end
				expicearea_file=['ComputedScalarsPaper/' group '/' simul '/hist_open/computed_icearea_AIS_' group '_' simul '_hist_open.nc'];
			end
			expicearea_file
			status=WriteNetCDFComputedOutputs(expicearea_file,'icearea','ice extent','m^2',time_vector,icearea_vector,...
				icearea_vector_sector1,icearea_vector_sector2,icearea_vector_sector3,icearea_vector_sector4,icearea_vector_sector5,...
				icearea_vector_sector6,icearea_vector_sector7,icearea_vector_sector8,icearea_vector_sector9,icearea_vector_sector10,...
				icearea_vector_sector11,icearea_vector_sector12,icearea_vector_sector13,icearea_vector_sector14,icearea_vector_sector15,...
				icearea_vector_sector16,icearea_vector_sector17,icearea_vector_sector18,...
				icearea_vector_region1,icearea_vector_region2,icearea_vector_region3,ice_density,ocean_density);

		end %end of model
	end %end of experiment

end %}}}
if step==3, % {{{Compute_icearea_exp_basins_minus_control

	for iexp=1:length(experiments_list),
			expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			expename=experiments_list{iexp};
			time_vector=[];
			icearea_vector=[];
			for isector=1:numsectors,
				eval(['icearea_vector_sector' num2str(isector) '=[];'])
			end
			for iregion=1:numregions,
				eval(['icearea_vector_region' num2str(iregion) '=[];'])
			end
			modelname=model_list{imodel};
			specifics;
			scalefac_model=scalefac(1:resolution:end,1:resolution:end);
			path=['' datapath '/' group '/' simul ''];
			eval(['exp_name=' expename '_name;'])
			exp_directory=[path '/' exp_name '/'];

			eval(['isexp=is' expename ';'])
			if strcmpi(expename,'hist_std') | strcmpi(expename,'hist_open') | strcmpi(expename,'ctrl_proj_std') | strcmpi(expename,'ctrl_proj_open'),
				%do nothing, we cannot substract the control for the historical or control runs
			elseif isexp,
				%Prepare the ctrl values
				if strcmpi(expename,'exp01') | strcmpi(expename,'exp02') | strcmpi(expename,'exp03') | strcmpi(expename,'exp04') | strcmpi(expename,'exp11') | strcmpi(expename,'expA1') | strcmpi(expename,'expA2') | strcmpi(expename,'expA3') | strcmpi(expename,'expA4'),
					ctrl_file=['ComputedScalarsPaper/' group '/' simul '/ctrl_proj_open/computed_icearea_AIS_' group '_' simul '_ctrl_proj_open.nc'];
				elseif strcmpi(expename,'exp05') | strcmpi(expename,'exp06') | strcmpi(expename,'exp07') | strcmpi(expename,'exp08') | strcmpi(expename,'exp09') | strcmpi(expename,'exp10') | strcmpi(expename,'exp12') | strcmpi(expename,'exp13') | strcmpi(expename,'expA5') | strcmpi(expename,'expA6') | strcmpi(expename,'expA7') | strcmpi(expename,'expA8'),
					ctrl_file=['ComputedScalarsPaper/' group '/' simul '/ctrl_proj_std/computed_icearea_AIS_' group '_' simul '_ctrl_proj_std.nc'];
				else error('experiment not supported yet');
				end

				if exist(exp_directory,'dir')==0,
					error(['directory ' exp_directory ' not found']);
				end

				field='sftgif';
				if is_sftgif==1,
					exp_file=[exp_directory '/' field '_AIS_' group '_' simul '_' expename '.nc'];
					data = double(ncread(exp_file,field)); 
					if(max(data(:))>1 | min(data(:))<0), disp(['data should be between 0 and 1 in model ' modelname '']); end
					pos=find(data>1);
					data(pos)=1;
					time = double(ncread(exp_file,'time'))/yearday_model;  
					if strcmpi(modelname,'VUW_PISM'),
						disp(['time should be in days in model ' modelname '']);
						time=time/(3600*24);
					end
					time=round(time,1);
					start_model_year=expe_model_year;
					time_vector(1:length(time))=start_model_year+time;
					if strcmpi(modelname,'ULB_FETISH32') | strcmpi(modelname,'ULB_FETISH16'),
						pos_time=find(time_vector>=2016 & time_vector<=2100);
					else
						pos_time=find(time_vector>=2016 & time_vector<=2101);
					end
					time_vector=time_vector(pos_time);
					for itime=1:length(time_vector),
						sftgif_i=data(:,:,itime);
						posnan=find(isnan(sftgif_i)); 
						sftgif_i(posnan)=0; 
						icearea_total=sum(sftgif_i(:).*scalefac_model(:))*(resolution*1000)^2; %in m^2
						icearea_vector(itime)=icearea_total; 
						for isector=1:numsectors,
							eval(['sectors=sectors_' num2str(resolution) 'km;'])
							pos_sector=find(sectors==isector);
							icearea_total_sector=sum(sftgif_i(pos_sector).*scalefac_model(pos_sector))*(resolution*1000)^2; %in m^2 
							eval(['icearea_total_sector_ctrlproj= double(ncread(ctrl_file,''icearea_sector_' num2str(isector) '''));'])
							eval(['icearea_vector_sector' num2str(isector) '(itime)=icearea_total_sector-icearea_total_sector_ctrlproj(itime);'])
						end
						for iregion=1:numregions,
							eval(['regions=regions_' num2str(resolution) 'km;'])
							pos_region=find(regions==iregion);
							icearea_total_region=sum(sftgif_i(pos_region).*scalefac_model(pos_region))*(resolution*1000)^2; %in m^2 
							eval(['icearea_total_region_ctrlproj= double(ncread(ctrl_file,''icearea_region_' num2str(iregion) '''));'])
							eval(['icearea_vector_region' num2str(iregion) '(itime)=icearea_total_region-icearea_total_region_ctrlproj(itime);'])
						end
					end %end of time
					icearea_total_ctrlproj= double(ncread(ctrl_file,'icearea'));
					if strcmpi(modelname,'DOE_MALI'),
						icearea_vector=icearea_vector-icearea_total_ctrlproj(1:end-1)';
					else
						icearea_vector=icearea_vector-icearea_total_ctrlproj';
					end

					if length(time)<85, error(['run is too short: ' int2str(length(time)) ' years']); end
					if (size(data,3)<85 | size(data,3)>122), error('field has the wrong size'); end
					if length(time)~=size(data,3), error(['lenght or time and data are not consistent in model ' modelname  ' for experiment ' exp_directory]); end
				end %end of lithk

				modelname
				size(start_model_year+time)
				size(time_vector)
				min(time_vector)
				max(time_vector)

				expicearea_file=['ComputedScalarsPaper/' group '/' simul '/' expename '/computed_icearea_minus_ctrl_proj_AIS_' group '_' simul '_' expename '.nc'];
				expicearea_file
				status=WriteNetCDFComputedOutputs(expicearea_file,'icearea','ice extent','m^2',time_vector,icearea_vector,...
					icearea_vector_sector1,icearea_vector_sector2,icearea_vector_sector3,icearea_vector_sector4,icearea_vector_sector5,...
					icearea_vector_sector6,icearea_vector_sector7,icearea_vector_sector8,icearea_vector_sector9,icearea_vector_sector10,...
					icearea_vector_sector11,icearea_vector_sector12,icearea_vector_sector13,icearea_vector_sector14,icearea_vector_sector15,...
					icearea_vector_sector16,icearea_vector_sector17,icearea_vector_sector18,...
					icearea_vector_region1,icearea_vector_region2,icearea_vector_region3,...
					ice_density,ocean_density);

			end %end of isexp
		end %end of model
	end %end of experiment

end %}}}
if step==4, % {{{Compute_iareagr_exp_basins

	for iexp=1:length(experiments_list),
			expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			expename=experiments_list{iexp};
			time_vector=[];
			iareagr_vector=[];
			for isector=1:numsectors,
				eval(['iareagr_vector_sector' num2str(isector) '=[];'])
			end
			for iregion=1:numregions,
				eval(['iareagr_vector_region' num2str(iregion) '=[];'])
			end
			modelname=model_list{imodel};
			specifics;
			scalefac_model=scalefac(1:resolution:end,1:resolution:end);
			path=['' datapath '/' group '/' simul ''];
			eval(['exp_name=' expename '_name;'])
			exp_directory=[path '/' exp_name '/'];

			eval(['isexp=is' expename ';'])
			if isexp,
				if exist(exp_directory,'dir')==0,
					error(['directory ' exp_directory ' not found']);
				end

				field='sftgrf';
				if is_sftgrf==1,
					exp_file=[exp_directory '/' field '_AIS_' group '_' simul '_' expename '.nc'];
					expmask_file=[exp_directory '/sftgif_AIS_' group '_' simul '_' expename '.nc'];
					data = double(ncread(exp_file,field)); 
					mask = double(ncread(expmask_file,'sftgif')); 
					if(max(mask(:))>1 | min(mask(:))<0), disp(['mask should be between 0 and 1 in model ' modelname '']); end
					pos=find(mask>1);
					mask(pos)=1;
					time = double(ncread(exp_file,'time'))/yearday_model;  
					if strcmpi(modelname,'VUW_PISM'),
						disp(['time should be in days in model ' modelname '']);
						time=time/(3600*24);
					end
					time=round(time,1);
					if strcmpi(expename,'ctrl') | strcmpi(expename,'asmb') |strcmpi(expename,'abmb') | strcmpi(expename,'hist') | strcmpi(expename,'ctrl_open') | strcmpi(expename,'hist_open'),
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
						sftgrf_i=data(:,:,itime);
						mask_i=mask(:,:,itime);
						pos=find(mask_i==0);
						sftgrf_i(pos)=0;
						posnan=find(isnan(mask_i));
						sftgrf_i(posnan)=0;
						mask_i(posnan)=0;
						posnan=find(isnan(sftgrf_i)); 
						sftgrf_i(posnan)=0; 
						mask_i(posnan)=0; 
						iareagr_total=sum(sftgrf_i(:).*mask_i(:).*scalefac_model(:))*(resolution*1000)^2; %in m^2
						iareagr_vector(itime)=iareagr_total; 
						for isector=1:numsectors,
							eval(['sectors=sectors_' num2str(resolution) 'km;'])
							pos_sector=find(sectors==isector);
							iareagr_total_sector=sum(sftgrf_i(pos_sector).*mask_i(pos_sector).*scalefac_model(pos_sector))*(resolution*1000)^2; %in m^2
							eval(['iareagr_vector_sector' num2str(isector) '(itime)=iareagr_total_sector;'])
						end
						for iregion=1:numregions,
							eval(['regions=regions_' num2str(resolution) 'km;'])
							pos_region=find(regions==iregion);
							iareagr_total_region=sum(sftgrf_i(pos_region).*mask_i(pos_region).*scalefac_model(pos_region))*(resolution*1000)^2; %in m^2
							eval(['iareagr_vector_region' num2str(iregion) '(itime)=iareagr_total_region;'])
						end
					end %end of time

					if ~strcmpi(expename,'hist_std') & ~strcmpi(expename,'hist_open'),
						if length(time)<85, error(['run is too short: ' int2str(length(time)) ' years']); end
						if (size(data,3)<85 | size(data,3)>122), error('field has the wrong size'); end
					end
					if length(time)~=size(data,3), error(['lenght or time and data are not consistent in model ' modelname  ' for experiment ' exp_directory]); end
				end %end of lithk

				modelname
				size(start_model_year+time)
				size(time_vector)
				min(time_vector)
				max(time_vector)

				if exist(['ComputedScalarsPaper/' group '/' simul '/' expename ''],'dir')==0,
					eval(['mkdir ComputedScalarsPaper/' group '/' simul '/' expename '']);
				end
				expiareagr_file=['ComputedScalarsPaper/' group '/' simul '/' expename '/computed_iareagr_AIS_' group '_' simul '_' expename '.nc'];
				expiareagr_file
				status=WriteNetCDFComputedOutputs(expiareagr_file,'iareagr','grounded ice extent','m^2',time_vector,iareagr_vector,...
					iareagr_vector_sector1,iareagr_vector_sector2,iareagr_vector_sector3,iareagr_vector_sector4,iareagr_vector_sector5,...
					iareagr_vector_sector6,iareagr_vector_sector7,iareagr_vector_sector8,iareagr_vector_sector9,iareagr_vector_sector10,...
					iareagr_vector_sector11,iareagr_vector_sector12,iareagr_vector_sector13,iareagr_vector_sector14,iareagr_vector_sector15,...
					iareagr_vector_sector16,iareagr_vector_sector17,iareagr_vector_sector18,...
					iareagr_vector_region1,iareagr_vector_region2,iareagr_vector_region3,...
					ice_density,ocean_density);

			end %end of isexp
		end %end of model
	end %end of experiment

end %}}}
if step==5, % {{{Compute_iareagr_exp_basins_init_noinit

	%For the models that do not have a historical run, take the first step of the control to be the historical
	experiments_list={'hist'}; %Just create one hist depending on what models use for the control
	model_list={'DOE_MALI','PIK_PISM2','UTAS_ElmerIce'}; 

	for iexp=1:length(experiments_list),
			expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			expename=experiments_list{iexp};
			time_vector=[];
			iareagr_vector=[];
			for isector=1:numsectors,
				eval(['iareagr_vector_sector' num2str(isector) '=[];'])
			end
			for iregion=1:numregions,
				eval(['iareagr_vector_region' num2str(iregion) '=[];'])
			end
			modelname=model_list{imodel};
			specifics;
			scalefac_model=scalefac(1:resolution:end,1:resolution:end);
			path=['' datapath '/' group '/' simul ''];

			if (strcmpi(modelname,'DOE_MALI') | strcmpi(modelname,'UTAS_ElmerIce')),
				expename='ctrl_proj_std'
			elseif strcmpi(modelname,'PIK_PISM2'),
				expename='ctrl_proj_open'
			end
			eval(['exp_name=' expename '_name;'])
			exp_directory=[path '/' exp_name '/'];
			if exist(exp_directory,'dir')==0,
				error(['directory ' exp_directory ' not found']);
			end

			field='sftgrf';
			if is_sftgrf==1,
				exp_file=[exp_directory '/' field '_AIS_' group '_' simul '_' expename '.nc'];
				expmask_file=[exp_directory '/sftgif_AIS_' group '_' simul '_' expename '.nc'];
				data = double(ncread(exp_file,field)); 
				mask = double(ncread(expmask_file,'sftgif')); 
				if(max(mask(:))>1 | min(mask(:))<0), disp(['mask should be between 0 and 1 in model ' modelname '']); end
				pos=find(mask>1);
				mask(pos)=1;
				time = double(ncread(exp_file,'time'))/yearday_model;  
				time=round(time,1);
				start_model_year=expe_model_year;
				time_vector(1:length(time))=start_model_year+time;
				pos_time=find(time_vector<=2015);
				if size(pos_time)>1, error('should not have several values before 2015'); end
				time_vector=time_vector(pos_time);
				for itime=1:length(time_vector),
					sftgrf_i=data(:,:,itime);
					mask_i=mask(:,:,itime);
					pos=find(mask_i==0);
					sftgrf_i(pos)=0;
					posnan=find(isnan(mask_i));
					sftgrf_i(posnan)=0;
					mask_i(posnan)=0;
					posnan=find(isnan(sftgrf_i)); 
					sftgrf_i(posnan)=0; 
					mask_i(posnan)=0; 
					iareagr_total=sum(sftgrf_i(:).*mask_i(:).*scalefac_model(:))*(resolution*1000)^2; %in m^2
					iareagr_vector(itime)=iareagr_total; 
					for isector=1:numsectors,
						eval(['sectors=sectors_' num2str(resolution) 'km;'])
						pos_sector=find(sectors==isector);
						iareagr_total_sector=sum(sftgrf_i(pos_sector).*mask_i(pos_sector).*scalefac_model(pos_sector))*(resolution*1000)^2; %in m^2 
						eval(['iareagr_vector_sector' num2str(isector) '(itime)=iareagr_total_sector;'])
					end
					for iregion=1:numregions,
						eval(['regions=regions_' num2str(resolution) 'km;'])
						pos_region=find(regions==iregion);
						iareagr_total_region=sum(sftgrf_i(pos_region).*mask_i(pos_region).*scalefac_model(pos_region))*(resolution*1000)^2; %in m^2 
						eval(['iareagr_vector_region' num2str(iregion) '(itime)=iareagr_total_region;'])
					end
				end %end of time

				if length(time)~=size(data,3), error(['lenght or time and data are not consistent in model ' modelname  ' for experiment ' exp_directory]); end
			end %end of lithk

			modelname
			size(start_model_year+time)
			size(time_vector)
			min(time_vector)
			max(time_vector)

			if (strcmpi(modelname,'DOE_MALI') | strcmpi(modelname,'UTAS_ElmerIce')),
				if exist(['ComputedScalarsPaper/' group '/' simul '/hist_std'],'dir')==0,
					eval(['mkdir ComputedScalarsPaper/' group '/' simul '/hist_std']);
				end
				expiareagr_file=['ComputedScalarsPaper/' group '/' simul '/hist_std/computed_iareagr_AIS_' group '_' simul '_hist_std.nc'];
			elseif (strcmpi(modelname,'PIK_PISM2')),
				if exist(['ComputedScalarsPaper/' group '/' simul '/hist_open'],'dir')==0,
					eval(['mkdir ComputedScalarsPaper/' group '/' simul '/hist_open']);
				end
				expiareagr_file=['ComputedScalarsPaper/' group '/' simul '/hist_open/computed_iareagr_AIS_' group '_' simul '_hist_open.nc'];
			end
			expiareagr_file
			status=WriteNetCDFComputedOutputs(expiareagr_file,'iareagr','grounded ice extent','m^2',time_vector,iareagr_vector,...
				iareagr_vector_sector1,iareagr_vector_sector2,iareagr_vector_sector3,iareagr_vector_sector4,iareagr_vector_sector5,...
				iareagr_vector_sector6,iareagr_vector_sector7,iareagr_vector_sector8,iareagr_vector_sector9,iareagr_vector_sector10,...
				iareagr_vector_sector11,iareagr_vector_sector12,iareagr_vector_sector13,iareagr_vector_sector14,iareagr_vector_sector15,...
				iareagr_vector_sector16,iareagr_vector_sector17,iareagr_vector_sector18,...
				iareagr_vector_region1,iareagr_vector_region2,iareagr_vector_region3,...
				ice_density,ocean_density);

		end %end of model
	end %end of experiment

end %}}}
if step==6, % {{{Compute_iareagr_exp_basins_minus_control

	for iexp=1:length(experiments_list),
			expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			expename=experiments_list{iexp};
			time_vector=[];
			iareagr_vector=[];
			for isector=1:numsectors,
				eval(['iareagr_vector_sector' num2str(isector) '=[];'])
			end
			for iregion=1:numregions,
				eval(['iareagr_vector_region' num2str(iregion) '=[];'])
			end
			modelname=model_list{imodel};
			specifics;
			scalefac_model=scalefac(1:resolution:end,1:resolution:end);
			path=['' datapath '/' group '/' simul ''];
			eval(['exp_name=' expename '_name;'])
			exp_directory=[path '/' exp_name '/'];

			eval(['isexp=is' expename ';'])
			if strcmpi(expename,'hist_std') | strcmpi(expename,'hist_open') | strcmpi(expename,'ctrl_proj_std') | strcmpi(expename,'ctrl_proj_open'),
				%do nothing, we cannot substract the control for the historical or control runs
			elseif isexp,
				%Prepare the ctrl values
				if strcmpi(expename,'exp01') | strcmpi(expename,'exp02') | strcmpi(expename,'exp03') | strcmpi(expename,'exp04') | strcmpi(expename,'exp11') | strcmpi(expename,'expA1') | strcmpi(expename,'expA2') | strcmpi(expename,'expA3') | strcmpi(expename,'expA4'),
					ctrl_file=['ComputedScalarsPaper/' group '/' simul '/ctrl_proj_open/computed_iareagr_AIS_' group '_' simul '_ctrl_proj_open.nc'];
				elseif strcmpi(expename,'exp05') | strcmpi(expename,'exp06') | strcmpi(expename,'exp07') | strcmpi(expename,'exp08') | strcmpi(expename,'exp09') | strcmpi(expename,'exp10') | strcmpi(expename,'exp12') | strcmpi(expename,'exp13') | strcmpi(expename,'expA5') | strcmpi(expename,'expA6') | strcmpi(expename,'expA7') | strcmpi(expename,'expA8'),
					ctrl_file=['ComputedScalarsPaper/' group '/' simul '/ctrl_proj_std/computed_iareagr_AIS_' group '_' simul '_ctrl_proj_std.nc'];
				else error('experiment not supported yet');
				end

				if exist(exp_directory,'dir')==0,
					error(['directory ' exp_directory ' not found']);
				end

				field='sftgrf';
				if is_sftgrf==1,
					exp_file=[exp_directory '/' field '_AIS_' group '_' simul '_' expename '.nc'];
					expmask_file=[exp_directory '/sftgif_AIS_' group '_' simul '_' expename '.nc'];
					data = double(ncread(exp_file,field)); 
					mask = double(ncread(expmask_file,'sftgif')); 
					if(max(mask(:))>1 | min(mask(:))<0), disp(['mask should be between 0 and 1 in model ' modelname '']); end
					pos=find(mask>1);
					mask(pos)=1;
					time = double(ncread(exp_file,'time'))/yearday_model;  
					if strcmpi(modelname,'VUW_PISM'),
						disp(['time should be in days in model ' modelname '']);
						time=time/(3600*24);
					end
					time=round(time,1);
					start_model_year=expe_model_year;
					time_vector(1:length(time))=start_model_year+time;
					if strcmpi(modelname,'ULB_FETISH32') | strcmpi(modelname,'ULB_FETISH16'),
						pos_time=find(time_vector>=2016 & time_vector<=2100);
					else
						pos_time=find(time_vector>=2016 & time_vector<=2101);
					end
					time_vector=time_vector(pos_time);
					for itime=1:length(time_vector),
						sftgrf_i=data(:,:,itime);
						mask_i=mask(:,:,itime);
						pos=find(mask_i==0);
						sftgrf_i(pos)=0;
						posnan=find(isnan(mask_i));
						sftgrf_i(posnan)=0;
						mask_i(posnan)=0;
						posnan=find(isnan(sftgrf_i)); 
						sftgrf_i(posnan)=0; 
						mask_i(posnan)=0; 
						iareagr_total=sum(sftgrf_i(:).*mask_i(:).*scalefac_model(:))*(resolution*1000)^2; %in m^2
						iareagr_vector(itime)=iareagr_total; 
						for isector=1:numsectors,
							eval(['sectors=sectors_' num2str(resolution) 'km;'])
							pos_sector=find(sectors==isector);
							iareagr_total_sector=sum(sftgrf_i(pos_sector).*mask_i(pos_sector).*scalefac_model(pos_sector))*(resolution*1000)^2; %in m^2
							eval(['iareagr_total_sector_ctrlproj= double(ncread(ctrl_file,''iareagr_sector_' num2str(isector) '''));'])
							eval(['iareagr_vector_sector' num2str(isector) '(itime)=iareagr_total_sector-iareagr_total_sector_ctrlproj(itime);'])
						end
						for iregion=1:numregions,
							eval(['regions=regions_' num2str(resolution) 'km;'])
							pos_region=find(regions==iregion);
							iareagr_total_region=sum(sftgrf_i(pos_region).*mask_i(pos_region).*scalefac_model(pos_region))*(resolution*1000)^2; %in m^2
							eval(['iareagr_total_region_ctrlproj= double(ncread(ctrl_file,''iareagr_region_' num2str(iregion) '''));'])
							eval(['iareagr_vector_region' num2str(iregion) '(itime)=iareagr_total_region-iareagr_total_region_ctrlproj(itime);'])
						end
					end %end of time
					iareagr_total_ctrlproj= double(ncread(ctrl_file,'iareagr'));
					if strcmpi(modelname,'DOE_MALI'),
						iareagr_vector=iareagr_vector-iareagr_total_ctrlproj(1:end-1)';
					else
						iareagr_vector=iareagr_vector-iareagr_total_ctrlproj';
					end

					if length(time)<85, error(['run is too short: ' int2str(length(time)) ' years']); end
					if (size(data,3)<85 | size(data,3)>122), error('field has the wrong size'); end
					if length(time)~=size(data,3), error(['lenght or time and data are not consistent in model ' modelname  ' for experiment ' exp_directory]); end
				end %end of lithk

				modelname
				size(start_model_year+time)
				size(time_vector)
				min(time_vector)
				max(time_vector)

				expiareagr_file=['ComputedScalarsPaper/' group '/' simul '/' expename '/computed_iareagr_minus_ctrl_proj_AIS_' group '_' simul '_' expename '.nc'];
				expiareagr_file
				status=WriteNetCDFComputedOutputs(expiareagr_file,'iareagr','grounded ice extent','m^2',time_vector,iareagr_vector,...
					iareagr_vector_sector1,iareagr_vector_sector2,iareagr_vector_sector3,iareagr_vector_sector4,iareagr_vector_sector5,...
					iareagr_vector_sector6,iareagr_vector_sector7,iareagr_vector_sector8,iareagr_vector_sector9,iareagr_vector_sector10,...
					iareagr_vector_sector11,iareagr_vector_sector12,iareagr_vector_sector13,iareagr_vector_sector14,iareagr_vector_sector15,...
					iareagr_vector_sector16,iareagr_vector_sector17,iareagr_vector_sector18,...
					iareagr_vector_region1,iareagr_vector_region2,iareagr_vector_region3,...
					ice_density,ocean_density);

			end %end of isexp
		end %end of model
	end %end of experiment

end %}}}
if step==7, % {{{Compute_iareafl_exp_basins

	for iexp=1:length(experiments_list),
			expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			expename=experiments_list{iexp};
			time_vector=[];
			iareafl_vector=[];
			for isector=1:numsectors,
				eval(['iareafl_vector_sector' num2str(isector) '=[];'])
			end
			for iregion=1:numregions,
				eval(['iareafl_vector_region' num2str(iregion) '=[];'])
			end
			modelname=model_list{imodel};
			specifics;
			scalefac_model=scalefac(1:resolution:end,1:resolution:end);
			path=['' datapath '/' group '/' simul ''];
			eval(['exp_name=' expename '_name;'])
			exp_directory=[path '/' exp_name '/'];

			eval(['isexp=is' expename ';'])
			if isexp,
				if exist(exp_directory,'dir')==0,
					error(['directory ' exp_directory ' not found']);
				end

				field='sftflf';
				if is_sftflf==1,
					exp_file=[exp_directory '/' field '_AIS_' group '_' simul '_' expename '.nc'];
					expmask_file=[exp_directory '/sftgif_AIS_' group '_' simul '_' expename '.nc'];
					data = double(ncread(exp_file,field)); 
					mask = double(ncread(expmask_file,'sftgif')); 
					if(max(mask(:))>1 | min(mask(:))<0), disp(['mask should be between 0 and 1 in model ' modelname '']); end
					pos=find(mask>1);
					mask(pos)=1;
					time = double(ncread(exp_file,'time'))/yearday_model;  
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
						sftflf_i=data(:,:,itime);
						mask_i=mask(:,:,itime);
						pos=find(mask_i==0);
						sftflf_i(pos)=0;
						posnan=find(isnan(mask_i));
						sftflf_i(posnan)=0;
						mask_i(posnan)=0;
						posnan=find(isnan(sftflf_i)); 
						sftflf_i(posnan)=0; 
						mask_i(posnan)=0; 
						iareafl_total=sum(sftflf_i(:).*mask_i(:).*scalefac_model(:))*(resolution*1000)^2; %in m^2
						iareafl_vector(itime)=iareafl_total; 
						for isector=1:numsectors,
							eval(['sectors=sectors_' num2str(resolution) 'km;'])
							pos_sector=find(sectors==isector);
							iareafl_total_sector=sum(sftflf_i(pos_sector).*mask_i(pos_sector).*scalefac_model(pos_sector))*(resolution*1000)^2; %in m^2
							eval(['iareafl_vector_sector' num2str(isector) '(itime)=iareafl_total_sector;'])
						end
						for iregion=1:numregions,
							eval(['regions=regions_' num2str(resolution) 'km;'])
							pos_region=find(regions==iregion);
							iareafl_total_region=sum(sftflf_i(pos_region).*mask_i(pos_region).*scalefac_model(pos_region))*(resolution*1000)^2; %in m^2
							eval(['iareafl_vector_region' num2str(iregion) '(itime)=iareafl_total_region;'])
						end
					end %end of time

					if ~strcmpi(expename,'hist_std') & ~strcmpi(expename,'hist_open'),
						if length(time)<85, error(['run is too short: ' int2str(length(time)) ' years']); end
						if (size(data,3)<85 | size(data,3)>122), error('field has the wrong size'); end
					end
					if length(time)~=size(data,3), error(['lenght or time and data are not consistent in model ' modelname  ' for experiment ' exp_directory]); end
				end %end of lithk

				modelname
				size(start_model_year+time)
				size(time_vector)
				min(time_vector)
				max(time_vector)

				if exist(['ComputedScalarsPaper/' group '/' simul '/' expename ''],'dir')==0,
					eval(['mkdir ComputedScalarsPaper/' group '/' simul '/' expename '']);
				end
				expiareafl_file=['ComputedScalarsPaper/' group '/' simul '/' expename '/computed_iareafl_AIS_' group '_' simul '_' expename '.nc'];
				expiareafl_file
				status=WriteNetCDFComputedOutputs(expiareafl_file,'iareafl','floating ice extent','m^2',time_vector,iareafl_vector,...
					iareafl_vector_sector1,iareafl_vector_sector2,iareafl_vector_sector3,iareafl_vector_sector4,iareafl_vector_sector5,...
					iareafl_vector_sector6,iareafl_vector_sector7,iareafl_vector_sector8,iareafl_vector_sector9,iareafl_vector_sector10,...
					iareafl_vector_sector11,iareafl_vector_sector12,iareafl_vector_sector13,iareafl_vector_sector14,iareafl_vector_sector15,...
					iareafl_vector_sector16,iareafl_vector_sector17,iareafl_vector_sector18,...
					iareafl_vector_region1,iareafl_vector_region2,iareafl_vector_region3,...
					ice_density,ocean_density);

			end %end of isexp
		end %end of model
	end %end of experiment

end %}}}
if step==8, % {{{Compute_iareafl_exp_basins_init_noinit

	%For the models that do not have a historical run, take the first step of the control to be the historical
	experiments_list={'hist'}; %Just create one hist depending on what models use for the control
	model_list={'DOE_MALI','PIK_PISM2','UTAS_ElmerIce'}; 

	for iexp=1:length(experiments_list),
			expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			expename=experiments_list{iexp};
			time_vector=[];
			iareafl_vector=[];
			for isector=1:numsectors,
				eval(['iareafl_vector_sector' num2str(isector) '=[];'])
			end
			for iregion=1:numregions,
				eval(['iareafl_vector_region' num2str(iregion) '=[];'])
			end
			modelname=model_list{imodel};
			specifics;
			scalefac_model=scalefac(1:resolution:end,1:resolution:end);
			path=['' datapath '/' group '/' simul ''];

			if (strcmpi(modelname,'DOE_MALI') | strcmpi(modelname,'UTAS_ElmerIce')),
				expename='ctrl_proj_std'
			elseif strcmpi(modelname,'PIK_PISM2'),
				expename='ctrl_proj_open'
			end
			eval(['exp_name=' expename '_name;'])
			exp_directory=[path '/' exp_name '/'];
			if exist(exp_directory,'dir')==0,
				error(['directory ' exp_directory ' not found']);
			end

			field='sftflf';
			if is_sftflf==1,
				exp_file=[exp_directory '/' field '_AIS_' group '_' simul '_' expename '.nc'];
				expmask_file=[exp_directory '/sftgif_AIS_' group '_' simul '_' expename '.nc'];
				data = double(ncread(exp_file,field)); 
				mask = double(ncread(expmask_file,'sftgif')); 
				if(max(mask(:))>1 | min(mask(:))<0), disp(['mask should be between 0 and 1 in model ' modelname '']); end
				pos=find(mask>1);
				mask(pos)=1;
				time = double(ncread(exp_file,'time'))/yearday_model;  
				time=round(time,1);
				start_model_year=expe_model_year;
				time_vector(1:length(time))=start_model_year+time;
				pos_time=find(time_vector<=2015);
				if size(pos_time)>1, error('should not have several values before 2015'); end
				time_vector=time_vector(pos_time);
				for itime=1:length(time_vector),
					sftflf_i=data(:,:,itime);
					mask_i=mask(:,:,itime);
					pos=find(mask_i==0);
					sftflf_i(pos)=0;
					posnan=find(isnan(mask_i));
					sftflf_i(posnan)=0;
					mask_i(posnan)=0;
					posnan=find(isnan(sftflf_i)); 
					sftflf_i(posnan)=0; 
					mask_i(posnan)=0; 
					iareafl_total=sum(sftflf_i(:).*mask_i(:).*scalefac_model(:))*(resolution*1000)^2; %in m^2
					iareafl_vector(itime)=iareafl_total; 
					for isector=1:numsectors,
						eval(['sectors=sectors_' num2str(resolution) 'km;'])
						pos_sector=find(sectors==isector);
						iareafl_total_sector=sum(sftflf_i(pos_sector).*mask_i(pos_sector).*scalefac_model(pos_sector))*(resolution*1000)^2; %in m^2
						eval(['iareafl_vector_sector' num2str(isector) '(itime)=iareafl_total_sector;'])
					end
					for iregion=1:numregions,
						eval(['regions=regions_' num2str(resolution) 'km;'])
						pos_region=find(regions==iregion);
						iareafl_total_region=sum(sftflf_i(pos_region).*mask_i(pos_region).*scalefac_model(pos_region))*(resolution*1000)^2; %in m^2
						eval(['iareafl_vector_region' num2str(iregion) '(itime)=iareafl_total_region;'])
					end
				end %end of time

				if ~strcmpi(expename,'hist') & ~strcmpi(expename,'hist_open'),
					if length(time)<85, error(['run is too short: ' int2str(length(time)) ' years']); end
					if (size(data,3)<85 | size(data,3)>122), error('field has the wrong size'); end
				end
				if length(time)~=size(data,3), error(['lenght or time and data are not consistent in model ' modelname  ' for experiment ' exp_directory]); end
			end %end of lithk

			modelname
			size(start_model_year+time)
			size(time_vector)
			min(time_vector)
			max(time_vector)

			if (strcmpi(modelname,'DOE_MALI') | strcmpi(modelname,'UTAS_ElmerIce')),
				if exist(['ComputedScalarsPaper/' group '/' simul '/hist_std'],'dir')==0,
					eval(['mkdir ComputedScalarsPaper/' group '/' simul '/hist_std']);
				end
				expiareafl_file=['ComputedScalarsPaper/' group '/' simul '/hist_std/computed_iareafl_AIS_' group '_' simul '_hist_std.nc'];
			elseif (strcmpi(modelname,'PIK_PISM2')),
				if exist(['ComputedScalarsPaper/' group '/' simul '/hist_open'],'dir')==0,
					eval(['mkdir ComputedScalarsPaper/' group '/' simul '/hist_open']);
				end
				expiareafl_file=['ComputedScalarsPaper/' group '/' simul '/hist_open/computed_iareafl_AIS_' group '_' simul '_hist_open.nc'];
			end
			expiareafl_file
			status=WriteNetCDFComputedOutputs(expiareafl_file,'iareafl','floating ice extent','m^2',time_vector,iareafl_vector,...
				iareafl_vector_sector1,iareafl_vector_sector2,iareafl_vector_sector3,iareafl_vector_sector4,iareafl_vector_sector5,...
				iareafl_vector_sector6,iareafl_vector_sector7,iareafl_vector_sector8,iareafl_vector_sector9,iareafl_vector_sector10,...
				iareafl_vector_sector11,iareafl_vector_sector12,iareafl_vector_sector13,iareafl_vector_sector14,iareafl_vector_sector15,...
				iareafl_vector_sector16,iareafl_vector_sector17,iareafl_vector_sector18,...
				iareafl_vector_region1,iareafl_vector_region2,iareafl_vector_region3,ice_density,ocean_density);

		end %end of model
	end %end of experiment

end %}}}
if step==9, % {{{Compute_iareafl_exp_basins_minus_control

	for iexp=1:length(experiments_list),
			expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			expename=experiments_list{iexp};
			time_vector=[];
			iareafl_vector=[];
			for isector=1:numsectors,
				eval(['iareafl_vector_sector' num2str(isector) '=[];'])
			end
			for iregion=1:numregions,
				eval(['iareafl_vector_region' num2str(iregion) '=[];'])
			end
			modelname=model_list{imodel};
			specifics;
			scalefac_model=scalefac(1:resolution:end,1:resolution:end);
			path=['' datapath '/' group '/' simul ''];
			eval(['exp_name=' expename '_name;'])
			exp_directory=[path '/' exp_name '/'];

			eval(['isexp=is' expename ';'])
			if strcmpi(expename,'hist_std') | strcmpi(expename,'hist_open') | strcmpi(expename,'ctrl_proj_std') | strcmpi(expename,'ctrl_proj_open'),
				%do nothing, we cannot substract the control for the historical or control runs
			elseif isexp,
				%Prepare the ctrl values
				if strcmpi(expename,'exp01') | strcmpi(expename,'exp02') | strcmpi(expename,'exp03') | strcmpi(expename,'exp04') | strcmpi(expename,'exp11') | strcmpi(expename,'expA1') | strcmpi(expename,'expA2') | strcmpi(expename,'expA3') | strcmpi(expename,'expA4'),
					ctrl_file=['ComputedScalarsPaper/' group '/' simul '/ctrl_proj_open/computed_iareafl_AIS_' group '_' simul '_ctrl_proj_open.nc'];
				elseif strcmpi(expename,'exp05') | strcmpi(expename,'exp06') | strcmpi(expename,'exp07') | strcmpi(expename,'exp08') | strcmpi(expename,'exp09') | strcmpi(expename,'exp10') | strcmpi(expename,'exp12') | strcmpi(expename,'exp13') | strcmpi(expename,'expA5') | strcmpi(expename,'expA6') | strcmpi(expename,'expA7') | strcmpi(expename,'expA8'),
					ctrl_file=['ComputedScalarsPaper/' group '/' simul '/ctrl_proj_std/computed_iareafl_AIS_' group '_' simul '_ctrl_proj_std.nc'];
				else error('experiment not supported yet');
				end

				if exist(exp_directory,'dir')==0,
					error(['directory ' exp_directory ' not found']);
				end

				field='sftflf';
				if is_sftflf==1,
					exp_file=[exp_directory '/' field '_AIS_' group '_' simul '_' expename '.nc'];
					expmask_file=[exp_directory '/sftgif_AIS_' group '_' simul '_' expename '.nc'];
					data = double(ncread(exp_file,field)); 
					mask = double(ncread(expmask_file,'sftgif')); 
					if(max(mask(:))>1 | min(mask(:))<0), disp(['mask should be between 0 and 1 in model ' modelname '']); end
					pos=find(mask>1);
					mask(pos)=1;
					time = double(ncread(exp_file,'time'))/yearday_model;  
					if strcmpi(modelname,'VUW_PISM'),
						disp(['time should be in days in model ' modelname '']);
						time=time/(3600*24);
					end
					time=round(time,1);
					start_model_year=expe_model_year;
					time_vector(1:length(time))=start_model_year+time;
					if strcmpi(modelname,'ULB_FETISH32') | strcmpi(modelname,'ULB_FETISH16'),
						pos_time=find(time_vector>=2016 & time_vector<=2100);
					else
						pos_time=find(time_vector>=2016 & time_vector<=2101);
					end
					time_vector=time_vector(pos_time);
					for itime=1:length(time_vector),
						sftflf_i=data(:,:,itime);
						mask_i=mask(:,:,itime);
						pos=find(mask_i==0);
						sftflf_i(pos)=0;
						posnan=find(isnan(mask_i));
						sftflf_i(posnan)=0;
						mask_i(posnan)=0;
						posnan=find(isnan(sftflf_i)); 
						sftflf_i(posnan)=0; 
						mask_i(posnan)=0; 
						iareafl_total=sum(sftflf_i(:).*mask_i(:).*scalefac_model(:))*(resolution*1000)^2; %in m^2
						iareafl_vector(itime)=iareafl_total; 
						for isector=1:numsectors,
							eval(['sectors=sectors_' num2str(resolution) 'km;'])
							pos_sector=find(sectors==isector);
							iareafl_total_sector=sum(sftflf_i(pos_sector).*mask_i(pos_sector).*scalefac_model(pos_sector))*(resolution*1000)^2; %in m^2
							eval(['iareafl_total_sector_ctrlproj= double(ncread(ctrl_file,''iareafl_sector_' num2str(isector) '''));'])
							eval(['iareafl_vector_sector' num2str(isector) '(itime)=iareafl_total_sector-iareafl_total_sector_ctrlproj(itime);'])
						end
						for iregion=1:numregions,
							eval(['regions=regions_' num2str(resolution) 'km;'])
							pos_region=find(regions==iregion);
							iareafl_total_region=sum(sftflf_i(pos_region).*mask_i(pos_region).*scalefac_model(pos_region))*(resolution*1000)^2; %in m^2
							eval(['iareafl_total_region_ctrlproj= double(ncread(ctrl_file,''iareafl_region_' num2str(iregion) '''));'])
							eval(['iareafl_vector_region' num2str(iregion) '(itime)=iareafl_total_region-iareafl_total_region_ctrlproj(itime);'])
						end
					end %end of time
					iareafl_total_ctrlproj= double(ncread(ctrl_file,'iareafl'));
					if strcmpi(modelname,'DOE_MALI'),
						iareafl_vector=iareafl_vector-iareafl_total_ctrlproj(1:end-1)';
					else
						iareafl_vector=iareafl_vector-iareafl_total_ctrlproj';
					end

					if length(time)~=size(data,3), error(['lenght or time and data are not consistent in model ' modelname  ' for experiment ' exp_directory]); end
				end %end of lithk

				modelname
				size(start_model_year+time)
				size(time_vector)
				min(time_vector)
				max(time_vector)

				expiareafl_file=['ComputedScalarsPaper/' group '/' simul '/' expename '/computed_iareafl_minus_ctrl_proj_AIS_' group '_' simul '_' expename '.nc'];
				expiareafl_file
				status=WriteNetCDFComputedOutputs(expiareafl_file,'iareafl','floating ice extent','m^2',time_vector,iareafl_vector,...
					iareafl_vector_sector1,iareafl_vector_sector2,iareafl_vector_sector3,iareafl_vector_sector4,iareafl_vector_sector5,...
					iareafl_vector_sector6,iareafl_vector_sector7,iareafl_vector_sector8,iareafl_vector_sector9,iareafl_vector_sector10,...
					iareafl_vector_sector11,iareafl_vector_sector12,iareafl_vector_sector13,iareafl_vector_sector14,iareafl_vector_sector15,...
					iareafl_vector_sector16,iareafl_vector_sector17,iareafl_vector_sector18,...
					iareafl_vector_region1,iareafl_vector_region2,iareafl_vector_region3,ice_density,ocean_density);

			end %end of isexp
		end %end of model
	end %end of experiment

end %}}}
if step==10, % {{{Compute_ivol_exp_basins

	for iexp=1:length(experiments_list),
			expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			expename=experiments_list{iexp};
			time_vector=[];
			lim_vector=[];
			for isector=1:numsectors,
				eval(['lim_vector_sector' num2str(isector) '=[];'])
			end
			for iregion=1:numregions,
				eval(['lim_vector_region' num2str(iregion) '=[];'])
			end
			modelname=model_list{imodel};
			specifics;
			scalefac_model=scalefac(1:resolution:end,1:resolution:end);
			path=['' datapath '/' group '/' simul ''];
			eval(['exp_name=' expename '_name;'])
			exp_directory=[path '/' exp_name '/'];

			eval(['isexp=is' expename ';'])
			if isexp,
				if exist(exp_directory,'dir')==0,
					error(['directory ' exp_directory ' not found']);
				end

				field='lithk';
				if is_lithk==1,
					exp_file=[exp_directory '/' field '_AIS_' group '_' simul '_' expename '.nc'];
					expmask_file=[exp_directory '/sftgif_AIS_' group '_' simul '_' expename '.nc'];
					data = double(ncread(exp_file,field)); 
					mask = double(ncread(expmask_file,'sftgif')); 
					if(max(mask(:))>1 | min(mask(:))<0), disp(['mask should be between 0 and 1 in model ' modelname '']); end
					pos=find(mask>1);
					mask(pos)=1;
					time = double(ncread(exp_file,'time'))/yearday_model;  
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
						thickness_i=data(:,:,pos_time(itime));
						mask_i=mask(:,:,pos_time(itime));
						pos=find(mask_i==0);
						thickness_i(pos)=0;
						posnan=find(isnan(mask_i));
						thickness_i(posnan)=0;
						mask_i(posnan)=0;
						posnan=find(isnan(thickness_i));
						thickness_i(posnan)=0;
						mask_i(posnan)=0;
						vol=sum(thickness_i(:).*mask_i(:).*scalefac_model(:))*(resolution*1000)^2;
						lim_total=vol; %in m^3
						lim_vector(itime)=lim_total;
						for isector=1:numsectors,
							eval(['sectors=sectors_' num2str(resolution) 'km;'])
							pos_sector=find(sectors==isector);
							lim_total_sector=sum(thickness_i(pos_sector).*mask_i(pos_sector).*scalefac_model(pos_sector))*(resolution*1000)^2; %in m^3 
							eval(['lim_vector_sector' num2str(isector) '(itime)=lim_total_sector;'])
						end
						for iregion=1:numregions,
							eval(['regions=regions_' num2str(resolution) 'km;'])
							pos_region=find(regions==iregion);
							lim_total_region=sum(thickness_i(pos_region).*mask_i(pos_region).*scalefac_model(pos_region))*(resolution*1000)^2; %in m^3 
							eval(['lim_vector_region' num2str(iregion) '(itime)=lim_total_region;'])
						end
					end %end of time

					if ~strcmpi(expename,'hist_std') & ~strcmpi(expename,'hist_open'),
						if length(time)<85, error(['run is too short: ' int2str(length(time)) ' years']); end
						if (size(data,3)<85 | size(data,3)>122), error('field has the wrong size'); end
					end
					if length(time)~=size(data,3), error(['lenght or time and data are not consistent in model ' modelname  ' for experiment ' exp_directory]); end
				end %end of lithk

				modelname
				size(start_model_year+time)
				size(time_vector)
				min(time_vector)
				max(time_vector)

				if exist(['ComputedScalarsPaper/' group '/' simul '/' expename ''],'dir')==0,
					eval(['mkdir ComputedScalarsPaper/' group '/' simul '/' expename '']);
				end
				explim_file=['ComputedScalarsPaper/' group '/' simul '/' expename '/computed_ivol_AIS_' group '_' simul '_' expename '.nc'];
				explim_file
				status=WriteNetCDFComputedOutputs(explim_file,'ivol','ice volume','m^3',time_vector,lim_vector,...
					lim_vector_sector1,lim_vector_sector2,lim_vector_sector3,lim_vector_sector4,lim_vector_sector5,...
					lim_vector_sector6,lim_vector_sector7,lim_vector_sector8,lim_vector_sector9,lim_vector_sector10,...
					lim_vector_sector11,lim_vector_sector12,lim_vector_sector13,lim_vector_sector14,lim_vector_sector15,...
					lim_vector_sector16,lim_vector_sector17,lim_vector_sector18,...
					lim_vector_region1,lim_vector_region2,lim_vector_region3,...
					ice_density,ocean_density);

			end %end of isexp
		end %end of model
	end %end of experiment

end %}}}
if step==11, % {{{Compute_ivol_exp_basins_init_noinit

	%For the models that do not have a historical run, take the first step of the control to be the historical
	experiments_list={'hist'}; %Just create one hist depending on what models use for the control
	model_list={'DOE_MALI','PIK_PISM2','UTAS_ElmerIce'}; 

	for iexp=1:length(experiments_list),
			expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			expename=experiments_list{iexp};
			time_vector=[];
			lim_vector=[];
			for isector=1:numsectors,
				eval(['lim_vector_sector' num2str(isector) '=[];'])
			end
			for iregion=1:numregions,
				eval(['lim_vector_region' num2str(iregion) '=[];'])
			end
			modelname=model_list{imodel};
			specifics;
			scalefac_model=scalefac(1:resolution:end,1:resolution:end);
			path=['' datapath '/' group '/' simul ''];

			if (strcmpi(modelname,'DOE_MALI') | strcmpi(modelname,'UTAS_ElmerIce')),
				expename='ctrl_proj_std'
			elseif strcmpi(modelname,'PIK_PISM2'),
				expename='ctrl_proj_open'
			end
			eval(['exp_name=' expename '_name;'])
			exp_directory=[path '/' exp_name '/'];
			if exist(exp_directory,'dir')==0,
				error(['directory ' exp_directory ' not found']);
			end

			field='lithk';
			if is_lithk==1,
				exp_file=[exp_directory '/' field '_AIS_' group '_' simul '_' expename '.nc'];
				expmask_file=[exp_directory '/sftgif_AIS_' group '_' simul '_' expename '.nc'];
				data = double(ncread(exp_file,field)); 
				mask = double(ncread(expmask_file,'sftgif')); 
				if(max(mask(:))>1 | min(mask(:))<0), disp(['mask should be between 0 and 1 in model ' modelname '']); end
				pos=find(mask>1);
				mask(pos)=1;
				time = double(ncread(exp_file,'time'))/yearday_model;  
				time=round(time,1);
				start_model_year=expe_model_year;
				time_vector(1:length(time))=start_model_year+time;
				pos_time=find(time_vector<=2015);
				if size(pos_time)>1, error('should not have several values before 2015'); end
				time_vector=time_vector(pos_time);
				for itime=1:length(time_vector),
					thickness_i=data(:,:,pos_time(itime));
					mask_i=mask(:,:,pos_time(itime));
					pos=find(mask_i==0);
					thickness_i(pos)=0;
					posnan=find(isnan(mask_i));
					thickness_i(posnan)=0;
					mask_i(posnan)=0;
					posnan=find(isnan(thickness_i));
					thickness_i(posnan)=0;
					mask_i(posnan)=0;
					vol=sum(thickness_i(:).*mask_i(:).*scalefac_model(:))*(resolution*1000)^2;
					lim_total=vol; %in m^3 
					lim_vector(itime)=lim_total;
					for isector=1:numsectors,
						eval(['sectors=sectors_' num2str(resolution) 'km;'])
						pos_sector=find(sectors==isector);
						lim_total_sector=sum(thickness_i(pos_sector).*mask_i(pos_sector).*scalefac_model(pos_sector))*(resolution*1000)^2; %in m^3
						eval(['lim_vector_sector' num2str(isector) '(itime)=lim_total_sector;'])
					end
					for iregion=1:numregions,
						eval(['regions=regions_' num2str(resolution) 'km;'])
						pos_region=find(regions==iregion);
						lim_total_region=sum(thickness_i(pos_region).*mask_i(pos_region).*scalefac_model(pos_region))*(resolution*1000)^2; %in m^3
						eval(['lim_vector_region' num2str(iregion) '(itime)=lim_total_region;'])
					end
				end %end of time

				if length(time)~=size(data,3), error(['lenght or time and data are not consistent in model ' modelname  ' for experiment ' exp_directory]); end
			end %end of lithk

			modelname
			size(start_model_year+time)
			size(time_vector)
			min(time_vector)
			max(time_vector)

			if (strcmpi(modelname,'DOE_MALI') | strcmpi(modelname,'UTAS_ElmerIce')),
				if exist(['ComputedScalarsPaper/' group '/' simul '/hist_std'],'dir')==0,
					eval(['mkdir ComputedScalarsPaper/' group '/' simul '/hist_std']);
				end
				explim_file=['ComputedScalarsPaper/' group '/' simul '/hist_std/computed_ivol_AIS_' group '_' simul '_hist_std.nc'];
			elseif (strcmpi(modelname,'PIK_PISM2')),
				if exist(['ComputedScalarsPaper/' group '/' simul '/hist_open'],'dir')==0,
					eval(['mkdir ComputedScalarsPaper/' group '/' simul '/hist_open']);
				end
				explim_file=['ComputedScalarsPaper/' group '/' simul '/hist_open/computed_ivol_AIS_' group '_' simul '_hist_open.nc'];
			end
			explim_file
			status=WriteNetCDFComputedOutputs(explim_file,'ivol','ice volume','m^3',time_vector,lim_vector,...
				lim_vector_sector1,lim_vector_sector2,lim_vector_sector3,lim_vector_sector4,lim_vector_sector5,...
				lim_vector_sector6,lim_vector_sector7,lim_vector_sector8,lim_vector_sector9,lim_vector_sector10,...
				lim_vector_sector11,lim_vector_sector12,lim_vector_sector13,lim_vector_sector14,lim_vector_sector15,...
				lim_vector_sector16,lim_vector_sector17,lim_vector_sector18,...
				lim_vector_region1,lim_vector_region2,lim_vector_region3,...
				ice_density,ocean_density);

		end %end of model
	end %end of experiment

end %}}}
if step==12, % {{{Compute_ivol_exp_basins_minus_control

	for iexp=1:length(experiments_list),
			expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			expename=experiments_list{iexp};
			time_vector=[];
			lim_vector=[];
			for isector=1:numsectors,
				eval(['lim_vector_sector' num2str(isector) '=[];'])
			end
			for iregion=1:numregions,
				eval(['lim_vector_region' num2str(iregion) '=[];'])
			end
			modelname=model_list{imodel};
			specifics;
			scalefac_model=scalefac(1:resolution:end,1:resolution:end);
			path=['' datapath '/' group '/' simul ''];
			eval(['exp_name=' expename '_name;'])
			exp_directory=[path '/' exp_name '/'];

			eval(['isexp=is' expename ';'])
			if strcmpi(expename,'hist_std') | strcmpi(expename,'hist_open') | strcmpi(expename,'ctrl_proj_std') | strcmpi(expename,'ctrl_proj_open'),
				%do nothing, we cannot substract the control for the historical or control runs
			elseif isexp,
				%Prepare the ctrl values
				if strcmpi(expename,'exp01') | strcmpi(expename,'exp02') | strcmpi(expename,'exp03') | strcmpi(expename,'exp04') | strcmpi(expename,'exp11') | strcmpi(expename,'expA1') | strcmpi(expename,'expA2') | strcmpi(expename,'expA3') | strcmpi(expename,'expA4'),
					ctrl_file=['ComputedScalarsPaper/' group '/' simul '/ctrl_proj_open/computed_ivol_AIS_' group '_' simul '_ctrl_proj_open.nc'];
				elseif strcmpi(expename,'exp05') | strcmpi(expename,'exp06') | strcmpi(expename,'exp07') | strcmpi(expename,'exp08') | strcmpi(expename,'exp09') | strcmpi(expename,'exp10') | strcmpi(expename,'exp12') | strcmpi(expename,'exp13') | strcmpi(expename,'expA5') | strcmpi(expename,'expA6') | strcmpi(expename,'expA7') | strcmpi(expename,'expA8'),
					ctrl_file=['ComputedScalarsPaper/' group '/' simul '/ctrl_proj_std/computed_ivol_AIS_' group '_' simul '_ctrl_proj_std.nc'];
				else error('experiment not supported yet');
				end

				if exist(exp_directory,'dir')==0,
					error(['directory ' exp_directory ' not found']);
				end

				field='lithk';
				if is_lithk==1,
					exp_file=[exp_directory '/' field '_AIS_' group '_' simul '_' expename '.nc'];
					expmask_file=[exp_directory '/sftgif_AIS_' group '_' simul '_' expename '.nc'];
					data = double(ncread(exp_file,field)); 
					mask = double(ncread(expmask_file,'sftgif')); 
					if(max(mask(:))>1 | min(mask(:))<0), disp(['mask should be between 0 and 1 in model ' modelname '']); end
					pos=find(mask>1);
					mask(pos)=1;
					time = double(ncread(exp_file,'time'))/yearday_model;  
					if strcmpi(modelname,'VUW_PISM'),
						disp(['time should be in days in model ' modelname '']);
						time=time/(3600*24);
					end
					time=round(time,1);
					start_model_year=expe_model_year;
					time_vector(1:length(time))=start_model_year+time;
					if strcmpi(modelname,'ULB_FETISH32') | strcmpi(modelname,'ULB_FETISH16'),
						pos_time=find(time_vector>=2016 & time_vector<=2100);
					else
						pos_time=find(time_vector>=2016 & time_vector<=2101);
					end
					time_vector=time_vector(pos_time);
					for itime=1:length(time_vector),
						thickness_i=data(:,:,pos_time(itime));
						mask_i=mask(:,:,pos_time(itime));
						pos=find(mask_i==0);
						thickness_i(pos)=0;
						posnan=find(isnan(mask_i));
						thickness_i(posnan)=0;
						mask_i(posnan)=0;
						posnan=find(isnan(thickness_i));
						thickness_i(posnan)=0;
						mask_i(posnan)=0;
						vol=sum(thickness_i(:).*mask_i(:).*scalefac_model(:))*(resolution*1000)^2;
						lim_total=vol; %in m^3
						lim_vector(itime)=lim_total;
						for isector=1:numsectors,
							eval(['sectors=sectors_' num2str(resolution) 'km;'])
							pos_sector=find(sectors==isector);
							lim_total_sector=sum(thickness_i(pos_sector).*mask_i(pos_sector).*scalefac_model(pos_sector))*(resolution*1000)^2; %in m^3
							eval(['lim_total_sector_ctrlproj= double(ncread(ctrl_file,''ivol_sector_' num2str(isector) '''));'])
							eval(['lim_vector_sector' num2str(isector) '(itime)=lim_total_sector-lim_total_sector_ctrlproj(itime);'])
						end
						for iregion=1:numregions,
							eval(['regions=regions_' num2str(resolution) 'km;'])
							pos_region=find(regions==iregion);
							lim_total_region=sum(thickness_i(pos_region).*mask_i(pos_region).*scalefac_model(pos_region))*(resolution*1000)^2; %in m^3
							eval(['lim_total_region_ctrlproj= double(ncread(ctrl_file,''ivol_region_' num2str(iregion) '''));'])
							eval(['lim_vector_region' num2str(iregion) '(itime)=lim_total_region-lim_total_region_ctrlproj(itime);'])
						end
					end %end of time
					lim_total_ctrlproj= double(ncread(ctrl_file,'ivol'));
					if strcmpi(modelname,'DOE_MALI'),
						lim_vector=lim_vector-lim_total_ctrlproj(1:end-1)';
					else
						lim_vector=lim_vector-lim_total_ctrlproj';
					end

					if length(time)<85, error(['run is too short: ' int2str(length(time)) ' years']); end
					if (size(data,3)<85 | size(data,3)>122), error('field has the wrong size'); end
					if length(time)~=size(data,3), error(['lenght or time and data are not consistent in model ' modelname  ' for experiment ' exp_directory]); end
				end %end of lithk

				modelname
				size(start_model_year+time)
				size(time_vector)
				min(time_vector)
				max(time_vector)

				explim_file=['ComputedScalarsPaper/' group '/' simul '/' expename '/computed_lim_minus_ctrl_proj_AIS_' group '_' simul '_' expename '.nc'];
				explim_file
				status=WriteNetCDFComputedOutputs(explim_file,'ivol','ice volume','m^3',time_vector,lim_vector,...
					lim_vector_sector1,lim_vector_sector2,lim_vector_sector3,lim_vector_sector4,lim_vector_sector5,...
					lim_vector_sector6,lim_vector_sector7,lim_vector_sector8,lim_vector_sector9,lim_vector_sector10,...
					lim_vector_sector11,lim_vector_sector12,lim_vector_sector13,lim_vector_sector14,lim_vector_sector15,...
					lim_vector_sector16,lim_vector_sector17,lim_vector_sector18,...
					lim_vector_region1,lim_vector_region2,lim_vector_region3,...
					ice_density,ocean_density);

			end %end of isexp
		end %end of model
	end %end of experiment

end %}}}
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
			modelname=model_list{imodel};
			specifics;
			scalefac_model=scalefac(1:resolution:end,1:resolution:end);
			path=['' datapath '/' group '/' simul ''];
			eval(['exp_name=' expename '_name;'])
			exp_directory=[path '/' exp_name '/'];

			eval(['isexp=is' expename ';'])
			if isexp,
				if exist(exp_directory,'dir')==0,
					error(['directory ' exp_directory ' not found']);
				end

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
if step==14, % {{{Compute_ivaf_exp_basins_init_noinit

	%For the models that do not have a historical run, take the first step of the control to be the historical
	experiments_list={'hist'}; %Just create one hist depending on what models use for the control
	model_list={'DOE_MALI','PIK_PISM2','UTAS_ElmerIce'}; 

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
			modelname=model_list{imodel};
			specifics;
			scalefac_model=scalefac(1:resolution:end,1:resolution:end);
			path=['' datapath '/' group '/' simul ''];

			if (strcmpi(modelname,'DOE_MALI') | strcmpi(modelname,'UTAS_ElmerIce')),
				expename='ctrl_proj_std'
			elseif strcmpi(modelname,'PIK_PISM2'),
				expename='ctrl_proj_open'
			end
			eval(['exp_name=' expename '_name;'])
			exp_directory=[path '/' exp_name '/'];
			if exist(exp_directory,'dir')==0,
				error(['directory ' exp_directory ' not found']);
			end

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
				time=round(time,1);
				start_model_year=expe_model_year;
				time_vector(1:length(time))=start_model_year+time;
				pos_time=find(time_vector<=2015);
				if size(pos_time)>1, error('should not have several values before 2015'); end
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
						ivaf_total_region=sum((thickness_i(pos_region)+ocean_density/ice_density*min(bed_i(pos_region),0)).*groundmask_i(pos_region).*mask_i(pos_region).*scalefac_model(pos_region))*(resolution*1000)^2;  %in m^3
						eval(['ivaf_vector_region' num2str(iregion) '(itime)=ivaf_total_region;'])
					end
				end %end of time

				if length(time)~=size(thicknessdata,3), error(['lenght or time and data are not consistent in model ' modelname  ' for experiment ' exp_directory]); end
			end %end of lithk

			modelname
			size(start_model_year+time)
			size(time_vector)
			min(time_vector)
			max(time_vector)

			if (strcmpi(modelname,'DOE_MALI') | strcmpi(modelname,'UTAS_ElmerIce')),
				if exist(['ComputedScalarsPaper/' group '/' simul '/hist_std'],'dir')==0,
					eval(['mkdir ComputedScalarsPaper/' group '/' simul '/hist_std']);
				end
				expivaf_file=['ComputedScalarsPaper/' group '/' simul '/hist_std/computed_ivaf_AIS_' group '_' simul '_hist_std.nc'];
			elseif (strcmpi(modelname,'PIK_PISM2')),
				if exist(['ComputedScalarsPaper/' group '/' simul '/hist_open'],'dir')==0,
					eval(['mkdir ComputedScalarsPaper/' group '/' simul '/hist_open']);
				end
				expivaf_file=['ComputedScalarsPaper/' group '/' simul '/hist_open/computed_ivaf_AIS_' group '_' simul '_hist_open.nc'];
			end
			expivaf_file
			status=WriteNetCDFComputedOutputs(expivaf_file,'ivaf','ice volume above floatation','m^3',time_vector,ivaf_vector,...
				ivaf_vector_sector1,ivaf_vector_sector2,ivaf_vector_sector3,ivaf_vector_sector4,ivaf_vector_sector5,...
				ivaf_vector_sector6,ivaf_vector_sector7,ivaf_vector_sector8,ivaf_vector_sector9,ivaf_vector_sector10,...
				ivaf_vector_sector11,ivaf_vector_sector12,ivaf_vector_sector13,ivaf_vector_sector14,ivaf_vector_sector15,...
				ivaf_vector_sector16,ivaf_vector_sector17,ivaf_vector_sector18,...
				ivaf_vector_region1,ivaf_vector_region2,ivaf_vector_region3,...
				ice_density,ocean_density);

		end %end of model
	end %end of experiment

end %}}}
if step==15, % {{{Compute_ivaf_exp_basins_minus_control

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
			modelname=model_list{imodel};
			specifics;
			scalefac_model=scalefac(1:resolution:end,1:resolution:end);
			path=['' datapath '/' group '/' simul ''];
			eval(['exp_name=' expename '_name;'])
			exp_directory=[path '/' exp_name '/'];

			eval(['isexp=is' expename ';'])
			if strcmpi(expename,'hist_std') | strcmpi(expename,'hist_open') | strcmpi(expename,'ctrl_proj_std') | strcmpi(expename,'ctrl_proj_open'),
				%do nothing, we cannot substract the control for the historical or control runs
			elseif isexp,
				if strcmpi(expename,'exp01') | strcmpi(expename,'exp02') | strcmpi(expename,'exp03') | strcmpi(expename,'exp04') | strcmpi(expename,'exp11') | strcmpi(expename,'expA1') | strcmpi(expename,'expA2') | strcmpi(expename,'expA3') | strcmpi(expename,'expA4'),
					ctrl_file=['ComputedScalarsPaper/' group '/' simul '/ctrl_proj_open/computed_ivaf_AIS_' group '_' simul '_ctrl_proj_open.nc'];
				elseif strcmpi(expename,'exp05') | strcmpi(expename,'exp06') | strcmpi(expename,'exp07') | strcmpi(expename,'exp08') | strcmpi(expename,'exp09') | strcmpi(expename,'exp10') | strcmpi(expename,'exp12') | strcmpi(expename,'exp13') | strcmpi(expename,'expA5') | strcmpi(expename,'expA6') | strcmpi(expename,'expA7') | strcmpi(expename,'expA8'),
					ctrl_file=['ComputedScalarsPaper/' group '/' simul '/ctrl_proj_std/computed_ivaf_AIS_' group '_' simul '_ctrl_proj_std.nc'];
				else error('experiment not supported yet');
				end
				if exist(exp_directory,'dir')==0,
					error(['directory ' exp_directory ' not found']);
				end

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
					start_model_year=expe_model_year;
					time_vector(1:length(time))=start_model_year+time;
					if strcmpi(modelname,'ULB_FETISH32') | strcmpi(modelname,'ULB_FETISH16'),
						pos_time=find(time_vector>=2016 & time_vector<=2100);
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
							eval(['ivaf_total_sector_ctrlproj= double(ncread(ctrl_file,''ivaf_sector_' num2str(isector) '''));'])
							eval(['ivaf_vector_sector' num2str(isector) '(itime)=ivaf_total_sector-ivaf_total_sector_ctrlproj(itime);'])
						end
						for iregion=1:numregions,
							eval(['regions=regions_' num2str(resolution) 'km;'])
							pos_region=find(regions==iregion);
							ivaf_total_region=sum((thickness_i(pos_region)+ocean_density/ice_density*min(bed_i(pos_region),0)).*groundmask_i(pos_region).*mask_i(pos_region).*scalefac_model(pos_region))*(resolution*1000)^2; %in m^3
							eval(['limnws_total_region_ctrlproj= double(ncread(ctrl_file,''ivaf_region_' num2str(iregion) '''));'])
							eval(['ivaf_vector_region' num2str(iregion) '(itime)=ivaf_total_region-limnws_total_region_ctrlproj(itime);'])
						end
					end %end of time
					ivaf_total_ctrlproj= double(ncread(ctrl_file,'ivaf'));
					if strcmpi(modelname,'DOE_MALI'),
						ivaf_vector=ivaf_vector-ivaf_total_ctrlproj(1:end-1)';
					else
						ivaf_vector=ivaf_vector-ivaf_total_ctrlproj';
					end

					if length(time)<85, error(['run is too short: ' int2str(length(time)) ' years']); end
					if (size(thicknessdata,3)<85 | size(thicknessdata,3)>122), error('field has the wrong size'); end
					if length(time)~=size(thicknessdata,3), error(['lenght or time and data are not consistent in model ' modelname  ' for experiment ' exp_directory]); end
				end %end of lithk

				modelname
				size(start_model_year+time)
				size(time_vector)
				min(time_vector)
				max(time_vector)

				expivaf_file=['ComputedScalarsPaper/' group '/' simul '/' expename '/computed_ivaf_minus_ctrl_proj_AIS_' group '_' simul '_' expename '.nc'];
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
if step==16, % {{{Compute_smb_exp_basins
	
	for iexp=1:length(experiments_list),
			expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			expename=experiments_list{iexp};
			time_vector=[];
			tendacabf_vector=[];
			for isector=1:numsectors,
				eval(['tendacabf_vector_sector' num2str(isector) '=[];'])
			end
			for iregion=1:numregions,
				eval(['tendacabf_vector_region' num2str(iregion) '=[];'])
			end
			modelname=model_list{imodel};
			specifics;
			path=['' datapath '/' group '/' simul ''];
			scalefac_model=scalefac(1:resolution:end,1:resolution:end);
			eval(['exp_name=' expename '_name;'])
			exp_directory=[path '/' exp_name '/'];

			eval(['isexp=is' expename ';'])
			if isexp,
				if exist(exp_directory,'dir')==0,
					error(['directory ' exp_directory ' not found']);
				end

				field='acabf';
				if is_acabf==1,
					exp_file=[exp_directory '/' field '_AIS_' group '_' simul '_' expename '.nc'];
					expmask_file=[exp_directory '/sftgif_AIS_' group '_' simul '_' expename '.nc'];
					data = double(ncread(exp_file,field)); 
					mask = double(ncread(expmask_file,'sftgif')); 
					if(max(mask(:))>1 | min(mask(:))<0), disp(['mask should be between 0 and 1 in model ' modelname '']); end
					pos=find(mask>1);
					mask(pos)=1;
					time = double(ncread(exp_file,'time'))/yearday_model;  
					if strcmpi(modelname,'VUW_PISM'),
						disp(['time should be in days in model ' modelname '']);
						time=time/(3600*24);
					end
					time=round(time,1);
					if strcmpi(expename,'ctrl_std') | strcmpi(expename,'asmb') |strcmpi(expename,'abmb') | strcmpi(expename,'hist_std') | strcmpi(expename,'ctrl_open') | strcmpi(expename,'hist_open'),
						start_model_year=initial_model_year;
					else
						start_model_year=expe_model_year;
					end
					time_vector(1:length(time))=start_model_year+time;
					%Forcings should be reported at the middle of the year so fix that for models that did not do it
					if strcmpi(modelname,'JPL1_ISSM') | strcmpi(modelname,'NCAR_CISM') | strcmpi(modelname,'PIK_PISM1') | strcmpi(modelname,'PIK_PISM2') | strcmpi(modelname,'UCIJPL_ISSM') | ...
							strcmpi(modelname,'ULB_FETISH32') | strcmpi(modelname,'ULB_FETISH16') | strcmpi(modelname,'UTAS_ElmerIce') | strcmpi(modelname,'VUB_AISMPALEO') | strcmpi(modelname,'VUW_PISM'),
						time_vector=time_vector-0.5;
					end
					if strcmpi(expename,'hist_std') | strcmpi(expename,'hist_open'),
						pos_time=find(time_vector<=2015);
					else
						pos_time=find(time_vector>=2015 & time_vector<=2101);
					end
					time_vector=time_vector(pos_time);
					for itime=1:length(time_vector),
						acabf_i=data(:,:,itime);
						mask_i=mask(:,:,itime);
						pos=find(mask_i==0);
						acabf_i(pos)=0;
						posnan=find(isnan(mask_i));
						acabf_i(posnan)=0;
						mask_i(posnan)=0;
						posnan=find(isnan(acabf_i)); 
						acabf_i(posnan)=0; 
						mask_i(posnan)=0; 
						tendacabf_total=sum(acabf_i(:).*mask_i(:).*scalefac_model(:))*(resolution*1000)^2; %in kg/s
						tendacabf_vector(itime)=tendacabf_total; 
						for isector=1:numsectors,
							eval(['sectors=sectors_' num2str(resolution) 'km;'])
							pos_sector=find(sectors==isector);
							tendacabf_total_sector=sum(acabf_i(pos_sector).*mask_i(pos_sector).*scalefac_model(pos_sector))*(resolution*1000)^2; %in kg/s
							eval(['tendacabf_vector_sector' num2str(isector) '(itime)=tendacabf_total_sector;'])
						end
						for iregion=1:numregions,
							eval(['regions=regions_' num2str(resolution) 'km;'])
							pos_region=find(regions==iregion);
							tendacabf_total_region=sum(acabf_i(pos_region).*mask_i(pos_region).*scalefac_model(pos_region))*(resolution*1000)^2; %in kg/s
							eval(['tendacabf_vector_region' num2str(iregion) '(itime)=tendacabf_total_region;'])
						end
					end %end of time

					if ~strcmpi(expename,'hist_std') & ~strcmpi(expename,'hist_open'),
						if length(time)<85, error(['run is too short: ' int2str(length(time)) ' years']); end
						if (size(data,3)<85 | size(data,3)>122), error('field has the wrong size'); end
					end
					if length(time)~=size(data,3), error(['lenght or time and data are not consistent in model ' modelname  ' for experiment ' exp_directory]); end
				end %end of lithk

				modelname
				size(start_model_year+time)
				size(time_vector)
				min(time_vector)
				max(time_vector)

				if exist(['ComputedScalarsPaper/' group '/' simul '/' expename ''],'dir')==0,
					eval(['mkdir ComputedScalarsPaper/' group '/' simul '/' expename '']);
				end
				exptendacabf_file=['ComputedScalarsPaper/' group '/' simul '/' expename '/computed_smb_AIS_' group '_' simul '_' expename '.nc'];
				exptendacabf_file
				status=WriteNetCDFComputedOutputs(exptendacabf_file,'smb','spatially integrated surface mass balance','kg/s',time_vector,tendacabf_vector,...
					tendacabf_vector_sector1,tendacabf_vector_sector2,tendacabf_vector_sector3,tendacabf_vector_sector4,tendacabf_vector_sector5,...
					tendacabf_vector_sector6,tendacabf_vector_sector7,tendacabf_vector_sector8,tendacabf_vector_sector9,tendacabf_vector_sector10,...
					tendacabf_vector_sector11,tendacabf_vector_sector12,tendacabf_vector_sector13,tendacabf_vector_sector14,tendacabf_vector_sector15,...
					tendacabf_vector_sector16,tendacabf_vector_sector17,tendacabf_vector_sector18,...
					tendacabf_vector_region1,tendacabf_vector_region2,tendacabf_vector_region3,...
					ice_density,ocean_density);

			end %end of isexp
		end %end of model
	end %end of experiment

end %}}}
if step==17, % {{{Compute_smb_exp_basins_minus_control

	for iexp=1:length(experiments_list),
			expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			expename=experiments_list{iexp};
			time_vector=[];
			tendacabf_vector=[];
			for isector=1:numsectors,
				eval(['tendacabf_vector_sector' num2str(isector) '=[];'])
			end
			for iregion=1:numregions,
				eval(['tendacabf_vector_region' num2str(iregion) '=[];'])
			end
			modelname=model_list{imodel};
			specifics;
			scalefac_model=scalefac(1:resolution:end,1:resolution:end);
			path=['' datapath '/' group '/' simul ''];
			eval(['exp_name=' expename '_name;'])
			exp_directory=[path '/' exp_name '/'];

			eval(['isexp=is' expename ';'])
			if strcmpi(expename,'hist_std') | strcmpi(expename,'hist_open') | strcmpi(expename,'ctrl_proj_std') | strcmpi(expename,'ctrl_proj_open'),
				%do nothing, we cannot substract the control for the historical or control runs
			elseif isexp,
				%Prepare the ctrl values
				if strcmpi(expename,'exp01') | strcmpi(expename,'exp02') | strcmpi(expename,'exp03') | strcmpi(expename,'exp04') | strcmpi(expename,'exp11') | strcmpi(expename,'expA1') | strcmpi(expename,'expA2') | strcmpi(expename,'expA3') | strcmpi(expename,'expA4'),
					ctrl_file=['ComputedScalarsPaper/' group '/' simul '/ctrl_proj_open/computed_smb_AIS_' group '_' simul '_ctrl_proj_open.nc'];
				elseif strcmpi(expename,'exp05') | strcmpi(expename,'exp06') | strcmpi(expename,'exp07') | strcmpi(expename,'exp08') | strcmpi(expename,'exp09') | strcmpi(expename,'exp10') | strcmpi(expename,'exp12') | strcmpi(expename,'exp13') | strcmpi(expename,'expA5') | strcmpi(expename,'expA6') | strcmpi(expename,'expA7') | strcmpi(expename,'expA8'),
					ctrl_file=['ComputedScalarsPaper/' group '/' simul '/ctrl_proj_std/computed_smb_AIS_' group '_' simul '_ctrl_proj_std.nc'];
				else error('experiment not supported yet');
				end

				if exist(exp_directory,'dir')==0,
					error(['directory ' exp_directory ' not found']);
				end

				field='acabf';
				if is_acabf==1,
					exp_file=[exp_directory '/' field '_AIS_' group '_' simul '_' expename '.nc'];
					expmask_file=[exp_directory '/sftgif_AIS_' group '_' simul '_' expename '.nc'];
					data = double(ncread(exp_file,field)); 
					mask = double(ncread(expmask_file,'sftgif')); 
					if(max(mask(:))>1 | min(mask(:))<0), disp(['mask should be between 0 and 1 in model ' modelname '']); end
					pos=find(mask>1);
					mask(pos)=1;
					time = double(ncread(exp_file,'time'))/yearday_model;  
					if strcmpi(modelname,'VUW_PISM'),
						disp(['time should be in days in model ' modelname '']);
						time=time/(3600*24);
					end
					time=round(time,1);
					start_model_year=expe_model_year;
					time_vector(1:length(time))=start_model_year+time;
					%Forcings should be reported at the middle of the year so fix that for models that did not do it
					if strcmpi(modelname,'JPL1_ISSM') | strcmpi(modelname,'NCAR_CISM') | strcmpi(modelname,'PIK_PISM1') | strcmpi(modelname,'PIK_PISM2') | strcmpi(modelname,'UCIJPL_ISSM') | ...
							strcmpi(modelname,'ULB_FETISH32') | strcmpi(modelname,'ULB_FETISH16') | strcmpi(modelname,'UTAS_ElmerIce') | strcmpi(modelname,'VUB_AISMPALEO') | strcmpi(modelname,'VUW_PISM'),
						time_vector=time_vector-0.5;
					end
					if strcmpi(modelname,'ULB_FETISH32') | strcmpi(modelname,'ULB_FETISH16'),
						pos_time=find(time_vector>=2015 & time_vector<=2100);
					else
						pos_time=find(time_vector>=2015 & time_vector<=2101);
					end
					time_vector=time_vector(pos_time);
					for itime=1:length(time_vector),
						acabf_i=data(:,:,itime);
						mask_i=mask(:,:,itime);
						pos=find(mask_i==0);
						acabf_i(pos)=0;
						posnan=find(isnan(mask_i));
						acabf_i(posnan)=0;
						mask_i(posnan)=0;
						posnan=find(isnan(acabf_i)); 
						acabf_i(posnan)=0; 
						mask_i(posnan)=0; 
						tendacabf_total=sum(acabf_i(:).*mask_i(:).*scalefac_model(:))*(resolution*1000)^2; %in kg/s
						tendacabf_vector(itime)=tendacabf_total; 
						for isector=1:numsectors,
							eval(['sectors=sectors_' num2str(resolution) 'km;'])
							pos_sector=find(sectors==isector);
							tendacabf_total_sector=sum(acabf_i(pos_sector).*mask_i(pos_sector).*scalefac_model(pos_sector))*(resolution*1000)^2; %in kg/s
							eval(['tendacabf_total_sector_ctrlproj= double(ncread(ctrl_file,''smb_sector_' num2str(isector) '''));'])
							eval(['tendacabf_vector_sector' num2str(isector) '(itime)=tendacabf_total_sector-tendacabf_total_sector_ctrlproj(itime);'])
						end
						for iregion=1:numregions,
							eval(['regions=regions_' num2str(resolution) 'km;'])
							pos_region=find(regions==iregion);
							tendacabf_total_region=sum(acabf_i(pos_region).*mask_i(pos_region).*scalefac_model(pos_region))*(resolution*1000)^2; %in kg/s
							eval(['tendacabf_total_region_ctrlproj= double(ncread(ctrl_file,''smb_region_' num2str(iregion) '''));'])
							eval(['tendacabf_vector_region' num2str(iregion) '(itime)=tendacabf_total_region-tendacabf_total_region_ctrlproj(itime);'])
						end
					end %end of time
					tendacabf_total_ctrlproj= double(ncread(ctrl_file,'smb'));
					tendacabf_vector=tendacabf_vector-tendacabf_total_ctrlproj';

					if length(time)<85, error(['run is too short: ' int2str(length(time)) ' years']); end
					if (size(data,3)<85 | size(data,3)>122), error('field has the wrong size'); end
					if length(time)~=size(data,3), error(['lenght or time and data are not consistent in model ' modelname  ' for experiment ' exp_directory]); end
				end %end of lithk

				modelname
				size(start_model_year+time)
				size(time_vector)
				min(time_vector)
				max(time_vector)

				exptendacabf_file=['ComputedScalarsPaper/' group '/' simul '/' expename '/computed_smb_minus_ctrl_proj_AIS_' group '_' simul '_' expename '.nc'];
				exptendacabf_file
				status=WriteNetCDFComputedOutputs(exptendacabf_file,'smb','spatially integrated surface mass balance','kg/s',time_vector,tendacabf_vector,...
					tendacabf_vector_sector1,tendacabf_vector_sector2,tendacabf_vector_sector3,tendacabf_vector_sector4,tendacabf_vector_sector5,...
					tendacabf_vector_sector6,tendacabf_vector_sector7,tendacabf_vector_sector8,tendacabf_vector_sector9,tendacabf_vector_sector10,...
					tendacabf_vector_sector11,tendacabf_vector_sector12,tendacabf_vector_sector13,tendacabf_vector_sector14,tendacabf_vector_sector15,...
					tendacabf_vector_sector16,tendacabf_vector_sector17,tendacabf_vector_sector18,...
					tendacabf_vector_region1,tendacabf_vector_region2,tendacabf_vector_region3,...
					ice_density,ocean_density);

			end %end of isexp
		end %end of model
	end %end of experiment

end %}}}
if step==18, % {{{Compute_smbgr_exp_basins

	for iexp=1:length(experiments_list),
			expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			expename=experiments_list{iexp};
			time_vector=[];
			smbgr_vector=[];
			for isector=1:numsectors,
				eval(['smbgr_vector_sector' num2str(isector) '=[];'])
			end
			for iregion=1:numregions,
				eval(['smbgr_vector_region' num2str(iregion) '=[];'])
			end
			modelname=model_list{imodel};
			specifics;
			scalefac_model=scalefac(1:resolution:end,1:resolution:end);
			path=['' datapath '/' group '/' simul ''];
			eval(['exp_name=' expename '_name;'])
			exp_directory=[path '/' exp_name '/'];

			eval(['isexp=is' expename ';'])
			if isexp,
				if exist(exp_directory,'dir')==0,
					error(['directory ' exp_directory ' not found']);
				end

				field='acabf';
				if is_acabf==1,
					exp_file=[exp_directory '/' field '_AIS_' group '_' simul '_' expename '.nc'];
					expmask_file=[exp_directory '/sftgif_AIS_' group '_' simul '_' expename '.nc'];
					expgrmask_file=[exp_directory '/sftgrf_AIS_' group '_' simul '_' expename '.nc'];
					data = double(ncread(exp_file,field)); 
					mask = double(ncread(expmask_file,'sftgif')); 
					grmask = double(ncread(expgrmask_file,'sftgrf')); 
					if(max(mask(:))>1 | min(mask(:))<0), disp(['mask should be between 0 and 1 in model ' modelname '']); end
					pos=find(mask>1);
					mask(pos)=1;
					pos=find(grmask>1);
					grmask(pos)=1;
					time = double(ncread(exp_file,'time'))/yearday_model;  
					if strcmpi(modelname,'VUW_PISM'),
						disp(['time should be in days in model ' modelname '']);
						time=time/(3600*24);
					end
					time=round(time,1);
					if strcmpi(expename,'ctrl_std') | strcmpi(expename,'asmb') |strcmpi(expename,'abmb') | strcmpi(expename,'hist_std') | strcmpi(expename,'ctrl_open') | strcmpi(expename,'hist_open'),
						start_model_year=initial_model_year;
					else
						start_model_year=expe_model_year;
					end
					time_vector(1:length(time))=start_model_year+time;
					%Forcings should be reported at the middle of the year so fix that for models that did not do it
					if strcmpi(modelname,'JPL1_ISSM') | strcmpi(modelname,'NCAR_CISM') | strcmpi(modelname,'PIK_PISM1') | strcmpi(modelname,'PIK_PISM2') | strcmpi(modelname,'UCIJPL_ISSM') | ...
							strcmpi(modelname,'ULB_FETISH32') | strcmpi(modelname,'ULB_FETISH16') | strcmpi(modelname,'UTAS_ElmerIce') | strcmpi(modelname,'VUB_AISMPALEO') | strcmpi(modelname,'VUW_PISM'),
						time_vector=time_vector-0.5;
					end
					if strcmpi(expename,'hist_std') | strcmpi(expename,'hist_open'),
						pos_time=find(time_vector<=2015);
					else
						pos_time=find(time_vector>=2015 & time_vector<=2101);
					end
					time_vector=time_vector(pos_time);
					for itime=1:length(time_vector),
						acabf_i=data(:,:,itime);
						mask_i=mask(:,:,itime);
						grmask_i=grmask(:,:,itime);
						pos=find(mask_i==0);
						acabf_i(pos)=0;
						grmask_i(pos)=0;
						posnan=find(isnan(mask_i));
						acabf_i(posnan)=0;
						mask_i(posnan)=0;
						grmask_i(posnan)=0;
						posnan=find(isnan(grmask_i));
						acabf_i(posnan)=0;
						mask_i(posnan)=0;
						grmask_i(posnan)=0;
						posnan=find(isnan(acabf_i));
						acabf_i(posnan)=0;
						mask_i(posnan)=0;
						grmask_i(posnan)=0;
						smbgr_total=sum(acabf_i(:).*grmask_i(:).*mask_i(:).*scalefac_model(:))*(resolution*1000)^2; %in kg/s
						smbgr_vector(itime)=smbgr_total;
						for isector=1:numsectors,
							eval(['sectors=sectors_' num2str(resolution) 'km;'])
							pos_sector=find(sectors==isector);
							smbgr_total_sector=sum(acabf_i(pos_sector).*grmask_i(pos_sector).*mask_i(pos_sector).*scalefac_model(pos_sector))*(resolution*1000)^2; %in kg/s 
							eval(['smbgr_vector_sector' num2str(isector) '(itime)=smbgr_total_sector;'])
						end
						for iregion=1:numregions,
							eval(['regions=regions_' num2str(resolution) 'km;'])
							pos_region=find(regions==iregion);
							smbgr_total_region=sum(acabf_i(pos_region).*grmask_i(pos_region).*mask_i(pos_region).*scalefac_model(pos_region))*(resolution*1000)^2; %in kg/s 
							eval(['smbgr_vector_region' num2str(iregion) '(itime)=smbgr_total_region;'])
						end
					end %end of time

					if ~strcmpi(expename,'hist_std') & ~strcmpi(expename,'hist_open'),
						if length(time)<85, error(['run is too short: ' int2str(length(time)) ' years']); end
						if (size(data,3)<85 | size(data,3)>122), error('field has the wrong size'); end
					end
					if length(time)~=size(data,3), error(['lenght or time and data are not consistent in model ' modelname  ' for experiment ' exp_directory]); end
				end %end of lithk

				modelname
				size(start_model_year+time)
				size(time_vector)
				min(time_vector)
				max(time_vector)

				if exist(['ComputedScalarsPaper/' group '/' simul '/' expename ''],'dir')==0,
					eval(['mkdir ComputedScalarsPaper/' group '/' simul '/' expename '']);
				end
				expsmbgr_file=['ComputedScalarsPaper/' group '/' simul '/' expename '/computed_smbgr_AIS_' group '_' simul '_' expename '.nc'];
				expsmbgr_file
				status=WriteNetCDFComputedOutputs(expsmbgr_file,'smbgr','spatially integrated grounded ice surface mass balance','kg/s',time_vector,smbgr_vector,...
					smbgr_vector_sector1,smbgr_vector_sector2,smbgr_vector_sector3,smbgr_vector_sector4,smbgr_vector_sector5,...
					smbgr_vector_sector6,smbgr_vector_sector7,smbgr_vector_sector8,smbgr_vector_sector9,smbgr_vector_sector10,...
					smbgr_vector_sector11,smbgr_vector_sector12,smbgr_vector_sector13,smbgr_vector_sector14,smbgr_vector_sector15,...
					smbgr_vector_sector16,smbgr_vector_sector17,smbgr_vector_sector18,...
					smbgr_vector_region1,smbgr_vector_region2,smbgr_vector_region3,...
					ice_density,ocean_density);
			end %end of isexp
		end %end of model
	end %end of experiment

end %}}}
if step==19, % {{{Compute_smbgr_exp_basins_minus_control

	for iexp=1:length(experiments_list),
			expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			expename=experiments_list{iexp};
			time_vector=[];
			smbgr_vector=[];
			for isector=1:numsectors,
				eval(['smbgr_vector_sector' num2str(isector) '=[];'])
			end
			for iregion=1:numregions,
				eval(['smbgr_vector_region' num2str(iregion) '=[];'])
			end
			modelname=model_list{imodel};
			specifics;
			scalefac_model=scalefac(1:resolution:end,1:resolution:end);
			path=['' datapath '/' group '/' simul ''];
			eval(['exp_name=' expename '_name;'])
			exp_directory=[path '/' exp_name '/'];

			eval(['isexp=is' expename ';'])
			if strcmpi(expename,'hist_std') | strcmpi(expename,'hist_open') | strcmpi(expename,'ctrl_proj_std') | strcmpi(expename,'ctrl_proj_open'),
				%do nothing, we cannot substract the control for the historical or control runs
			elseif isexp,
				%Prepare the ctrl values
				if strcmpi(expename,'exp01') | strcmpi(expename,'exp02') | strcmpi(expename,'exp03') | strcmpi(expename,'exp04') | strcmpi(expename,'exp11') | strcmpi(expename,'expA1') | strcmpi(expename,'expA2') | strcmpi(expename,'expA3') | strcmpi(expename,'expA4'),
					ctrl_file=['ComputedScalarsPaper/' group '/' simul '/ctrl_proj_open/computed_smbgr_AIS_' group '_' simul '_ctrl_proj_open.nc'];
				elseif strcmpi(expename,'exp05') | strcmpi(expename,'exp06') | strcmpi(expename,'exp07') | strcmpi(expename,'exp08') | strcmpi(expename,'exp09') | strcmpi(expename,'exp10') | strcmpi(expename,'exp12') | strcmpi(expename,'exp13') | strcmpi(expename,'expA5') | strcmpi(expename,'expA6') | strcmpi(expename,'expA7') | strcmpi(expename,'expA8'),
					ctrl_file=['ComputedScalarsPaper/' group '/' simul '/ctrl_proj_std/computed_smbgr_AIS_' group '_' simul '_ctrl_proj_std.nc'];
				else error('experiment not supported yet');
				end

				if exist(exp_directory,'dir')==0,
					error(['directory ' exp_directory ' not found']);
				end

				field='acabf';
				if is_acabf==1,
					exp_file=[exp_directory '/' field '_AIS_' group '_' simul '_' expename '.nc'];
					expmask_file=[exp_directory '/sftgif_AIS_' group '_' simul '_' expename '.nc'];
					expgrmask_file=[exp_directory '/sftgrf_AIS_' group '_' simul '_' expename '.nc'];
					data = double(ncread(exp_file,field)); 
					mask = double(ncread(expmask_file,'sftgif')); 
					grmask = double(ncread(expgrmask_file,'sftgrf')); 
					if(max(mask(:))>1 | min(mask(:))<0), disp(['mask should be between 0 and 1 in model ' modelname '']); end
					pos=find(mask>1);
					mask(pos)=1;
					pos=find(grmask>1);
					grmask(pos)=1;
					time = double(ncread(exp_file,'time'))/yearday_model;  
					if strcmpi(modelname,'VUW_PISM'),
						disp(['time should be in days in model ' modelname '']);
						time=time/(3600*24);
					end
					time=round(time,1);
					start_model_year=expe_model_year;
					time_vector(1:length(time))=start_model_year+time;
					%Forcings should be reported at the middle of the year so fix that for models that did not do it
					if strcmpi(modelname,'JPL1_ISSM') | strcmpi(modelname,'NCAR_CISM') | strcmpi(modelname,'PIK_PISM1') | strcmpi(modelname,'PIK_PISM2') | strcmpi(modelname,'UCIJPL_ISSM') | ...
							strcmpi(modelname,'ULB_FETISH32') | strcmpi(modelname,'ULB_FETISH16') | strcmpi(modelname,'UTAS_ElmerIce') | strcmpi(modelname,'VUB_AISMPALEO') | strcmpi(modelname,'VUW_PISM'),
						time_vector=time_vector-0.5;
					end
					if strcmpi(modelname,'ULB_FETISH32') | strcmpi(modelname,'ULB_FETISH16'),
						pos_time=find(time_vector>=2015 & time_vector<=2100);
					else
						pos_time=find(time_vector>=2015 & time_vector<=2101);
					end
					time_vector=time_vector(pos_time);
					for itime=1:length(time_vector),
						acabf_i=data(:,:,itime);
						mask_i=mask(:,:,itime);
						grmask_i=grmask(:,:,itime);
						pos=find(mask_i==0);
						acabf_i(pos)=0;
						grmask_i(pos)=0;
						posnan=find(isnan(mask_i));
						acabf_i(posnan)=0;
						mask_i(posnan)=0;
						grmask_i(posnan)=0;
						posnan=find(isnan(grmask_i));
						acabf_i(posnan)=0;
						mask_i(posnan)=0;
						grmask_i(posnan)=0;
						posnan=find(isnan(acabf_i));
						acabf_i(posnan)=0;
						mask_i(posnan)=0;
						grmask_i(posnan)=0;
						smbgr_total=sum(acabf_i(:).*grmask_i(:).*mask_i(:).*scalefac_model(:))*(resolution*1000)^2; %in kg/s
						smbgr_vector(itime)=smbgr_total;
						for isector=1:numsectors,
							eval(['sectors=sectors_' num2str(resolution) 'km;'])
							pos_sector=find(sectors==isector);
							smbgr_total_sector=sum(acabf_i(pos_sector).*grmask_i(pos_sector).*mask_i(pos_sector).*scalefac_model(pos_sector))*(resolution*1000)^2; %in kg/s
							eval(['smbgr_total_sector_ctrlproj= double(ncread(ctrl_file,''smbgr_sector_' num2str(isector) '''));'])
							eval(['smbgr_vector_sector' num2str(isector) '(itime)=smbgr_total_sector-smbgr_total_sector_ctrlproj(itime);'])
						end
						for iregion=1:numregions,
							eval(['regions=regions_' num2str(resolution) 'km;'])
							pos_region=find(regions==iregion);
							smbgr_total_region=sum(acabf_i(pos_region).*grmask_i(pos_region).*mask_i(pos_region).*scalefac_model(pos_region))*(resolution*1000)^2; %in kg/s
							eval(['smbgr_total_region_ctrlproj= double(ncread(ctrl_file,''smbgr_region_' num2str(iregion) '''));'])
							eval(['smbgr_vector_region' num2str(iregion) '(itime)=smbgr_total_region-smbgr_total_region_ctrlproj(itime);'])
						end
					end %end of time
					smbgr_total_ctrlproj= double(ncread(ctrl_file,'smbgr'));
					smbgr_vector=smbgr_vector-smbgr_total_ctrlproj';

					if length(time)<85, error(['run is too short: ' int2str(length(time)) ' years']); end
					if (size(data,3)<85 | size(data,3)>122), error('field has the wrong size'); end
					if length(time)~=size(data,3), error(['lenght or time and data are not consistent in model ' modelname  ' for experiment ' exp_directory]); end
				end %end of lithk

				modelname
				size(start_model_year+time)
				size(time_vector)
				min(time_vector)
				max(time_vector)

				expsmbgr_file=['ComputedScalarsPaper/' group '/' simul '/' expename '/computed_smbgr_minus_ctrl_proj_AIS_' group '_' simul '_' expename '.nc'];
				expsmbgr_file
				status=WriteNetCDFComputedOutputs(expsmbgr_file,'smbgr','spatially integrated grounded surface mass balance','kg/s',time_vector,smbgr_vector,...
					smbgr_vector_sector1,smbgr_vector_sector2,smbgr_vector_sector3,smbgr_vector_sector4,smbgr_vector_sector5,...
					smbgr_vector_sector6,smbgr_vector_sector7,smbgr_vector_sector8,smbgr_vector_sector9,smbgr_vector_sector10,...
					smbgr_vector_sector11,smbgr_vector_sector12,smbgr_vector_sector13,smbgr_vector_sector14,smbgr_vector_sector15,...
					smbgr_vector_sector16,smbgr_vector_sector17,smbgr_vector_sector18,...
					smbgr_vector_region1,smbgr_vector_region2,smbgr_vector_region3,...
					ice_density,ocean_density);

			end %end of isexp
		end %end of model
	end %end of experiment

end %}}}
if step==20, % {{{Compute_bmbfl_exp_basins

	for iexp=1:length(experiments_list),
			expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			expename=experiments_list{iexp};
			time_vector=[];
			bmbfl_vector=[];
			for isector=1:numsectors,
				eval(['bmbfl_vector_sector' num2str(isector) '=[];'])
			end
			for iregion=1:numregions,
				eval(['bmbfl_vector_region' num2str(iregion) '=[];'])
			end
			modelname=model_list{imodel};
			specifics;
			scalefac_model=scalefac(1:resolution:end,1:resolution:end);
			path=['' datapath '/' group '/' simul ''];
			eval(['exp_name=' expename '_name;'])
			exp_directory=[path '/' exp_name '/'];

			eval(['isexp=is' expename ';'])
			if isexp,
				if exist(exp_directory,'dir')==0,
					error(['directory ' exp_directory ' not found']);
				end

				field='libmassbffl';
				if is_libmassbffl==1,
					exp_file=[exp_directory '/' field '_AIS_' group '_' simul '_' expename '.nc'];
					expmask_file=[exp_directory '/sftgif_AIS_' group '_' simul '_' expename '.nc'];
					expshelfmask_file=[exp_directory '/sftflf_AIS_' group '_' simul '_' expename '.nc'];
					data = double(ncread(exp_file,field)); 
					mask = double(ncread(expmask_file,'sftgif')); 
					shelfmask = double(ncread(expshelfmask_file,'sftflf')); 
					if(max(mask(:))>1 | min(mask(:))<0), disp(['mask should be between 0 and 1 in model ' modelname '']); end
					pos=find(mask>1);
					mask(pos)=1;
					pos=find(shelfmask>1);
					shelfmask(pos)=1;
					time = double(ncread(exp_file,'time'))/yearday_model;  
					if strcmpi(modelname,'VUW_PISM'),
						disp(['time should be in days in model ' modelname '']);
						time=time/(3600*24);
					end
					time=round(time,1);
					if strcmpi(expename,'ctrl_std') | strcmpi(expename,'asmb') |strcmpi(expename,'abmb') | strcmpi(expename,'hist_std') | strcmpi(expename,'ctrl_open') | strcmpi(expename,'hist_open'),
						start_model_year=initial_model_year;
					else
						start_model_year=expe_model_year;
					end
					time_vector(1:length(time))=start_model_year+time;
					%Forcings should be reported at the middle of the year so fix that for models that did not do it
					if strcmpi(modelname,'JPL1_ISSM') | strcmpi(modelname,'NCAR_CISM') | strcmpi(modelname,'PIK_PISM1') | strcmpi(modelname,'PIK_PISM2') | strcmpi(modelname,'UCIJPL_ISSM') | ...
							strcmpi(modelname,'ULB_FETISH32') | strcmpi(modelname,'ULB_FETISH16') | strcmpi(modelname,'UTAS_ElmerIce') | strcmpi(modelname,'VUB_AISMPALEO') | strcmpi(modelname,'VUW_PISM'),
						time_vector=time_vector-0.5;
					end
					if strcmpi(expename,'hist_std') | strcmpi(expename,'hist_open'),
						pos_time=find(time_vector<=2015);
					else
						pos_time=find(time_vector>=2015 & time_vector<=2101);
					end
					time_vector=time_vector(pos_time);
					for itime=1:length(time_vector),
						libmassbffl_i=data(:,:,itime);
						mask_i=mask(:,:,itime);
						shelfmask_i=shelfmask(:,:,itime);
						pos=find(mask_i==0);
						libmassbffl_i(pos)=0;
						shelfmask_i(pos)=0;
						posnan=find(isnan(mask_i));
						libmassbffl_i(posnan)=0;
						mask_i(posnan)=0;
						shelfmask_i(posnan)=0;
						posnan=find(isnan(libmassbffl_i));
						libmassbffl_i(posnan)=0;
						mask_i(posnan)=0;
						shelfmask_i(posnan)=0;
						posnan=find(isnan(shelfmask_i));
						libmassbffl_i(posnan)=0;
						mask_i(posnan)=0;
						shelfmask_i(posnan)=0;
						bmbfl_total=sum(libmassbffl_i(:).*shelfmask_i(:).*mask_i(:).*scalefac_model(:))*(resolution*1000)^2; %in kg/s
						if strcmpi(modelname,'VUW_PISM'),
							bmbfl_total=-bmbfl_total; %wrong sign in libmassbffl
						end
						bmbfl_vector(itime)=bmbfl_total;
						for isector=1:numsectors,
							eval(['sectors=sectors_' num2str(resolution) 'km;'])
							pos_sector=find(sectors==isector);
							bmbfl_total_sector=sum(libmassbffl_i(pos_sector).*shelfmask_i(pos_sector).*mask_i(pos_sector).*scalefac_model(pos_sector))*(resolution*1000)^2; %in kg/s 
							eval(['bmbfl_vector_sector' num2str(isector) '(itime)=bmbfl_total_sector;'])
						end
						for iregion=1:numregions,
							eval(['regions=regions_' num2str(resolution) 'km;'])
							pos_region=find(regions==iregion);
							bmbfl_total_region=sum(libmassbffl_i(pos_region).*shelfmask_i(pos_region).*mask_i(pos_region).*scalefac_model(pos_region))*(resolution*1000)^2; %in kg/s
							eval(['bmbfl_vector_region' num2str(iregion) '(itime)=bmbfl_total_region;'])
						end
					end %end of time

					if ~strcmpi(expename,'hist_std') & ~strcmpi(expename,'hist_open'),
						if length(time)<85, error(['run is too short: ' int2str(length(time)) ' years']); end
						if (size(data,3)<85 | size(data,3)>122), error('field has the wrong size'); end
					end
					if length(time)~=size(data,3), error(['lenght or time and data are not consistent in model ' modelname  ' for experiment ' exp_directory]); end
				end %end of lithk

				modelname
				size(start_model_year+time)
				size(time_vector)
				min(time_vector)
				max(time_vector)

				if exist(['ComputedScalarsPaper/' group '/' simul '/' expename ''],'dir')==0,
					eval(['mkdir ComputedScalarsPaper/' group '/' simul '/' expename '']);
				end
				expbmbfl_file=['ComputedScalarsPaper/' group '/' simul '/' expename '/computed_bmbfl_AIS_' group '_' simul '_' expename '.nc'];
				expbmbfl_file
				status=WriteNetCDFComputedOutputs(expbmbfl_file,'bmbfl','spatially integrated floating ice basal melt','kg/s',time_vector,bmbfl_vector,...
					bmbfl_vector_sector1,bmbfl_vector_sector2,bmbfl_vector_sector3,bmbfl_vector_sector4,bmbfl_vector_sector5,...
					bmbfl_vector_sector6,bmbfl_vector_sector7,bmbfl_vector_sector8,bmbfl_vector_sector9,bmbfl_vector_sector10,...
					bmbfl_vector_sector11,bmbfl_vector_sector12,bmbfl_vector_sector13,bmbfl_vector_sector14,bmbfl_vector_sector15,...
					bmbfl_vector_sector16,bmbfl_vector_sector17,bmbfl_vector_sector18,...
					bmbfl_vector_region1,bmbfl_vector_region2,bmbfl_vector_region3,...
					ice_density,ocean_density);

			end %end of isexp
		end %end of model
	end %end of experiment

end %}}}
if step==21, % {{{Compute_bmbfl_exp_basins_minus_control

	for iexp=1:length(experiments_list),
			expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			expename=experiments_list{iexp};
			time_vector=[];
			bmbfl_vector=[];
			for isector=1:numsectors,
				eval(['bmbfl_vector_sector' num2str(isector) '=[];'])
			end
			for iregion=1:numregions,
				eval(['bmbfl_vector_region' num2str(iregion) '=[];'])
			end
			modelname=model_list{imodel};
			specifics;
			scalefac_model=scalefac(1:resolution:end,1:resolution:end);
			path=['' datapath '/' group '/' simul ''];
			eval(['exp_name=' expename '_name;'])
			exp_directory=[path '/' exp_name '/'];

			eval(['isexp=is' expename ';'])
			if strcmpi(expename,'hist_std') | strcmpi(expename,'hist_open') | strcmpi(expename,'ctrl_proj_std') | strcmpi(expename,'ctrl_proj_open'),
				%do nothing, we cannot substract the control for the historical or control runs
			elseif isexp,
				%Prepare the ctrl values
				if strcmpi(expename,'exp01') | strcmpi(expename,'exp02') | strcmpi(expename,'exp03') | strcmpi(expename,'exp04') | strcmpi(expename,'exp11') | strcmpi(expename,'expA1') | strcmpi(expename,'expA2') | strcmpi(expename,'expA3') | strcmpi(expename,'expA4'),
					ctrl_file=['ComputedScalarsPaper/' group '/' simul '/ctrl_proj_open/computed_bmbfl_AIS_' group '_' simul '_ctrl_proj_open.nc'];
				elseif strcmpi(expename,'exp05') | strcmpi(expename,'exp06') | strcmpi(expename,'exp07') | strcmpi(expename,'exp08') | strcmpi(expename,'exp09') | strcmpi(expename,'exp10') | strcmpi(expename,'exp12') | strcmpi(expename,'exp13') | strcmpi(expename,'expA5') | strcmpi(expename,'expA6') | strcmpi(expename,'expA7') | strcmpi(expename,'expA8'),
					ctrl_file=['ComputedScalarsPaper/' group '/' simul '/ctrl_proj_std/computed_bmbfl_AIS_' group '_' simul '_ctrl_proj_std.nc'];
				else error('experiment not supported yet');
				end

				if exist(exp_directory,'dir')==0,
					error(['directory ' exp_directory ' not found']);
				end

				field='libmassbffl';
				if is_libmassbffl==1,
					exp_file=[exp_directory '/' field '_AIS_' group '_' simul '_' expename '.nc'];
					expmask_file=[exp_directory '/sftgif_AIS_' group '_' simul '_' expename '.nc'];
					expshelfmask_file=[exp_directory '/sftflf_AIS_' group '_' simul '_' expename '.nc'];
					data = double(ncread(exp_file,field)); 
					mask = double(ncread(expmask_file,'sftgif')); 
					shelfmask = double(ncread(expshelfmask_file,'sftflf')); 
					if(max(mask(:))>1 | min(mask(:))<0), disp(['mask should be between 0 and 1 in model ' modelname '']); end
					pos=find(mask>1);
					mask(pos)=1;
					pos=find(shelfmask>1);
					shelfmask(pos)=1;
					time = double(ncread(exp_file,'time'))/yearday_model;  
					if strcmpi(modelname,'VUW_PISM'),
						disp(['time should be in days in model ' modelname '']);
						time=time/(3600*24);
					end
					time=round(time,1);
					start_model_year=expe_model_year;
					time_vector(1:length(time))=start_model_year+time;
					%Forcings should be reported at the middle of the year so fix that for models that did not do it
					if strcmpi(modelname,'JPL1_ISSM') | strcmpi(modelname,'NCAR_CISM') | strcmpi(modelname,'PIK_PISM1') | strcmpi(modelname,'PIK_PISM2') | strcmpi(modelname,'UCIJPL_ISSM') | ...
							strcmpi(modelname,'ULB_FETISH32') | strcmpi(modelname,'ULB_FETISH16') | strcmpi(modelname,'UTAS_ElmerIce') | strcmpi(modelname,'VUB_AISMPALEO') | strcmpi(modelname,'VUW_PISM'),
						time_vector=time_vector-0.5;
					end
					if strcmpi(modelname,'ULB_FETISH32') | strcmpi(modelname,'ULB_FETISH16'),
						pos_time=find(time_vector>=2015 & time_vector<=2100);
					else
						pos_time=find(time_vector>=2015 & time_vector<=2101);
					end
					time_vector=time_vector(pos_time);
					for itime=1:length(time_vector),
						libmassbffl_i=data(:,:,itime);
						mask_i=mask(:,:,itime);
						shelfmask_i=shelfmask(:,:,itime);
						pos=find(mask_i==0);
						libmassbffl_i(pos)=0;
						shelfmask_i(pos)=0;
						posnan=find(isnan(mask_i));
						libmassbffl_i(posnan)=0;
						mask_i(posnan)=0;
						shelfmask_i(posnan)=0;
						posnan=find(isnan(libmassbffl_i));
						libmassbffl_i(posnan)=0;
						mask_i(posnan)=0;
						shelfmask_i(posnan)=0;
						posnan=find(isnan(shelfmask_i));
						libmassbffl_i(posnan)=0;
						mask_i(posnan)=0;
						shelfmask_i(posnan)=0;
						bmbfl_total=sum(libmassbffl_i(:).*shelfmask_i(:).*mask_i(:).*scalefac_model(:))*(resolution*1000)^2; %in kg/s
						if strcmpi(modelname,'VUW_PISM'),
							bmbfl_total=-bmbfl_total; %wrong sign in libmassbffl
						end
						bmbfl_vector(itime)=bmbfl_total;
						for isector=1:numsectors,
							eval(['sectors=sectors_' num2str(resolution) 'km;'])
							pos_sector=find(sectors==isector);
							bmbfl_total_sector=sum(libmassbffl_i(pos_sector).*shelfmask_i(pos_sector).*mask_i(pos_sector).*scalefac_model(pos_sector))*(resolution*1000)^2; %in kg/s
							eval(['bmbfl_total_sector_ctrlproj= double(ncread(ctrl_file,''bmbfl_sector_' num2str(isector) '''));'])
							eval(['bmbfl_vector_sector' num2str(isector) '(itime)=bmbfl_total_sector-bmbfl_total_sector_ctrlproj(itime);'])
						end
						for iregion=1:numregions,
							eval(['regions=regions_' num2str(resolution) 'km;'])
							pos_region=find(regions==iregion);
							bmbfl_total_region=sum(libmassbffl_i(pos_region).*shelfmask_i(pos_region).*mask_i(pos_region).*scalefac_model(pos_region))*(resolution*1000)^2; %in kg/s
							eval(['bmbfl_total_region_ctrlproj= double(ncread(ctrl_file,''bmbfl_region_' num2str(iregion) '''));'])
							eval(['bmbfl_vector_region' num2str(iregion) '(itime)=bmbfl_total_region-bmbfl_total_region_ctrlproj(itime);'])
						end
					end %end of time
					bmbfl_total_ctrlproj= double(ncread(ctrl_file,'bmbfl'));
					bmbfl_vector=bmbfl_vector-bmbfl_total_ctrlproj';

					if length(time)<85, error(['run is too short: ' int2str(length(time)) ' years']); end
					if (size(data,3)<85 | size(data,3)>122), error('field has the wrong size'); end
					if length(time)~=size(data,3), error(['lenght or time and data are not consistent in model ' modelname  ' for experiment ' exp_directory]); end
				end %end of lithk

				modelname
				size(start_model_year+time)
				size(time_vector)
				min(time_vector)
				max(time_vector)

				expbmbfl_file=['ComputedScalarsPaper/' group '/' simul '/' expename '/computed_bmbfl_minus_ctrl_proj_AIS_' group '_' simul '_' expename '.nc'];
				expbmbfl_file
				status=WriteNetCDFComputedOutputs(expbmbfl_file,'bmbfl','spatially integrated floating ice basal melt','kg/s',time_vector,bmbfl_vector,...
					bmbfl_vector_sector1,bmbfl_vector_sector2,bmbfl_vector_sector3,bmbfl_vector_sector4,bmbfl_vector_sector5,...
					bmbfl_vector_sector6,bmbfl_vector_sector7,bmbfl_vector_sector8,bmbfl_vector_sector9,bmbfl_vector_sector10,...
					bmbfl_vector_sector11,bmbfl_vector_sector12,bmbfl_vector_sector13,bmbfl_vector_sector14,bmbfl_vector_sector15,...
					bmbfl_vector_sector16,bmbfl_vector_sector17,bmbfl_vector_sector18,...
					bmbfl_vector_region1,bmbfl_vector_region2,bmbfl_vector_region3,...
					ice_density,ocean_density);

			end %end of isexp
		end %end of model
	end %end of experiment

end %}}}
