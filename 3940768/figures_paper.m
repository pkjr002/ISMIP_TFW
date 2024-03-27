%Routine to reproduce figures of manuscript https://tc.copernicus.org/articles/14/3033/2020/
%Contact: Helene Seroussi helene.seroussi@jpl.nasa.gov

step=1;
error('Add the correct paths where the scalars outputs and regridded datasets can be found on the lines below')
scalarpath='./ComputedScalarsPaper/'; %Change path for ISMIP6 regridded outputs
gridpath='./ISMIP6GriddedData/'; %Change path for ISMIP6 gridded outputs

model_list={'AWI_PISM1','DOE_MALI','ILTS_PIK_SICOPOLIS','IMAU_IMAUICE1','IMAU_IMAUICE2','JPL1_ISSM','LSCE_GRISLI','NCAR_CISM','PIK_PISM1','PIK_PISM2','UCIJPL_ISSM','ULB_FETISH32','ULB_FETISH16','UTAS_ElmerIce','VUB_AISMPALEO','VUW_PISM'};

model_list2=model_list;
for i=1:numel(model_list2)
	   model_list2{i} = strrep(model_list2{i},'_','\_');
end

%Figure 1: evolution of historical and unforced control 
if step==1 % {{{Figure 1a

	figure(1); set(gcf,'color','w'); set(gcf,'Position',[400 400 900 450]);
	colors = distinguishable_colors(length(model_list)*2);
	experiments_list={'ctrl_proj_std','ctrl_proj_open'};

	number=0;
	results_model={};
	max_acabf=-10; min_acabf=10^6;

	if 1,
		for imodel=1:length(model_list),
			modelname=model_list{imodel};

			for iexp=1:length(experiments_list),
				expename=experiments_list{iexp};
				specifics;

				eval(['isexp=is' expename ';'])
				if isexp,
					number=number+1;
					exptendacabf_file=['' scalarpath '/' group '/' simul '/' expename '/computed_smb_AIS_' group '_' simul '_' expename '.nc'];
					time_model=ncread(exptendacabf_file,'time');
					tendacabf_model=ncread(exptendacabf_file,'smb');
					if strcmpi(expename,'ctrl_proj_std') 
						if ishist_std,
							exptendacabf_file=['' scalarpath '/' group '/' simul '/hist_std/computed_smb_AIS_' group '_' simul '_hist_std.nc'];
							time_model=[ncread(exptendacabf_file,'time');time_model];
							tendacabf_model=[ncread(exptendacabf_file,'smb'); tendacabf_model];
							if strcmpi(modelname,'PIK_PISM1'), %Take every other step before 1950 to focus mostly on recent past
								pos=find(time_model<1950);
								time_model(1:length(pos)/2)=[];
								tendacabf_model(pos(1:2:end))=[];
							end
						end
						plot(time_model,tendacabf_model*yearday_model*3600*24/(10^9*1000),'color',colors(number,:)); hold on % from kg/s to Gt/yr
						results_model{end+1}=[model_list2{imodel} '\_std'] ;
					elseif strcmpi(expename,'ctrl_proj_open'),
						if ishist_open,
							exptendacabf_file=['' scalarpath '/' group '/' simul '/hist_open/computed_smb_AIS_' group '_' simul '_hist_open.nc'];
							time_model=[ncread(exptendacabf_file,'time');time_model];
							tendacabf_model=[ncread(exptendacabf_file,'smb');tendacabf_model];
						end
						plot(time_model,tendacabf_model*yearday_model*3600*24/(10^9*1000),'color',colors(number,:)); hold on % from kg/s to Gr/yr
						results_model{end+1}=[model_list2{imodel} '\_open'] ;
					end
				end

			end %end of model
		end %end of isexp

		x = [1900.5 2015 2015 1900.5];
		y = [2010 2010 3290 3290];
		text(1950,2050,'historical','VerticalAlignment','middle','HorizontalAlignment','center','color','k','fontweight','b');
		text(2055,2050,'future','VerticalAlignment','middle','HorizontalAlignment','center','color','k','fontweight','b');
		patch(x,y,[-1 -1 -1 -1],[0.9 0.9 0.9],'edgecolor',[0.9 0.9 0.9]); hold on
		xticks([1900 1925 1950 2000 2050 2100])
		xticklabels({'1850','1900','1950','2000','2050','2100'})
		xlim([1900 2100])
		ylim([2000 3300])
		xlabel('Time (yr)');
		ylabel('SMB (Gt/yr)');
		set(gcf,'Position',[400 400 900 450]);
		text(1880,1900,'a','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',18,'fontweight','b');
		hPlots = flip(findall(gcf,'Type','Line'));
		plot([1950 1950],[2010 3290],'--k');
		legend_str =results_model;
		legend(results_model,'location','EastOutside'); 
		legend boxoff
		h = gcf;
		set(h,'Units','Inches');
		pos = get(h,'Position');
		set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
		print(gcf, '-dpdf', '-painters', 'FiguresPaperFinal/Figure1a.pdf');
	end

end %}}}
if step==2 % {{{Figure 1b

	figure(1); set(gcf,'color','w'); set(gcf,'Position',[400 400 900 450]);
	colors = distinguishable_colors(length(model_list)*2);
	experiments_list={'ctrl_proj_std','ctrl_proj_open'};

	number=0;
	results_model={};
	%Ctr + hist special case
	if 1,
		for imodel=1:length(model_list),
			modelname=model_list{imodel};

			for iexp=1:length(experiments_list),
				expename=experiments_list{iexp};
				specifics;

				eval(['isexp=is' expename ';'])
				if isexp,
					number=number+1;
					exptendlibmassbffl_file=['' scalarpath '/' group '/' simul '/' expename '/computed_bmbfl_AIS_' group '_' simul '_' expename '.nc'];
					time_model=ncread(exptendlibmassbffl_file,'time');
					tendlibmassbffl_model=ncread(exptendlibmassbffl_file,'bmbfl');
					if strcmpi(expename,'ctrl_proj_std') 
						if ishist_std,
							exptendlibmassbffl_file=['' scalarpath '/' group '/' simul '/hist_std/computed_bmbfl_AIS_' group '_' simul '_hist_std.nc'];
							time_model=[ncread(exptendlibmassbffl_file,'time');time_model];
							tendlibmassbffl_model=[ncread(exptendlibmassbffl_file,'bmbfl'); tendlibmassbffl_model];
							if strcmpi(modelname,'PIK_PISM1'), %Take every other step before 1950 to focus mostly on recent past
								pos=find(time_model<1950);
								time_model(1:length(pos)/2)=[];
								tendlibmassbffl_model(pos(1:2:end))=[];
							end
						end
						plot(time_model,-tendlibmassbffl_model*yearday_model*3600*24/(10^9*1000),'color',colors(number,:)); hold on
						results_model{end+1}=[model_list2{imodel} '\_std'] ;
					elseif strcmpi(expename,'ctrl_proj_open'),
						if ishist_open,
							exptendlibmassbffl_file=['' scalarpath '/' group '/' simul '/hist_open/computed_bmbfl_AIS_' group '_' simul '_hist_open.nc'];
							time_model=[ncread(exptendlibmassbffl_file,'time');time_model];
							tendlibmassbffl_model=[ncread(exptendlibmassbffl_file,'bmbfl');tendlibmassbffl_model];
						end
						plot(time_model,-tendlibmassbffl_model*yearday_model*3600*24/(10^9*1000),'color',colors(number,:)); hold on
						results_model{end+1}=[model_list2{imodel} '\_open'] ;
					end
				end

			end %end of model
		end %end of isexp
		x = [1900.5 2015 2015 1900.5];
		y = [4480 4480 10 10];
		text(1950,100,'historical','VerticalAlignment','middle','HorizontalAlignment','center','color','k','fontweight','b');
		text(2055,100,'future','VerticalAlignment','middle','HorizontalAlignment','center','color','k','fontweight','b');
		xticks([1900 1925 1950 2000 2050 2100])
		xticklabels({'1850','1900','1950','2000','2050','2100'})
		patch(x,y,[-1 -1 -1 -1],[0.9 0.9 0.9],'edgecolor',[0.9 0.9 0.9]); hold on
		xlim([1900 2100])
		ylim([-10 4500])
		xlabel('Time (yr)');
		ylabel('Basal Melt (Gt/yr)');
		set(gcf,'Position',[400 400 900 450]);
		plot([1950 1950],[0 2290],'--k');
		text(1880,-350,'b','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',18,'fontweight','b');
		hPlots = flip(findall(gcf,'Type','Line'));
		plot([1950 1950],[0 4500],'--k');
		legend_str =results_model;
		legend(results_model,'location','EastOutside'); 
		legend boxoff
		h = gcf;
		set(h,'Units','Inches');
		pos = get(h,'Position');
		set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
		print(gcf, '-dpdf', '-painters', 'FiguresPaperFinal/Figure1b.pdf');
	end

end %}}}
if step==3 % {{{Figure 1c

	figure(1); set(gcf,'color','w'); set(gcf,'Position',[400 400 900 450]);
	colors = distinguishable_colors(length(model_list)*2);
	experiments_list={'ctrl_proj_std','ctrl_proj_open'};

	number=0;
	results_model={};
	%Ctr + hist special case
	if 1,
		for imodel=1:length(model_list),
			modelname=model_list{imodel};

			for iexp=1:length(experiments_list),
				expename=experiments_list{iexp};
				specifics;

				eval(['isexp=is' expename ';'])
				if isexp,
					number=number+1;
					explimnsw_file=['' scalarpath '/' group '/' simul '/' expename '/computed_ivaf_AIS_' group '_' simul '_' expename '.nc'];
					time_model=ncread(explimnsw_file,'time');
					limnsw_model=ncread(explimnsw_file,'ivaf');
					if strcmpi(expename,'ctrl_proj_std') 
						if ishist_std,
							explimnsw_file=['' scalarpath '/' group '/' simul '/hist_std/computed_ivaf_AIS_' group '_' simul '_hist_std.nc'];
							time_model=[ncread(explimnsw_file,'time');time_model];
							limnsw_model=[ncread(explimnsw_file,'ivaf'); limnsw_model];
							if strcmpi(modelname,'PIK_PISM1'), %Take every other step before 1950 to focus mostly on recent past
								pos=find(time_model<1950);
								time_model(1:length(pos)/2)=[];
								limnsw_model(pos(1:2:end))=[];
							end
						end
						plot(time_model,limnsw_model*ice_density/(10^9*1000),'color',colors(number,:)); hold on %from m^3 to Gt
						results_model{end+1}=[model_list2{imodel} '\_std'] ;
					elseif strcmpi(expename,'ctrl_proj_open'),
						if ishist_open,
							explimnsw_file=['' scalarpath '/' group '/' simul '/hist_open/computed_ivaf_AIS_' group '_' simul '_hist_open.nc'];
							time_model=[ncread(explimnsw_file,'time');time_model];
							limnsw_model=[ncread(explimnsw_file,'ivaf');limnsw_model];
						end
						plot(time_model,limnsw_model*ice_density/(10^9*1000),'color',colors(number,:)); hold on %from m^3 to Gt
						results_model{end+1}=[model_list2{imodel} '\_open'] ;
					end
				end

			end %end of model
		end %end of isexp
		x = [1900.5 2015 2015 1900.5];
		y = [1.981*10^7 1.981*10^7 2.159*10^7 2.159*10^7];
		text(1950,1.985*10^7,'historical','VerticalAlignment','middle','HorizontalAlignment','center','color','k','fontweight','b');
		text(2055,1.985*10^7,'future','VerticalAlignment','middle','HorizontalAlignment','center','color','k','fontweight','b');
		patch(x,y,[-1 -1 -1 -1],[0.9 0.9 0.9],'edgecolor',[0.9 0.9 0.9]); hold on
		h2=plot([1950 1950],[1.981*10^7 2.159*10^7],'--k');
		xticks([1900 1925 1950 2000 2050 2100])
		xticklabels({'1850','1900','1950','2000','2050','2100'})
		xlim([1900 2100])
		ylim([1.98 2.16]*10^7)
		xlabel('Time (yr)');
		ylabel('Ice Volume Above Floatation (Gt)');
		set(gcf,'Position',[400 400 900 450]);
		text(1880,1.97*10^7,'c','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',18,'fontweight','b');
		hPlots = flip(findall(gcf,'Type','Line'));
		legend_str =results_model;
		legend(results_model,'location','EastOutside'); 
		legend boxoff
		h = gcf;
		set(h,'Units','Inches');
		pos = get(h,'Position');
		set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
		print(gcf, '-dpdf', '-painters', 'FiguresPaperFinal/Figure1c.pdf');
	end

end %}}}

%Figure 2: ice and ice shelf extent
if step==4 % {{{Figure 2a

	experiments_list={'ctrl_proj_std','ctrl_proj_open'};

	numcases=0;
	ice_extent=zeros(761,761);

	for imodel=1:length(model_list),
		modelname=model_list{imodel};
		
		for iexp=1:length(experiments_list),
			expename=experiments_list{iexp};
			specifics;
			if strcmpi(expename,'ctrl_proj_std') 
				ctrl_directory=['' gridpath '/' group '/' simul '/' ctrl_proj_std_regrid '/'];
			elseif strcmpi(expename,'ctrl_proj_open') 
				ctrl_directory=['' gridpath '/' group '/' simul '/' ctrl_proj_open_regrid '/'];
			end

			eval(['isexp=is' expename ';'])
			if isexp,
				if is_sftgif==1,
					field='sftgif';
					if strcmpi(expename,'ctrl_proj_std') 
						ctrl_file=[ctrl_directory '/' field '_AIS_' group '_' simul '_ctrl_proj_std.nc'];
					elseif strcmpi(expename,'ctrl_proj_open') 
						ctrl_file=[ctrl_directory '/' field '_AIS_' group '_' simul '_ctrl_proj_open.nc'];
					end
					data = rot90(double(ncread(ctrl_file,field)));
					data_init=data(:,:,1);
					[data_nan]=find(isnan(data_init)); %Make sure there is no NaN data
					data_init(data_nan)=0;
					ice_extent=ice_extent+data_init;
					numcases=numcases+1;
				end
			end
		end
	end
	
	set(gcf,'color','w'); set(gcf,'Position',[400 400 700 500]);
	[pos_nani pos_nanj]=find(ice_extent==0);
	data_min=0; data_max=numcases;
	colorm   =flipud( parula(numcases));
	image_rgb = ind2rgb(uint16((max(data_min,min(data_max,ice_extent)) - data_min)*(size(colorm,1)/(data_max-data_min))),colorm);
	image_rgb(sub2ind(size(image_rgb),repmat(pos_nani,1,3),repmat(pos_nanj,1,3),repmat(1:3,size(pos_nani,1),1))) = repmat([1 1 1],size(pos_nani,1),1);
	imagesc(-3040:6080/size(ice_extent,2):3040,-3040:6080/size(ice_extent,1):3040,image_rgb); colorbar('off');
	set(gca,'fontsize',14)
	axis('equal','off'); xlim([-3040 3040]); ylim([-3040 3040]);
	caxis([data_min data_max]); colorbar
	text(-2400,2400,'a','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',18,'fontweight','b');
	x0=-1.70*10^3; y0=1.80*10^3;
	lengthscale=5*10^2; widthscale=2*10^1;
	patch([x0 x0+lengthscale x0+lengthscale x0],[y0 y0 y0+widthscale y0+widthscale],2*ones(1,4),'k','Edgecolor','k');
	text(x0,y0-10*widthscale,'500 km','fontsize',12)
	h = gcf;
	set(h,'Units','Inches');
	pos = get(h,'Position');
	set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
	print(gcf, '-dpdf', '-painters', 'FiguresPaperFinal/Figure2a.pdf');

end %}}}
if step==5 % {{{Figure 2b

	experiments_list={'ctrl_proj_std','ctrl_proj_open'};

	ice_extent=zeros(761,761);
	numcases=0;

	for imodel=1:length(model_list),
		modelname=model_list{imodel};
		
		for iexp=1:length(experiments_list),
			expename=experiments_list{iexp};
			specifics;
			if strcmpi(expename,'ctrl_proj_std') 
				ctrl_directory=['' gridpath '/' group '/' simul '/' ctrl_proj_std_regrid '/'];
			elseif strcmpi(expename,'ctrl_proj_open') 
				ctrl_directory=['' gridpath '/' group '/' simul '/' ctrl_proj_open_regrid '/'];
			end

			eval(['isexp=is' expename ';'])
			if isexp,
				if is_sftflf==1,
					field='sftflf';
					if strcmpi(expename,'ctrl_proj_std') 
						ctrl_file=[ctrl_directory '/' field '_AIS_' group '_' simul '_ctrl_proj_std.nc'];
						if strcmpi(modelname,'UCIJPL_ISSM'),
							mask_file=[ctrl_directory '/sftgif_AIS_' group '_' simul '_ctrl_proj_std.nc'];
						end
					elseif strcmpi(expename,'ctrl_proj_open') 
						ctrl_file=[ctrl_directory '/' field '_AIS_' group '_' simul '_ctrl_proj_open.nc'];
						if strcmpi(modelname,'UCIJPL_ISSM'),
							mask_file=[ctrl_directory '/sftgif_AIS_' group '_' simul '_ctrl_proj_open.nc'];
						end
					end
					data = rot90(double(ncread(ctrl_file,field)));
					data_init=data(:,:,1);
					if strcmpi(modelname,'UCIJPL_ISSM'),
						mask = rot90(double(ncread(mask_file,'sftgif')));
						mask_init=mask(:,:,1);
						data_init=data_init.*mask_init;
					end
					[data_nan]=find(isnan(data_init));
					data_init(data_nan)=0;
					ice_extent=ice_extent+data_init;
					numcases=numcases+1;
				end
			end
		end
	end
	
	set(gcf,'color','w'); set(gcf,'Position',[400 400 700 500]);
	[pos_nani pos_nanj]=find(ice_extent==0 | isnan(ice_extent));
	data_min=0; data_max=numcases;
	colorm   =flipud( parula(numcases));
	image_rgb = ind2rgb(uint16((max(data_min,min(data_max,ice_extent)) - data_min)*(size(colorm,1)/(data_max-data_min))),colorm);
	image_rgb(sub2ind(size(image_rgb),repmat(pos_nani,1,3),repmat(pos_nanj,1,3),repmat(1:3,size(pos_nani,1),1))) = repmat([1 1 1],size(pos_nani,1),1);
	imagesc(-3040:6080/size(ice_extent,2):3040,-3040:6080/size(ice_extent,1):3040,image_rgb); colorbar('off');
	set(gca,'fontsize',14)
	axis('equal','off'); xlim([-3040 3040]); ylim([-3040 3040]);
	caxis([data_min data_max]); colorbar
	text(-2400,2400,'b','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',18,'fontweight','b');
	x0=-1.70*10^3; y0=1.80*10^3;
	lengthscale=5*10^2; widthscale=2*10^1;
	patch([x0 x0+lengthscale x0+lengthscale x0],[y0 y0 y0+widthscale y0+widthscale],2*ones(1,4),'k','Edgecolor','k');
	text(x0,y0-10*widthscale,'500 km','fontsize',12)
	h = gcf;
	set(h,'Units','Inches');
	pos = get(h,'Position');
	set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
	print(gcf, '-dpdf', '-painters', 'FiguresPaperFinal/Figure2b.pdf');

end %}}}

%Figure 3: errors in vel and thickness 
if step==6 % {{{Figure 3a

	experiments_list={'ctrl_proj_std','ctrl_proj_open'};
	rms_thickness=[];
	thickness_grid=[];
	error('Need to create thickness grid observations on the 8 km standard ISMIP6-Antarctic grid, for example using BedMarchineAntarctica (https://nsidc.org/data/nsidc-0756) and the associated Matlab toolbox (https://www.mathworks.com/matlabcentral/fileexchange/69159-bedmachine) ');
	if size(thickness_grid,1)~=761 | size(thickness_grid,1)~=761, error('size of the thickness 8 km grid should be 761*761'); end
	number=0;
	results_model={};

	for imodel=1:length(model_list),
		thickness_mod=[];
		thickness_obs=[];
		modelname=model_list{imodel};
		
		for iexp=1:length(experiments_list),
			expename=experiments_list{iexp};
			specifics;

			eval(['isexp=is' expename ';'])
			if isexp,
				number=number+1;
				if strcmpi(expename,'ctrl_proj_std') 
					ctrl_directory=['' gridpath '/' group '/' simul '/' ctrl_proj_std_regrid '/'];
				elseif strcmpi(expename,'ctrl_proj_open') 
					ctrl_directory=['' gridpath '/' group '/' simul '/' ctrl_proj_open_regrid '/'];
				end

				field='lithk';
				if strcmpi(expename,'ctrl_proj_std') 
					ctrl_file=[ctrl_directory '/' field '_AIS_' group '_' simul '_ctrl_proj_std.nc'];
					results_model{end+1}=[model_list2{imodel} '\_std'] ;
				elseif strcmpi(expename,'ctrl_proj_open') 
					ctrl_file=[ctrl_directory '/' field '_AIS_' group '_' simul '_ctrl_proj_open.nc'];
					results_model{end+1}=[model_list2{imodel} '\_open'] ;
				end
				data = rot90(double(ncread(ctrl_file,field)));
				data_init=data(:,:,1);
				pos=find(data_init~=0 & ~isnan(data_init) & ~isnan(thickness_grid));
				thickness_mod=data_init(pos);
				thickness_obs=thickness_grid(pos);
				rms_thickness(number)=rms(thickness_obs-thickness_mod);
			end
		end
	end

	close all
	figure(1); set(gcf,'color','w'); set(gcf,'Position',[400 400 470 800]); hold on;
	ylim([0 number+1]);
	xlim([0 400]);
	for i=0:100:500,
		plot([i i],[0 number+1],':k')
	end
	b=barh(rms_thickness);  set(gca, 'XAxisLocation', 'top'); xlabel('RMSE Thickness (m)'); box('on');
	set(gca, 'YTickLabel', ''); set(gca,'Ydir','Reverse')
	colors = distinguishable_colors(number);
	for imodel=1:number,
		bi=barh(imodel, rms_thickness(imodel),'facecolor',colors(imodel,:));
		text(-8,imodel,results_model{imodel},'VerticalAlignment','middle','HorizontalAlignment','right','color',colors(imodel,:),'fontweight','b');
	 end
   set(gca,'position',[0.42 0.1 0.52 0.8])
	ylim([0 number+1]);
	text(-250,-1,'a','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',18,'fontweight','b');
	axis on

	h = gcf;
	set(h,'Units','Inches');
	pos = get(h,'Position');
	set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
	print(gcf, '-dpdf', '-painters', 'FiguresPaperFinal/Figure3a.pdf');
	
end %}}}
if step==7 % {{{Figure 3b

	experiments_list={'ctrl_proj_std','ctrl_proj_open'};
	rms_velocity=[];
	velocity_grid=[];
	error('Need to create velocity grid observations on the 8 km standard ISMIP6-Antarctic grid, for example using MEaSUREs Antarctic velocities https://nsidc.org/data/nsidc-0484');
	if size(velocity_grid,1)~=761 | size(velocity_grid,1)~=761, error('size of the velocity 8 km grid should be 761*761'); end
	velocity_grid=flipud(velocity_grid);
	number=0;
	results_model={};

	for imodel=1:length(model_list),
		velocity_mod=[];
		velocity_obs=[];
		modelname=model_list{imodel};

		for iexp=1:length(experiments_list),
			expename=experiments_list{iexp};
			specifics;

			eval(['isexp=is' expename ';'])
			if isexp,
				number=number+1;
				if strcmpi(expename,'ctrl_proj_std') 
					ctrl_directory=['' gridpath '/' group '/' simul '/' ctrl_proj_std_regrid '/'];
				elseif strcmpi(expename,'ctrl_proj_open') 
					ctrl_directory=['' gridpath '/' group '/' simul '/' ctrl_proj_open_regrid '/'];
				end

				if is_xvelsurf==1,
					field='xvelsurf'; field2='yvelsurf';
				elseif is_xvelmean==1, 
					field='xvelmean'; field2='yvelmean';
				else error('should have some velocity fields');
				end
				if strcmpi(expename,'ctrl_proj_std') 
					ctrl_file=[ctrl_directory '/' field '_AIS_' group '_' simul '_ctrl_proj_std.nc'];
					ctrl_file2=[ctrl_directory '/' field2 '_AIS_' group '_' simul '_ctrl_proj_std.nc'];
					results_model{end+1}=[model_list2{imodel} '\_std'] ;
				elseif strcmpi(expename,'ctrl_proj_open') 
					ctrl_file=[ctrl_directory '/' field '_AIS_' group '_' simul '_ctrl_proj_open.nc'];
					ctrl_file2=[ctrl_directory '/' field2 '_AIS_' group '_' simul '_ctrl_proj_open.nc'];
					results_model{end+1}=[model_list2{imodel} '\_open'] ;
				end
				datau = rot90(double(ncread(ctrl_file,field)));
				datav = rot90(double(ncread(ctrl_file2,field2)));
				data=sqrt(datau.^2+datav.^2)*31556926;
				if size(data)~=[761,761,21];
					error(['warming: file ' ctrl_file ' has the wrong size']);
				end
				data_init=data(:,:,1);
				pos=find(data_init~=0 & ~isnan(data_init) & ~isnan(velocity_grid));
				velocity_mod=data_init(pos);
				velocity_obs=velocity_grid(pos);
				rms_velocity(number)=rms(velocity_obs-velocity_mod);
				rms_velocity_log(number)=rms(log(velocity_obs)-log(velocity_mod));
			end

		end
	end

	figure(1); set(gcf,'color','w'); set(gcf,'Position',[400 400 450 800]); hold on;
	for i=0:100:500,
		plot([i i],[0 number+1],':k')
	end
	b=barh(rms_velocity);  set(gca, 'XAxisLocation', 'top'); xlabel('RMSE Velocity (m/yr)'); box('on'); xlim([0 400]);
	set(gca, 'YTickLabel', ''); set(gca,'Ydir','Reverse')
	colors = distinguishable_colors(number+10);
	for imodel=1:number,
		bi=barh(imodel, rms_velocity(imodel),'facecolor',colors(imodel,:));
		text(-10,imodel,results_model{imodel},'VerticalAlignment','middle','HorizontalAlignment','right','color',colors(imodel,:),'fontweight','b');
	end
	set(gca,'position',[0.43 0.1 0.52 0.8])
	ylim([0 number+1]);
	xlim([0 500]);
	text(-335,-1,'b','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',18,'fontweight','b');

	h = gcf;
	set(h,'Units','Inches');
	pos = get(h,'Position');
	set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
	print(gcf, '-dpdf', '-painters', 'FiguresPaperFinal/Figure3b.pdf');

end %}}}

%Figures 4: NorESM open and standard
if step==8 % {{{Figure 4

	figure(1); set(gcf,'color','w'); set(gcf,'Position',[400 500 800 450]);
	colors = distinguishable_colors(length(model_list)+10);
	experiments_list={'exp05','exp01'}; %open and standard NorESM

	number=0;
	results_model={};

	if 1,
		for imodel=1:length(model_list),
			modelname=model_list{imodel};
			for iexp=1:length(experiments_list),
				expename=experiments_list{iexp};
				specifics;

				eval(['isexp=is' expename ';'])
				if isexp,
					number=number+1;
					explimnsw_file=['' scalarpath '/' group '/' simul '/' expename '/computed_ivaf_minus_ctrl_proj_AIS_' group '_' simul '_' expename '.nc'];
					time_model=ncread(explimnsw_file,'time');
					limnsw_model=ncread(explimnsw_file,'ivaf');
					if strcmpi(expename,'exp05'),
						results_model{end+1}=[model_list2{imodel} '\_std'] ;
					elseif strcmpi(expename,'exp01'),
						results_model{end+1}=[model_list2{imodel} '\_open'] ;
					end
					plot([2015;time_model],[0;-limnsw_model/362.5*ice_density/(10^9*1000)],'color',colors(number,:)); hold on % in mm SLE
				end

			end %end of model
		end %end of isexp
		legend(results_model,'location','EastOutside'); 
		legend boxoff
		xlim([2015 2100])
		ylim([-50 200])
		xlabel('Time (yr)');
		ylabel('Sea Level Contribution (mm SLE)');
		set(gcf,'Position',[400 500 800 450]);
		h = gcf;
		set(h,'Units','Inches');
		pos = get(h,'Position');
		set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
		print(gcf, '-dpdf', '-painters', 'FiguresPaperFinal/Figure4.pdf');
	end

end %}}}

%Figure 5: NorESM by region
if step==9 % {{{Figure 5

	colors = distinguishable_colors(length(model_list)*2);
	experiments_list={'exp05','exp01'}; %open and standard NorESM

	results_west=[]; results_east=[]; results_penin=[]; 
	resultstendacabfgr_west=[]; resultstendacabfgr_east=[]; resultstendacabfgr_penin=[]; 
	results_model={}; results_meltparam={};

	if 1,
		for imodel=1:length(model_list),
			modelname=model_list{imodel};
			for iexp=1:length(experiments_list),
				expename=experiments_list{iexp};
				specifics;

				eval(['isexp=is' expename ';'])
				if isexp,
					explimnsw_file=['' scalarpath '/' group '/' simul '/' expename '/computed_ivaf_minus_ctrl_proj_AIS_' group '_' simul '_' expename '.nc'];
					exptendacabfgr_file=['' scalarpath '/' group '/' simul '/' expename '/computed_smbgr_minus_ctrl_proj_AIS_' group '_' simul '_' expename '.nc'];
					time_model=ncread(explimnsw_file,'time');
					limnsw_west=ncread(explimnsw_file,'ivaf_region_1');
					limnsw_east=ncread(explimnsw_file,'ivaf_region_2');
					limnsw_penin=ncread(explimnsw_file,'ivaf_region_3');
					tendacabfgr_west=ncread(exptendacabfgr_file,'smbgr_region_1');
					tendacabfgr_east=ncread(exptendacabfgr_file,'smbgr_region_2');
					tendacabfgr_penin=ncread(exptendacabfgr_file,'smbgr_region_3');
					results_west(end+1)=limnsw_west(85);
					results_east(end+1)=limnsw_east(85);
					results_penin(end+1)=limnsw_penin(85);
					resultstendacabfgr_west(end+1)=sum(tendacabfgr_west(1:85));
					resultstendacabfgr_east(end+1)=sum(tendacabfgr_east(1:85));
					resultstendacabfgr_penin(end+1)=sum(tendacabfgr_penin(1:85));
					if strcmpi(expename,'exp05'),
						results_model{end+1}=[model_list2{imodel} '\_std'] ;
						results_meltparam{end+1}='standard';
					elseif strcmpi(expename,'exp01'),
						results_model{end+1}=[model_list2{imodel} '\_open'] ;
						results_meltparam{end+1}='open';
					end
				end
			end %end of model
		end %end of isexp

		results=[results_west;results_east;results_penin];
		%Plot results
		figure(1); set(gcf,'color','w'); set(gcf,'Position',[400 500 900 400]);
		hb=bar(results*(-1/362.5)*ice_density/(10^9*1000),'grouped'); hold on  %in mm SLE
		for iresults=1:length(results_model),
			hb(iresults).FaceColor=colors(iresults,:);
		end
		set(gca, 'XTickLabel', '');
		h=text(1,-48,'WAIS','VerticalAlignment','middle','HorizontalAlignment','center','color','k','fontweight','b');
		h=text(2,-48,'EAIS','VerticalAlignment','middle','HorizontalAlignment','center','color','k','fontweight','b');
		h=text(3,-48,'Peninsula','VerticalAlignment','middle','HorizontalAlignment','center','color','k','fontweight','b');
		ylabel('Sea Level Contribution (mm SLE)');
		xlim([0.5 3.5]);
		ylim([-40 160]);
		hold on

		%Add SMB results
		resultstendacabfgr=[resultstendacabfgr_west,resultstendacabfgr_east,resultstendacabfgr_penin];
		for ib = 1:numel(hb)
			%XData property is the tick labels/group centers; XOffset is the offset of each distinct group
			xData = hb(ib).XData+hb(ib).XOffset;
			plot(xData(1),resultstendacabfgr_west(ib)*(-1/362.5)*yearday_model*3600*24/(10^9*1000),'d','MarkerFaceColor',colors(ib,:),'MarkerSize',6,'MarkerEdgeColor','k') %from kg/s to mm SLE
			plot(xData(2),resultstendacabfgr_east(ib)*(-1/362.5)*yearday_model*3600*24/(10^9*1000),'d','MarkerFaceColor',colors(ib,:),'MarkerSize',6,'MarkerEdgeColor','k')
			plot(xData(3),resultstendacabfgr_penin(ib)*(-1/362.5)*yearday_model*3600*24/(10^9*1000),'d','MarkerFaceColor',colors(ib,:),'MarkerSize',6,'MarkerEdgeColor','k')
		end
		legend(results_model,'location','EastOutside')

		h = gcf;
		set(h,'Units','Inches');
		pos = get(h,'Position');
		set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
		print(gcf, '-dpdf', '-painters', 'FiguresPaperFinal/Figure5.pdf');
	end

end %}}}

%Figure 6: DeltaH and DeltaVel NorESM 
if step==10 % {{{Figure 6a

	nummodels=zeros(761,761);
	totalthicknesschange=zeros(761,761);
	thicknessstd=zeros(761,761);

	experiments_list={'exp01','exp05'}; %NorESM RCP8.5

	for iexp=1:length(experiments_list),
		expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			modelname=model_list{imodel};
			specifics;
			eval(['exp_name=' expename '_regrid;'])
			exp_directory=['' gridpath '/' group '/' simul '/' exp_name '/'];
			if strcmp(expename,'exp01'),
				ctrl_directory=['' gridpath '/' group '/' simul '/' ctrl_proj_open_regrid '/'];
			elseif strcmp(expename,'exp05'),
				ctrl_directory=['' gridpath '/' group '/' simul '/' ctrl_proj_std_regrid '/'];
			end

			eval(['isexp=is' expename ';'])
			if isexp,
				if exist(exp_directory,'dir')==0,
					error(['directory ' exp_directory ' not found']);
				end

				if is_lithk==1,
					field='lithk';
					exp_file=[exp_directory '/' field '_AIS_' group '_' simul '_' expename '.nc'];
					mask_file=[exp_directory '/sftgif_AIS_' group '_' simul '_' expename '.nc'];
					if strcmp(expename,'exp01'),
						ctrl_file=[ctrl_directory '/' field '_AIS_' group '_' simul '_ctrl_proj_open.nc'];
						maskctrl_file=[ctrl_directory '/sftgif_AIS_' group '_' simul '_ctrl_proj_open.nc'];
					elseif strcmp(expename,'exp05'),
						ctrl_file=[ctrl_directory '/' field '_AIS_' group '_' simul '_ctrl_proj_std.nc'];
						maskctrl_file=[ctrl_directory '/sftgif_AIS_' group '_' simul '_ctrl_proj_std.nc'];
					end

					data = rot90(double(ncread(exp_file,field)));
					datac = rot90(double(ncread(ctrl_file,field)));
					mask = rot90(double(ncread(mask_file,'sftgif')));
					maskc = rot90(double(ncread(maskctrl_file,'sftgif')));
					if size(data)~=[761,761,21];
						error(['warming: file ' exp_file ' has the wrong size']);
					end
					data_init=data(:,:,1).*mask(:,:,1);
					data_end=data(:,:,86).*mask(:,:,86);
					datac_init=datac(:,:,1).*maskc(:,:,1);
					datac_end=datac(:,:,86).*maskc(:,:,86);
					pos=find(data_init~=0 & ~isnan(data_init) & data_end~=0 & ~isnan(data_end) & datac_init~=0 & ~isnan(datac_init) & datac_end~=0 & ~isnan(datac_end));
					nummodels(pos)=nummodels(pos)+1;
					totalthicknesschange(pos)=totalthicknesschange(pos)+data_end(pos)-data_init(pos)-(datac_end(pos)-datac_init(pos));

				end
			end
		end
	end
	mean_thicknesschange=totalthicknesschange./nummodels;
	pos=find(nummodels<5);
	mean_thicknesschange(pos)=NaN;

	if 1,
		set(gcf,'color','w');
		[pos_nani pos_nanj]=find(isnan(mean_thicknesschange));
		data_min=-100; data_max=100;
		colorm = jet(100);
		%colormap used in the paper can be found here: https://www.mathworks.com/matlabcentral/fileexchange/17555-light-bartlein-color-maps
		image_rgb = ind2rgb(uint16((max(data_min,min(data_max,mean_thicknesschange)) - data_min)*(size(colorm,1)/(data_max-data_min))),colorm);
		image_rgb(sub2ind(size(image_rgb),repmat(pos_nani,1,3),repmat(pos_nanj,1,3),repmat(1:3,size(pos_nani,1),1))) = repmat([1 1 1],size(pos_nani,1),1);
		imagesc(-3040:6080/size(mean_thicknesschange,2):3040,-3040:6080/size(mean_thicknesschange,1):3040,image_rgb);
		set(gca,'fontsize',14); caxis([data_min data_max]); colormap(jet); hcb=colorbar; title(hcb,'m');
		axis('equal','off'); xlim([-3040 3040]); ylim([-3040 3040]);
		text(-2400,2400,'a','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',18,'fontweight','b');

		h = gcf;
		set(h,'Units','Inches');
		pos = get(h,'Position');
		set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
		print(gcf, '-dpdf', '-painters', 'FiguresPaperFinal/Figure6a.pdf');
	end

	for iexp=1:length(experiments_list),
		expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			modelname=model_list{imodel};
			specifics;
			eval(['exp_name=' expename '_regrid;'])
			exp_directory=[path '/' exp_name '/'];
			exp_directory=['' gridpath '/' group '/' simul '/' exp_name '/'];
			if strcmp(expename,'exp01'),
				ctrl_directory=['' gridpath '/' group '/' simul '/' ctrl_proj_open_regrid '/'];
			elseif strcmp(expename,'exp05'),
				ctrl_directory=['' gridpath '/' group '/' simul '/' ctrl_proj_std_regrid '/'];
			end

			eval(['isexp=is' expename ';'])
			if isexp,
				if exist(exp_directory,'dir')==0,
					error(['directory ' exp_directory ' not found']);
				end

				if is_lithk==1,
					field='lithk';
					exp_file=[exp_directory '/' field '_AIS_' group '_' simul '_' expename '.nc'];
					if strcmp(expename,'exp01'),
						ctrl_file=[ctrl_directory '/' field '_AIS_' group '_' simul '_ctrl_proj_open.nc'];
					elseif strcmp(expename,'exp05'),
						ctrl_file=[ctrl_directory '/' field '_AIS_' group '_' simul '_ctrl_proj_std.nc'];
					end

					data = rot90(double(ncread(exp_file,field)));
					datac = rot90(double(ncread(ctrl_file,field)));
					if size(data)~=[761,761,21];
						error(['warming: file ' exp_file ' has the wrong size']);
					end
					data_init=data(:,:,end-85);
					data_end=data(:,:,end);
					datac_init=datac(:,:,end-85);
					datac_end=datac(:,:,end);
					pos=find(data_init~=0 & ~isnan(data_init) & data_end~=0 & ~isnan(data_end) & datac_init~=0 & ~isnan(datac_init) & datac_end~=0 & ~isnan(datac_end));
					thicknessstd(pos)=thicknessstd(pos)+((data_end(pos)-data_init(pos))-(datac_end(pos)-datac_init(pos))-mean_thicknesschange(pos)).^2;

				end
			end
		end
	end
	thicknessstd=sqrt(thicknessstd./(nummodels-1));
	pos=find(nummodels<5);
	thicknessstd(pos)=NaN;

	close; set(gcf,'color','w'); set(gcf,'Position',[400 400 700 500]);
	[pos_nani pos_nanj]=find(isnan(thicknessstd));
	data_min=0; data_max=200;
	colorm = flipud(hot(100));
	image_rgb = ind2rgb(uint16((max(data_min,min(data_max,thicknessstd)) - data_min)*(size(colorm,1)/(data_max-data_min))),colorm);
	image_rgb(sub2ind(size(image_rgb),repmat(pos_nani,1,3),repmat(pos_nanj,1,3),repmat(1:3,size(pos_nani,1),1))) = repmat([1 1 1],size(pos_nani,1),1);
	imagesc(-3040:6080/size(thicknessstd,2):3040,-3040:6080/size(thicknessstd,1):3040,image_rgb); 
	set(gca,'fontsize',14); colormap(flipud(hot));  hcb=colorbar; title(hcb,'m'); 
	axis('equal','off'); caxis([data_min data_max]); xlim([-3040 3040]); ylim([-3040 3040]);
	text(-2400,2400,'c','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',18,'fontweight','b');

	h = gcf;
	set(h,'Units','Inches');
	pos = get(h,'Position');
	set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
	print(gcf, '-dpdf', '-painters', 'FiguresPaperFinal/Figure6c.pdf');
	
end %}}}
if step==11 % {{{Figure 6b

	experiments_list={'exp01','exp05'}; %NorESM RCP8.5

	nummodels=zeros(761,761);
	totalvelocitychange=zeros(761,761);
	velocitystd=zeros(761,761);

	for iexp=1:length(experiments_list),
		expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			modelname=model_list{imodel};
			specifics;
			eval(['exp_name=' expename '_regrid;'])
			exp_directory=['' gridpath '/' group '/' simul '/' exp_name '/'];
			if strcmp(expename,'exp01'),
				ctrl_directory=['' gridpath '/' group '/' simul '/' ctrl_proj_open_regrid '/'];
			elseif strcmp(expename,'exp05'),
				ctrl_directory=['' gridpath '/' group '/' simul '/' ctrl_proj_std_regrid '/'];
			end

			eval(['isexp=is' expename ';'])
			if isexp,
				if exist(exp_directory,'dir')==0,
					error(['directory ' exp_directory ' not found']);
				end

				if is_xvelmean==1, 
					field='xvelmean'; field2='yvelmean';
				else 
					error('should have some velocity fields');
				end
				exp_file=[exp_directory '/' field '_AIS_' group '_' simul '_' expename '.nc'];
				exp_file2=[exp_directory '/' field2 '_AIS_' group '_' simul '_' expename '.nc'];
				if strcmp(expename,'exp01'),
					ctrl_file=[ctrl_directory '/' field '_AIS_' group '_' simul '_ctrl_proj_open.nc'];
					ctrl_file2=[ctrl_directory '/' field2 '_AIS_' group '_' simul '_ctrl_proj_open.nc'];
				elseif strcmp(expename,'exp05'),
					ctrl_file=[ctrl_directory '/' field '_AIS_' group '_' simul '_ctrl_proj_std.nc'];
					ctrl_file2=[ctrl_directory '/' field2 '_AIS_' group '_' simul '_ctrl_proj_std.nc'];
				end

				datau = rot90(double(ncread(exp_file,field)));
				datav = rot90(double(ncread(exp_file2,field2)));
				data=sqrt(datau.^2+datav.^2)*31556926; %in m/yr
				datacu = rot90(double(ncread(ctrl_file,field)));
				datacv = rot90(double(ncread(ctrl_file2,field2)));
				datac=sqrt(datacu.^2+datacv.^2)*31556926; % in m/yr
				if size(data)~=[761,761,21];
					error(['warming: file ' abmb_file ' has the wrong size']);
				end
				data_init=data(:,:,1);
				data_end=data(:,:,end);
				datac_init=datac(:,:,1);
				datac_end=datac(:,:,end);
				[data_nan]=find(isnan(data_init)); data_init(data_nan)=0;
				[data_nan]=find(isnan(data_end)); data_end(data_nan)=0;
				[datac_nan]=find(isnan(datac_init)); datac_init(datac_nan)=0;
				[datac_nan]=find(isnan(datac_end)); datac_end(datac_nan)=0;
				pos=find(datac_init~=0 & ~isnan(datac_init) & datac_end~=0 & ~isnan(datac_end) & datac_init~=0 & ~isnan(datac_init) & datac_end~=0 & ~isnan(datac_end)); 
				nummodels(pos)=nummodels(pos)+1;
				totalvelocitychange(pos)=totalvelocitychange(pos)+(data_end(pos)-data_init(pos)) - (datac_end(pos)-datac_init(pos));

			end
		end
	end
	mean_velocitychange=totalvelocitychange./nummodels;
	pos=find(nummodels<5);
	mean_velocitychange(pos)=NaN;

	if 1,
		set(gcf,'color','w');
		[pos_nani pos_nanj]=find(isnan(mean_velocitychange));
		data_min=-100; data_max=100;
		%colormap used in the paper can be found here: https://www.mathworks.com/matlabcentral/fileexchange/17555-light-bartlein-color-maps
		colorm = jet(100);
		image_rgb = ind2rgb(uint16((max(data_min,min(data_max,mean_velocitychange)) - data_min)*(size(colorm,1)/(data_max-data_min))),colorm);
		image_rgb(sub2ind(size(image_rgb),repmat(pos_nani,1,3),repmat(pos_nanj,1,3),repmat(1:3,size(pos_nani,1),1))) = repmat([1 1 1],size(pos_nani,1),1);
		imagesc(-3040:6080/size(mean_velocitychange,2):3040,-3040:6080/size(mean_velocitychange,1):3040,image_rgb);
		axis('equal','off'); xlim([-3040 3040]); ylim([-3040 3040]);
		colormap(jet); caxis([data_min data_max]); hcb=colorbar; title(hcb,'m yr^{-1}');
		text(-2400,2400,'b','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',18,'fontweight','b');

		h = gcf;
		set(h,'Units','Inches');
		pos = get(h,'Position');
		set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
		print(gcf, '-dpdf', '-painters', 'FiguresPaperFinal/Figure6b.pdf');
	end

	for iexp=1:length(experiments_list),
		expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			modelname=model_list{imodel};
			specifics;
			eval(['exp_name=' expename '_regrid;'])
			exp_directory=['' gridpath '/' group '/' simul '/' exp_name '/'];
			if strcmp(expename,'exp01'),
				ctrl_directory=['' gridpath '/' group '/' simul '/' ctrl_proj_open_regrid '/'];
			elseif strcmp(expename,'exp05'),
				ctrl_directory=['' gridpath '/' group '/' simul '/' ctrl_proj_std_regrid '/'];
			end

			eval(['isexp=is' expename ';'])
			if isexp,
				if exist(exp_directory,'dir')==0,
					error(['directory ' exp_directory ' not found']);
				end

				if is_xvelmean==1, 
					field='xvelmean'; field2='yvelmean';
				else error('should have some velocity fields');
				end
				exp_file=[exp_directory '/' field '_AIS_' group '_' simul '_' expename '.nc'];
				exp_file2=[exp_directory '/' field2 '_AIS_' group '_' simul '_' expename '.nc'];
				mask_file2=[exp_directory '/sftgif_AIS_' group '_' simul '_' expename '.nc'];
				if strcmp(expename,'exp01'),
						ctrl_file=[ctrl_directory '/' field '_AIS_' group '_' simul '_ctrl_proj_open.nc'];
						ctrl_file2=[ctrl_directory '/' field2 '_AIS_' group '_' simul '_ctrl_proj_open.nc'];
						maskctrl_file=[ctrl_directory '/sftgif_AIS_' group '_' simul '_ctrl_proj_open.nc'];
				elseif strcmp(expename,'exp05'),
					ctrl_file=[ctrl_directory '/' field '_AIS_' group '_' simul '_ctrl_proj_std.nc'];
					ctrl_file2=[ctrl_directory '/' field2 '_AIS_' group '_' simul '_ctrl_proj_std.nc'];
						maskctrl_file=[ctrl_directory '/sftgif_AIS_' group '_' simul '_ctrl_proj_std.nc'];
				end

				datau = rot90(double(ncread(exp_file,field)));
				datav = rot90(double(ncread(exp_file2,field2)));
				data=sqrt(datau.^2+datav.^2)*31556926;
				datacu = rot90(double(ncread(ctrl_file,field)));
				datacv = rot90(double(ncread(ctrl_file2,field2)));
				datac=sqrt(datacu.^2+datacv.^2)*31556926;
				if size(data)~=[761,761,21];
					error(['warming: file ' abmb_file ' has the wrong size']);
				end
				data_init=data(:,:,1); %.*mask(:,:,1);
				data_end=data(:,:,end); %.*mask(:,:,end);
				datac_init=datac(:,:,1); %.*maskctrl(:,:,1);
				datac_end=datac(:,:,end); %.*maskctrl(:,:,end);
				[data_nan]=find(isnan(data_init)); data_init(data_nan)=0;
				[data_nan]=find(isnan(data_end)); data_end(data_nan)=0;
				[datac_nan]=find(isnan(datac_init)); datac_init(datac_nan)=0;
				[datac_nan]=find(isnan(datac_end)); datac_end(datac_nan)=0;
				pos=find(datac_init~=0 & ~isnan(datac_init) & datac_end~=0 & ~isnan(datac_end) & datac_init~=0 & ~isnan(datac_init) & datac_end~=0 & ~isnan(datac_end)); 
				velocitystd(pos)=velocitystd(pos)+(data_end(pos)-data_init(pos) -(datac_end(pos)-datac_init(pos)) -mean_velocitychange(pos)).^2;
			end
		end
	end

	velocitystd=sqrt(velocitystd./(nummodels));
	pos=find(nummodels<5);

	close; set(gcf,'color','w'); set(gcf,'Position',[400 400 700 500]);
	[pos_nani pos_nanj]=find(isnan(velocitystd));
	data_min=0; data_max=200;
	colorm = flipud(hot(100));
	image_rgb = ind2rgb(uint16((max(data_min,min(data_max,velocitystd)) - data_min)*(size(colorm,1)/(data_max-data_min))),colorm);
	image_rgb(sub2ind(size(image_rgb),repmat(pos_nani,1,3),repmat(pos_nanj,1,3),repmat(1:3,size(pos_nani,1),1))) = repmat([1 1 1],size(pos_nani,1),1);
	imagesc(-3040:6080/size(velocitystd,2):3040,-3040:6080/size(velocitystd,1):3040,image_rgb);
	axis('equal','off'); xlim([-3040 3040]); ylim([-3040 3040]);
	colormap(flipud(hot)); caxis([data_min data_max]); hcb=colorbar; title(hcb,'m yr^{-1}');
	text(-2400,2400,'d','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',18,'fontweight','b');

	h = gcf;
	set(h,'Units','Inches');
	pos = get(h,'Position');
	set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
	print(gcf, '-dpdf', '-painters', 'FiguresPaperFinal/Figure6d.pdf');
end %}}}

%Figure 7: all RCP8.5 evolution 
if step==12 % {{{Figure 7

	figure(1); set(gcf,'color','w'); set(gcf,'Position',[400 500 800 450]);
	colors_order = get(gca,'colororder');
	colors_long_list = [colors_order(1:3,:);colors_order(1:3,:);colors_order(4:6,:);colors_order(4:6,:)];
	colors = colors_long_list + 0.75*(1-colors_long_list);
	colors_exp = colors_order(1:6,:);

	experiments_list={'exp01','exp02','exp04','exp05','exp06','exp08','expA1','expA2','expA3','expA5','expA6','expA7'}; %open and standard RCP8.5

	for i=1:6,
		plot([0 1],[0 1],'color',colors_exp(i,:),'linewidth',1); hold on
	end
	mean_vaf=zeros(86,length(experiments_list));
	min_vaf=9999*ones(86,length(experiments_list));
	max_vaf=-9999*ones(86,length(experiments_list));
	num_models=zeros(length(experiments_list),1);
	if 1,

		for iexp=1:length(experiments_list),
			expename=experiments_list{iexp};

			for imodel=1:length(model_list),
				modelname=model_list{imodel};
				specifics;

				eval(['isexp=is' expename ';'])
				if isexp,
					num_models(iexp)=num_models(iexp)+1;
					explimnsw_file=['' scalarpath '/' group '/' simul '/' expename '/computed_ivaf_minus_ctrl_proj_AIS_' group '_' simul '_' expename '.nc'];
					time_model=ncread(explimnsw_file,'time');
					limnsw_model=ncread(explimnsw_file,'ivaf')*ice_density/(10^9*1000); %in Gt
					mean_vaf(:,iexp)=mean_vaf(:,iexp)+[0;limnsw_model(1:85)];
					min_vaf(:,iexp)=min(min_vaf(:,iexp),[0;limnsw_model(1:85)]);
					max_vaf(:,iexp)=max(max_vaf(:,iexp),[0;limnsw_model(1:85)]);
					limnsw_model=[0;limnsw_model];
					plot([2015;time_model(1:85)],[0;-limnsw_model(1:85)/362.5],'color',colors(iexp,:),'linewidth',1); hold on %in mm SLE
				else
				end

			end %end of model
		end %end of isexp
		mean_exp=zeros(86,6);
		min_exp=zeros(86,6);
		max_exp=zeros(86,6);
		exp_models=zeros(6,1);
		mean_exp(:,1)=mean_vaf(:,1)+mean_vaf(:,4); exp_models(1)=num_models(1)+num_models(4);
		mean_exp(:,2)=mean_vaf(:,2)+mean_vaf(:,5); exp_models(2)=num_models(2)+num_models(5);
		mean_exp(:,3)=mean_vaf(:,3)+mean_vaf(:,6); exp_models(3)=num_models(3)+num_models(6);
		mean_exp(:,4)=mean_vaf(:,7)+mean_vaf(:,10); exp_models(4)=num_models(7)+num_models(10);
		mean_exp(:,5)=mean_vaf(:,8)+mean_vaf(:,11); exp_models(5)=num_models(8)+num_models(11);
		mean_exp(:,6)=mean_vaf(:,9)+mean_vaf(:,12); exp_models(6)=num_models(9)+num_models(12);
		mean_exp=mean_exp./exp_models';
		for i=1:6,
			plot([2015;time_model(1:85)],-mean_exp(:,i)/362.5,'color',colors_exp(i,:),'linewidth',2); hold on
		end
		ax = gca; ax.Clipping = 'off';
		x = [2100.5 2100.5 2101 2101];
		y = [min(min_vaf(85,1),min_vaf(85,4))  max(max_vaf(85,1),max_vaf(85,4)) max(max_vaf(85,1),max_vaf(85,4)) min(min_vaf(85,1),min_vaf(85,4))]/-362.5;
		patch(x,y,[-1 -1 -1 -1],colors(1,:),'edgecolor',colors(1,:)); hold on
		y = mean_exp(85,1)/-362.5+[0.4 -0.4 -0.4 0.4];
		patch(x,y,[-1 -1 -1 -1],colors_exp(1,:),'edgecolor',colors_exp(1,:)); hold on
		x = [2101.5 2101.5 2102 2102];
		y = [min(min_vaf(85,2),min_vaf(85,5))  max(max_vaf(85,2),max_vaf(85,5)) max(max_vaf(85,2),max_vaf(85,5)) min(min_vaf(85,2),min_vaf(85,5))]/-362.5;
		patch(x,y,[-1 -1 -1 -1],colors(2,:),'edgecolor',colors(2,:)); hold on
		y = mean_exp(85,2)/-362.5+[0.4 -0.4 -0.4 0.4];
		patch(x,y,[-1 -1 -1 -1],colors_exp(2,:),'edgecolor',colors_exp(2,:)); hold on
		x = [2102.5 2102.5 2103 2103];
		y = [min(min_vaf(85,3),min_vaf(85,6))  max(max_vaf(85,3),max_vaf(85,6)) max(max_vaf(85,3),max_vaf(85,6)) min(min_vaf(85,3),min_vaf(85,6))]/-362.5;
		patch(x,y,[-1 -1 -1 -1],colors(3,:),'edgecolor',colors(3,:)); hold on
		y = mean_exp(85,3)/-362.5+[0.4 -0.4 -0.4 0.4];
		patch(x,y,[-1 -1 -1 -1],colors_exp(3,:),'edgecolor',colors_exp(3,:)); hold on
		x = [2103.5 2103.5 2104 2104];
		y = [min(min_vaf(85,7),min_vaf(85,10))  max(max_vaf(85,7),max_vaf(85,10)) max(max_vaf(85,7),max_vaf(85,10)) min(min_vaf(85,7),min_vaf(85,10))]/-362.5;
		patch(x,y,[-1 -1 -1 -1],colors(7,:),'edgecolor',colors(7,:)); hold on
		y = mean_exp(85,4)/-362.5+[0.4 -0.4 -0.4 0.4];
		patch(x,y,[-1 -1 -1 -1],colors_exp(4,:),'edgecolor',colors_exp(4,:)); hold on
		x = [2104.5 2104.5 2105 2105];
		y = [min(min_vaf(85,8),min_vaf(85,11))  max(max_vaf(85,8),max_vaf(85,11)) max(max_vaf(85,8),max_vaf(85,11)) min(min_vaf(85,8),min_vaf(85,11))]/-362.5;
		patch(x,y,[-1 -1 -1 -1],colors(8,:),'edgecolor',colors(8,:)); hold on
		y = mean_exp(85,5)/-362.5+[0.4 -0.4 -0.4 0.4];
		patch(x,y,[-1 -1 -1 -1],colors_exp(5,:),'edgecolor',colors_exp(5,:)); hold on
		x = [2105.5 2105.5 2106 2106];
		y = [min(min_vaf(85,9),min_vaf(85,12))  max(max_vaf(85,9),max_vaf(85,12)) max(max_vaf(85,9),max_vaf(85,12)) min(min_vaf(85,9),min_vaf(85,12))]/-362.5;
		patch(x,y,[-1 -1 -1 -1],colors(9,:),'edgecolor',colors(9,:)); hold on
		y = mean_exp(85,6)/-362.5+[0.4 -0.4 -0.4 0.4];
		patch(x,y,[-1 -1 -1 -1],colors_exp(6,:),'edgecolor',colors_exp(6,:)); hold on
		legend({'NorESM1','MIROC','CCSM4','HadGEM2','CSIRO','IPSL'},'location','NorthWest')
		legend boxoff
		xlim([2015 2100])
		xlabel('Time (yr)');
		ylabel('Sea Level Contribution (mm SLE)');
		set(gcf,'Position',[400 500 800 450]);

		h = gcf;
		set(h,'Units','Inches');
		pos = get(h,'Position');
		set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
		pos = get(h,'Position');
		print(gcf, '-dpdf', '-painters', 'FiguresPaperFinal/Figure7.pdf');
	end

end %}}}

%Figure 8: all RCP8.5 by region 
if step==13 % {{{Figure 8

	experiments_list={'exp01','exp02','exp04','exp05','exp06','exp08','expA1','expA2','expA3','expA5','expA6','expA7'}; %open and standard RCP8.5

	results_west_1=[]; results_east_1=[]; results_penin_1=[]; 
	results_west_2=[]; results_east_2=[]; results_penin_2=[]; 
	results_west_3=[]; results_east_3=[]; results_penin_3=[]; 
	results_west_4=[]; results_east_4=[]; results_penin_4=[]; 
	results_west_5=[]; results_east_5=[]; results_penin_5=[]; 
	results_west_6=[]; results_east_6=[]; results_penin_6=[]; 

	if 1,
		for imodel=1:length(model_list),
			modelname=model_list{imodel};
			for iexp=1:length(experiments_list),
				expename=experiments_list{iexp};
				specifics;

				eval(['isexp=is' expename ';'])
				if isexp,
					explimnsw_file=['' scalarpath '/' group '/' simul '/' expename '/computed_ivaf_minus_ctrl_proj_AIS_' group '_' simul '_' expename '.nc'];
					time_model=ncread(explimnsw_file,'time');
					limnsw_west=ncread(explimnsw_file,'ivaf_region_1')*ice_density/(10^9*1000); %in Gt
					limnsw_east=ncread(explimnsw_file,'ivaf_region_2')*ice_density/(10^9*1000);
					limnsw_penin=ncread(explimnsw_file,'ivaf_region_3')*ice_density/(10^9*1000);
					if strcmpi(expename,'exp01') |strcmpi(expename,'exp05'),
						results_west_1(end+1)=limnsw_west(85);
						results_east_1(end+1)=limnsw_east(85);
						results_penin_1(end+1)=limnsw_penin(85);
					elseif strcmpi(expename,'exp02') |strcmpi(expename,'exp06'),
						results_west_2(end+1)=limnsw_west(85);
						results_east_2(end+1)=limnsw_east(85);
						results_penin_2(end+1)=limnsw_penin(85);
					elseif strcmpi(expename,'exp04') |strcmpi(expename,'exp08'),
						results_west_3(end+1)=limnsw_west(85);
						results_east_3(end+1)=limnsw_east(85);
						results_penin_3(end+1)=limnsw_penin(85);
					elseif strcmpi(expename,'expA1') |strcmpi(expename,'expA5'),
						results_west_4(end+1)=limnsw_west(85);
						results_east_4(end+1)=limnsw_east(85);
						results_penin_4(end+1)=limnsw_penin(85);
					elseif strcmpi(expename,'expA2') |strcmpi(expename,'expA6'),
						results_west_5(end+1)=limnsw_west(85);
						results_east_5(end+1)=limnsw_east(85);
						results_penin_5(end+1)=limnsw_penin(85);
					elseif strcmpi(expename,'expA3') |strcmpi(expename,'expA7'),
						results_west_6(end+1)=limnsw_west(85);
						results_east_6(end+1)=limnsw_east(85);
						results_penin_6(end+1)=limnsw_penin(85);
					else error('exp not supported');
					end
				end

			end %end of model
		end %end of isexp

		results=[mean(results_west_1), mean(results_west_2),mean(results_west_3) mean(results_west_4),mean(results_west_5), mean(results_west_6);...
			mean(results_east_1),mean(results_east_2),mean(results_east_3),mean(results_east_4),mean(results_east_5),mean(results_east_6);...
			mean(results_penin_1),mean(results_penin_2),mean(results_penin_3),mean(results_penin_4),mean(results_penin_5),mean(results_penin_6)];
		results_std=[std(results_west_1), std(results_west_2),std(results_west_3) std(results_west_4),std(results_west_5), std(results_west_6);...
			std(results_east_1),std(results_east_2),std(results_east_3),std(results_east_4),std(results_east_5),std(results_east_6);...
			std(results_penin_1),std(results_penin_2),std(results_penin_3),std(results_penin_4),std(results_penin_5),std(results_penin_6)];
		%Plot results
		close all; clf;
		colors = get(gca,'colororder');
		figure(1); set(gcf,'color','w'); set(gcf,'Position',[400 500 900 400]);
		hb=bar(results*(-1/362.5)); hold on 
		for ib = 1:numel(hb)
			%XData property is the tick labels/group centers; XOffset is the offset of each distinct group
			xData = hb(ib).XData+hb(ib).XOffset;
		   errorbar(xData(1),results(1,ib)*(-1/362.5),results_std(1,ib)*(-1/362.5),'color','k'); 
		   errorbar(xData(2),results(2,ib)*(-1/362.5),results_std(2,ib)*(-1/362.5),'color','k'); 
		   errorbar(xData(3),results(3,ib)*(-1/362.5),results_std(3,ib)*(-1/362.5),'color','k'); 
		end
		set(gca, 'XTickLabel', '');
		h=text(1,-78,'WAIS','VerticalAlignment','middle','HorizontalAlignment','center','color','k','fontweight','b');
		h=text(2,-78,'EAIS','VerticalAlignment','middle','HorizontalAlignment','center','color','k','fontweight','b');
		h=text(3,-78,'Peninsula','VerticalAlignment','middle','HorizontalAlignment','center','color','k','fontweight','b');
		for i=1:6,
			hb(i).FaceColor = colors(i,:);
		end
		ylabel('Sea Level Contribution (mm SLE)');
		xlim([0.5 3.5]);
		ylim([-70 150]);
		legend({'NorESM1','MIROC','CCSM','HadGEM2','CSIRO','IPSL'},'location','NorthEast')

		h = gcf;
		set(h,'Units','Inches');
		pos = get(h,'Position');
		set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
		print(gcf, '-dpdf', '-painters', 'FiguresPaperFinal/Figure8.pdf');
	end

end %}}}

%Figure 9: RCP8.5 vs RCP2.6
if step==14 % {{{Figure 9

	%Remove UTAS and expA7 in VUB to have models that did both 8.5 and 2.6
	model_list={'AWI_PISM1','DOE_MALI','ILTS_PIK_SICOPOLIS','IMAU_IMAUICE1','IMAU_IMAUICE2','JPL1_ISSM','LSCE_GRISLI','NCAR_CISM','PIK_PISM1','PIK_PISM2','UCIJPL_ISSM','ULB_FETISH32','ULB_FETISH16','VUB_AISMPALEO','VUW_PISM'};
	figure(1); set(gcf,'color','w'); set(gcf,'Position',[400 500 800 450]);

	experiments_list={'exp01','exp03','exp05','exp07','expA3','expA4','expA7','expA8'}; %open and standard NorESM and IPSL (8.5 vs 2.6)

	mean_vaf=zeros(86,length(experiments_list));
	std_total=zeros(86,length(experiments_list));
	num_models=zeros(length(experiments_list),1);
	time=[2015:2100];
	if 1,

		for iexp=1:length(experiments_list),
			expename=experiments_list{iexp};

			for imodel=1:length(model_list),
				modelname=model_list{imodel};
				specifics;

				eval(['isexp=is' expename ';'])
				if strcmp(expename,'expA7') & strcmp(modelname,'VUB_AISMPALEO'),
					%do nothing as RCP 2.6 is missing
				elseif isexp,
					explimnsw_file=['' scalarpath '/' group '/' simul '/' expename '/computed_ivaf_minus_ctrl_proj_AIS_' group '_' simul '_' expename '.nc'];
					time_model=ncread(explimnsw_file,'time');
					limnsw_model=ncread(explimnsw_file,'ivaf')*ice_density/(10^9*1000); %from m^3 to Gt
					mean_vaf(:,iexp)=mean_vaf(:,iexp)+[0;limnsw_model(1:85)];
					num_models(iexp)=num_models(iexp)+1;
				else
					%Experiment not done by this model, do nothing
				end

			end %end of model
		end %end of isexp
		mean_exp=zeros(86,4);
		mean_exp_extended=zeros(86,8);
		mean_exp(:,1)=(mean_vaf(:,1)+mean_vaf(:,3))/(num_models(1)+num_models(3)); 
		mean_exp(:,2)=(mean_vaf(:,2)+mean_vaf(:,4))/(num_models(2)+num_models(4)); 
		mean_exp(:,3)=(mean_vaf(:,5)+mean_vaf(:,7))/(num_models(5)+num_models(7)); 
		mean_exp(:,4)=(mean_vaf(:,6)+mean_vaf(:,8))/(num_models(6)+num_models(8)); 
		mean_exp_extended(:,1)=mean_exp(:,1);
		mean_exp_extended(:,2)=mean_exp(:,2);
		mean_exp_extended(:,3)=mean_exp(:,1);
		mean_exp_extended(:,4)=mean_exp(:,2);
		mean_exp_extended(:,5)=mean_exp(:,3);
		mean_exp_extended(:,6)=mean_exp(:,4);
		mean_exp_extended(:,7)=mean_exp(:,3);
		mean_exp_extended(:,8)=mean_exp(:,4);

		for iexp=1:length(experiments_list),
			expename=experiments_list{iexp};

			for imodel=1:length(model_list),
				modelname=model_list{imodel};
				specifics;

				eval(['isexp=is' expename ';'])
				if strcmp(expename,'expA7') & strcmp(modelname,'VUB_AISMPALEO'),
					%do nothing as RCP 2.6 is missing
				elseif isexp,
					explimnsw_file=['' scalarpath '/' group '/' simul '/' expename '/computed_ivaf_minus_ctrl_proj_AIS_' group '_' simul '_' expename '.nc'];
					time_model=ncread(explimnsw_file,'time');
					limnsw_model=ncread(explimnsw_file,'ivaf')*ice_density/(10^9*1000); %from m^3 to Gt
					std_total(:,iexp)=std_total(:,iexp)+([0;limnsw_model(1:85)]-mean_exp_extended(:,iexp)).^2;
				else
					%Experiment not done by this model, do nothing
				end

			end %end of model
		end %end of isexp

		std_exp=zeros(86,3);
		std_exp(:,1)=sqrt((std_total(:,1)+std_total(:,3))/(num_models(1)+num_models(3)));
		std_exp(:,2)=sqrt((std_total(:,2)+std_total(:,4))/(num_models(2)+num_models(4)));
		std_exp(:,3)=sqrt((std_total(:,5)+std_total(:,7))/(num_models(5)+num_models(7)));
		std_exp(:,4)=sqrt((std_total(:,6)+std_total(:,8))/(num_models(6)+num_models(8)));

		mean_exp=-mean_exp/362.5; %vaf in SLE mm
		std_exp=std_exp/362.5; %in SLE mm
		if 1, %Figure NorESM
			plot(time,mean_exp(:,1),'color','r','linewidth',2); hold on
			plot(time,mean_exp(:,2),'color','b','linewidth',2); hold on
			patch([time';flipud(time')],[mean_exp(:,1)-std_exp(:,1);flipud(mean_exp(:,1)+std_exp(:,1))],[1 0 0]+(1-[1 0 0])*0.25,'FaceAlpha',.3,'EdgeColor','None')
			patch([time';flipud(time')],[mean_exp(:,2)-std_exp(:,2);flipud(mean_exp(:,2)+std_exp(:,2))],[0 0 1]+(1-[0 0 1])*0.25,'FaceAlpha',.3,'EdgeColor','None')
			legend({'NorESM1 RCP 8.5 (exp01 & exp05)','NorESM1 RCP 2.6 (exp03 & exp07)'},'location','NorthWest')
			text(2010,-45,'a','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',18,'fontweight','b');
			legend boxoff
			xlim([2015 2100])
			xlabel('Time (yr)');
			ylabel('Sea Level Contribution (mm SLE)');
			set(gcf,'Position',[400 500 800 450]);

			h = gcf;
			set(h,'Units','Inches');
			pos = get(h,'Position');
			set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
			print(gcf, '-dpdf', '-painters', 'FiguresPaperFinal/Figure9a.pdf');
		end

		if 1, %Figure IPSL
			close
			plot(time,mean_exp(:,3),'color','r','linewidth',2); hold on
			plot(time,mean_exp(:,4),'color','b','linewidth',2); hold on
			patch([time';flipud(time')],[mean_exp(:,3)-std_exp(:,3);flipud(mean_exp(:,3)+std_exp(:,3))],[1 0 0]+(1-[1 0 0])*0.25,'FaceAlpha',.3,'EdgeColor','None')
			patch([time';flipud(time')],[mean_exp(:,4)-std_exp(:,4);flipud(mean_exp(:,4)+std_exp(:,4))],[0 0 1]+(1-[0 0 1])*0.25,'FaceAlpha',.3,'EdgeColor','None')
			legend({'IPSL RCP 8.5 (expA3 & expA7)','IPSL RCP 2.6 (expA4 & expA8)'},'location','NorthWest')
			text(2010,-55,'b','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',18,'fontweight','b');
			legend boxoff
			xlim([2015 2100])
			xlabel('Time (yr)');
			ylabel('Sea Level Contribution (mm SLE)');
			set(gcf,'Position',[400 500 800 450]);

			h = gcf;
			set(h,'Units','Inches');
			pos = get(h,'Position');
			set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
			print(gcf, '-dpdf', '-painters', 'FiguresPaperFinal/Figure9b.pdf');
		end

	end

end %}}}

%Figure 10: RCP8.5 vs RCP2.6 by region
if step==15 % {{{Figure 10a

	%Remove UTAS exp05 have models that did both 8.5 and 2.6
	model_list={'AWI_PISM1','DOE_MALI','ILTS_PIK_SICOPOLIS','IMAU_IMAUICE1','IMAU_IMAUICE2','JPL1_ISSM','LSCE_GRISLI','NCAR_CISM','PIK_PISM1','PIK_PISM2','UCIJPL_ISSM','ULB_FETISH32','ULB_FETISH16','VUB_AISMPALEO','VUW_PISM'};
	colors = distinguishable_colors(length(model_list));
	experiments_list={'exp05','exp01','exp07','exp03'}; %open and standard NorESM

	results_west_85=[]; results_east_85=[]; results_penin_85=[]; 
	resultstendacabfgr_west_85=[]; resultstendacabfgr_east_85=[]; resultstendacabfgr_penin_85=[]; 
	results_west_26=[]; results_east_26=[]; results_penin_26=[]; 
	resultstendacabfgr_west_26=[]; resultstendacabfgr_east_26=[]; resultstendacabfgr_penin_26=[]; 
	results_model={}; 

	if 1,
		for imodel=1:length(model_list),
			modelname=model_list{imodel};
			for iexp=1:length(experiments_list),
				expename=experiments_list{iexp};
				specifics;

				eval(['isexp=is' expename ';'])
				if isexp,
					explimnsw_file=['' scalarpath '/' group '/' simul '/' expename '/computed_ivaf_minus_ctrl_proj_AIS_' group '_' simul '_' expename '.nc'];
					exptendacabfgr_file=['' scalarpath '/' group '/' simul '/' expename '/computed_smbgr_minus_ctrl_proj_AIS_' group '_' simul '_' expename '.nc'];
					time_model=ncread(explimnsw_file,'time');
					limnsw_west=ncread(explimnsw_file,'ivaf_region_1')*ice_density/(10^9*1000); %from m^3 to Gt
					limnsw_east=ncread(explimnsw_file,'ivaf_region_2')*ice_density/(10^9*1000); %from m^3 to Gt
					limnsw_penin=ncread(explimnsw_file,'ivaf_region_3')*ice_density/(10^9*1000); %from m^3 to Gt
					tendacabfgr_west=ncread(exptendacabfgr_file,'smbgr_region_1')*yearday_model*3600*24/(10^9*1000); %from kg/s to Gt/yr
					tendacabfgr_east=ncread(exptendacabfgr_file,'smbgr_region_2')*yearday_model*3600*24/(10^9*1000); %from kg/s to Gt/yr
					tendacabfgr_penin=ncread(exptendacabfgr_file,'smbgr_region_3')*yearday_model*3600*24/(10^9*1000); %from kg/s to Gt/yr
					if strcmpi(expename,'exp05'),
						results_west_85(end+1)=limnsw_west(85);
						results_east_85(end+1)=limnsw_east(85);
						results_penin_85(end+1)=limnsw_penin(85);
						resultstendacabfgr_west_85(end+1)=sum(tendacabfgr_west(1:85));
						resultstendacabfgr_east_85(end+1)=sum(tendacabfgr_east(1:85));
						resultstendacabfgr_penin_85(end+1)=sum(tendacabfgr_penin(1:85));
						results_model{end+1}=[model_list2{imodel} '\_std'] ;
					elseif strcmpi(expename,'exp01'),
						results_west_85(end+1)=limnsw_west(85);
						results_east_85(end+1)=limnsw_east(85);
						results_penin_85(end+1)=limnsw_penin(85);
						resultstendacabfgr_west_85(end+1)=sum(tendacabfgr_west(1:85));
						resultstendacabfgr_east_85(end+1)=sum(tendacabfgr_east(1:85));
						resultstendacabfgr_penin_85(end+1)=sum(tendacabfgr_penin(1:85));
						results_model{end+1}=[model_list2{imodel} '\_open'] ;
					elseif strcmpi(expename,'exp07') | strcmpi(expename,'exp03'),
						results_west_26(end+1)=limnsw_west(85);
						results_east_26(end+1)=limnsw_east(85);
						results_penin_26(end+1)=limnsw_penin(85);
						resultstendacabfgr_west_26(end+1)=sum(tendacabfgr_west(1:85));
						resultstendacabfgr_east_26(end+1)=sum(tendacabfgr_east(1:85));
						resultstendacabfgr_penin_26(end+1)=sum(tendacabfgr_penin(1:85));
					end
				end
			end

		end %end of model
	end %end of isexp

	results_west=[results_west_85; results_west_26]; results_west=results_west(:)';
	results_east=[results_east_85; results_east_26]; results_east=results_east(:)';
	results_penin=[results_penin_85; results_penin_26]; results_penin=results_penin(:)';
	results=[results_west;results_east;results_penin];
	resultstendacabfgr_west=[resultstendacabfgr_west_85; resultstendacabfgr_west_26];
	resultstendacabfgr_west=resultstendacabfgr_west(:)';
	resultstendacabfgr_east=[resultstendacabfgr_east_85; resultstendacabfgr_east_26];
	resultstendacabfgr_east=resultstendacabfgr_east(:)';
	resultstendacabfgr_penin=[resultstendacabfgr_penin_85; resultstendacabfgr_penin_26];
	resultstendacabfgr_penin=resultstendacabfgr_penin(:)';
	resultstendacabfgr=[resultstendacabfgr_west;resultstendacabfgr_east;resultstendacabfgr_penin];

	%Plot results
	figure(1); set(gcf,'color','w'); set(gcf,'Position',[400 500 1200 550]);
	hb=bar(results*(-1/362.5),'grouped'); hold on 
	for iresults=1:numel(hb),
		if mod(iresults,2)==1,
			hb(iresults).FaceColor=[1 0 0]+(1-[1 0 0])*0.5;
		else
			hb(iresults).FaceColor=[0 0 1];
		end
	end
	set(gca, 'XTickLabel', '');
	h=text(1,170,'WAIS','VerticalAlignment','middle','HorizontalAlignment','center','color','k','fontweight','b');
	h=text(2,170,'EAIS','VerticalAlignment','middle','HorizontalAlignment','center','color','k','fontweight','b');
	h=text(3,170,'Peninsula','VerticalAlignment','middle','HorizontalAlignment','center','color','k','fontweight','b');
	ylabel('Sea Level Contribution (mm SLE)');
	text(0.4,-160,'a','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',18,'fontweight','b');
	xlim([0.5 3.5]);
	ylim([-50 160]);
	set(gca,'position',[0.10 0.35 0.835 0.585]);
	hold on

	%Add SMB results
	color_scenario=[0 0 1;[1 0 0]+(1-[1 0 0])*0.5];
	for ib = 1:numel(hb)
		%XData property is the tick labels/group centers; XOffset is the offset of each distinct group
		xData = hb(ib).XData+hb(ib).XOffset;
		plot(xData(1),resultstendacabfgr_west(ib)*(-1/362.5),'d','MarkerFaceColor',color_scenario(mod(ib,2)+1,:),'MarkerSize',6,'MarkerEdgeColor','k')
		plot(xData(2),resultstendacabfgr_east(ib)*(-1/362.5),'d','MarkerFaceColor',color_scenario(mod(ib,2)+1,:),'MarkerSize',6,'MarkerEdgeColor','k')
		plot(xData(3),resultstendacabfgr_penin(ib)*(-1/362.5),'d','MarkerFaceColor',color_scenario(mod(ib,2)+1,:),'MarkerSize',6,'MarkerEdgeColor','k')
	end
	for ib = 1:numel(hb)
		if mod(ib,2)==1,
			xData = (hb(ib).XData+hb(ib).XOffset + hb(ib+1).XData+hb(ib+1).XOffset)/2;
			h=text(xData(1),-50,results_model{(ib+1)/2},'VerticalAlignment','middle','HorizontalAlignment','right','color','k','fontweight','n');
			set(h,'Rotation',90);
			h=text(xData(2),-50,results_model{(ib+1)/2},'VerticalAlignment','middle','HorizontalAlignment','right','color','k','fontweight','n');
			set(h,'Rotation',90);
			h=text(xData(3),-50,results_model{(ib+1)/2},'VerticalAlignment','middle','HorizontalAlignment','right','color','k','fontweight','n');
			set(h,'Rotation',90);
		end
	end
	yticks([-50 0 50 100 150])
	legend({'RCP 8.5','RCP 2.6'},'location','NorthEast')
	h = gcf;
	set(h,'Units','Inches');
	pos = get(h,'Position');
	set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
	print(gcf, '-dpdf', '-painters', 'FiguresPaperFinal/Figure10a.pdf');

end %}}}
if step==16 % {{{Figure 10b

	%Remove expA7 in VUB to have models that did both 8.5 and 2.6
	colors = distinguishable_colors(length(model_list));
	experiments_list={'expA7','expA3','expA8','expA4'}; %open and standard NorESM

	results_west_85=[]; results_east_85=[]; results_penin_85=[]; 
	resultstendacabfgr_west_85=[]; resultstendacabfgr_east_85=[]; resultstendacabfgr_penin_85=[]; 
	results_west_26=[]; results_east_26=[]; results_penin_26=[]; 
	resultstendacabfgr_west_26=[]; resultstendacabfgr_east_26=[]; resultstendacabfgr_penin_26=[]; 
	results_model={}; results_meltparam={};

	if 1,
		for imodel=1:length(model_list),
			modelname=model_list{imodel};
			for iexp=1:length(experiments_list),
				expename=experiments_list{iexp};
				specifics;

				eval(['isexp=is' expename ';'])
				if strcmp(expename,'expA7') & strcmp(modelname,'VUB_AISMPALEO'),
					%do nothing as RCP 2.6 is missing
				elseif isexp,
					explimnsw_file=['' scalarpath '/' group '/' simul '/' expename '/computed_ivaf_minus_ctrl_proj_AIS_' group '_' simul '_' expename '.nc'];
					exptendacabfgr_file=['' scalarpath '/' group '/' simul '/' expename '/computed_smbgr_minus_ctrl_proj_AIS_' group '_' simul '_' expename '.nc'];
					time_model=ncread(explimnsw_file,'time');
					limnsw_west=ncread(explimnsw_file,'ivaf_region_1')*ice_density/(10^9*1000); %from m^3 to Gt
					limnsw_east=ncread(explimnsw_file,'ivaf_region_2')*ice_density/(10^9*1000); %from m^3 to Gt
					limnsw_penin=ncread(explimnsw_file,'ivaf_region_3')*ice_density/(10^9*1000); %from m^3 to Gt 
					tendacabfgr_west=ncread(exptendacabfgr_file,'smbgr_region_1')*yearday_model*3600*24/(10^9*1000); %from kg/s to Gt/yr
					tendacabfgr_east=ncread(exptendacabfgr_file,'smbgr_region_2')*yearday_model*3600*24/(10^9*1000); %from kg/s to Gt/yr
					tendacabfgr_penin=ncread(exptendacabfgr_file,'smbgr_region_3')*yearday_model*3600*24/(10^9*1000); %from kg/s to Gt/yr
					if strcmpi(expename,'expA7'),
						results_west_85(end+1)=limnsw_west(85);
						results_east_85(end+1)=limnsw_east(85);
						results_penin_85(end+1)=limnsw_penin(85);
						resultstendacabfgr_west_85(end+1)=sum(tendacabfgr_west(1:85));
						resultstendacabfgr_east_85(end+1)=sum(tendacabfgr_east(1:85));
						resultstendacabfgr_penin_85(end+1)=sum(tendacabfgr_penin(1:85));
						results_model{end+1}=[model_list2{imodel} '\_std'] ;
					elseif strcmpi(expename,'expA3'),
						results_west_85(end+1)=limnsw_west(85);
						results_east_85(end+1)=limnsw_east(85);
						results_penin_85(end+1)=limnsw_penin(85);
						resultstendacabfgr_west_85(end+1)=sum(tendacabfgr_west(1:85));
						resultstendacabfgr_east_85(end+1)=sum(tendacabfgr_east(1:85));
						resultstendacabfgr_penin_85(end+1)=sum(tendacabfgr_penin(1:85));
						results_model{end+1}=[model_list2{imodel} '\_open'] ;
					elseif strcmpi(expename,'expA8'),
						results_west_26(end+1)=limnsw_west(85);
						results_east_26(end+1)=limnsw_east(85);
						results_penin_26(end+1)=limnsw_penin(85);
						resultstendacabfgr_west_26(end+1)=sum(tendacabfgr_west(1:85));
						resultstendacabfgr_east_26(end+1)=sum(tendacabfgr_east(1:85));
						resultstendacabfgr_penin_26(end+1)=sum(tendacabfgr_penin(1:85));
					elseif strcmpi(expename,'expA4'),
						results_west_26(end+1)=limnsw_west(85);
						results_east_26(end+1)=limnsw_east(85);
						results_penin_26(end+1)=limnsw_penin(85);
						resultstendacabfgr_west_26(end+1)=sum(tendacabfgr_west(1:85));
						resultstendacabfgr_east_26(end+1)=sum(tendacabfgr_east(1:85));
						resultstendacabfgr_penin_26(end+1)=sum(tendacabfgr_penin(1:85));
					end
				end

			end %end of model
		end %end of isexp

		results_west=[results_west_85; results_west_26]; results_west=results_west(:)';
		results_east=[results_east_85; results_east_26]; results_east=results_east(:)';
		results_penin=[results_penin_85; results_penin_26]; results_penin=results_penin(:)';
		results=[results_west;results_east;results_penin];
		resultstendacabfgr_west=[resultstendacabfgr_west_85; resultstendacabfgr_west_26];
		resultstendacabfgr_west=resultstendacabfgr_west(:)';
		resultstendacabfgr_east=[resultstendacabfgr_east_85; resultstendacabfgr_east_26];
		resultstendacabfgr_east=resultstendacabfgr_east(:)';
		resultstendacabfgr_penin=[resultstendacabfgr_penin_85; resultstendacabfgr_penin_26];
		resultstendacabfgr_penin=resultstendacabfgr_penin(:)';
		resultstendacabfgr=[resultstendacabfgr_west;resultstendacabfgr_east;resultstendacabfgr_penin];

		%Plot results
		figure(1); set(gcf,'color','w'); set(gcf,'Position',[400 500 1200 550]);
		hb=bar(results*(-1/362.5),'grouped'); hold on 
		for iresults=1:numel(hb),
			if mod(iresults,2)==1,
				hb(iresults).FaceColor=[1 0 0]+(1-[1 0 0])*0.5;
			else
				hb(iresults).FaceColor=[0 0 1];
			end
		end
		set(gca, 'XTickLabel', '');
		h=text(1,170,'WAIS','VerticalAlignment','middle','HorizontalAlignment','center','color','k','fontweight','b');
		h=text(2,170,'EAIS','VerticalAlignment','middle','HorizontalAlignment','center','color','k','fontweight','b');
		h=text(3,170,'Peninsula','VerticalAlignment','middle','HorizontalAlignment','center','color','k','fontweight','b');
		ylabel('Sea Level Contribution (mm SLE)');
		xlim([0.5 3.5]);
		ylim([-50 160]);
		set(gca,'position',[0.10 0.32 0.835 0.555]);
		text(0.4,-100,'b','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',18,'fontweight','b');
		hold on

		%Add SMB results
		color_scenario=[0 0 1;[1 0 0]+(1-[1 0 0])*0.5];
		for ib = 1:numel(hb)
			%XData property is the tick labels/group centers; XOffset is the offset of each distinct group
			xData = hb(ib).XData+hb(ib).XOffset;
			plot(xData(1),resultstendacabfgr_west(ib)*(-1/362.5),'d','MarkerFaceColor',color_scenario(mod(ib,2)+1,:),'MarkerSize',6,'MarkerEdgeColor','k')
			plot(xData(2),resultstendacabfgr_east(ib)*(-1/362.5),'d','MarkerFaceColor',color_scenario(mod(ib,2)+1,:),'MarkerSize',6,'MarkerEdgeColor','k')
			plot(xData(3),resultstendacabfgr_penin(ib)*(-1/362.5),'d','MarkerFaceColor',color_scenario(mod(ib,2)+1,:),'MarkerSize',6,'MarkerEdgeColor','k')
		end
		for ib = 1:numel(hb)
			if mod(ib,2)==1,
				xData = (hb(ib).XData+hb(ib).XOffset + hb(ib+1).XData+hb(ib+1).XOffset)/2;
				h=text(xData(1),-50,results_model{(ib+1)/2},'VerticalAlignment','middle','HorizontalAlignment','right','color','k','fontweight','n');
				set(h,'Rotation',90);
				h=text(xData(2),-50,results_model{(ib+1)/2},'VerticalAlignment','middle','HorizontalAlignment','right','color','k','fontweight','n');
				set(h,'Rotation',90);
				h=text(xData(3),-50,results_model{(ib+1)/2},'VerticalAlignment','middle','HorizontalAlignment','right','color','k','fontweight','n');
				set(h,'Rotation',90);
			end
		end
		legend({'RCP 8.5','RCP 2.6'},'location','NorthEast')

		h = gcf;
		set(h,'Units','Inches');
		pos = get(h,'Position');
		set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
		print(gcf, '-dpdf', '-painters', 'FiguresPaperFinal/Figure10b.pdf');
	end

end %}}}

%Figure 11: Open vs standard melt parameterization
if step==17 % {{{Figure 11a

	colors = get(gca,'colororder');
	experiments_list={'exp01','exp05','exp02','exp06','exp04','exp08','expA1','expA5','expA2','expA6','expA3','expA7'}; %open and standard RCP8.5
 
	for iexp=1:length(experiments_list),
		eval(['results_west_' int2str(iexp) '=[];'])
		eval(['results_east_' int2str(iexp) '=[];'])
		eval(['results_penin_' int2str(iexp) '=[];'])
	end

	for imodel=1:length(model_list),
		modelname=model_list{imodel};
		for iexp=1:length(experiments_list),
			expename=experiments_list{iexp};
			specifics;

			eval(['isexp=is' expename ';'])
			if isexp,
				exptendlibmassbffl_file=['' scalarpath '/' group '/' simul '/' expename '/computed_bmbfl_minus_ctrl_proj_AIS_' group '_' simul '_' expename '.nc'];
				time_model=ncread(exptendlibmassbffl_file,'time');
				tendlibmassbffl_west=ncread(exptendlibmassbffl_file,'bmbfl_region_1')*yearday_model*3600*24/(10^9*1000); %from kg/s to Gt/yr
				tendlibmassbffl_east=ncread(exptendlibmassbffl_file,'bmbfl_region_2')*yearday_model*3600*24/(10^9*1000); %from kg/s to Gt/yr
				tendlibmassbffl_penin=ncread(exptendlibmassbffl_file,'bmbfl_region_3')*yearday_model*3600*24/(10^9*1000); %from kg/s to Gt/yr
				eval(['results_west_' int2str(iexp) '(end+1)=sum(tendlibmassbffl_west(1:85));'])
				eval(['results_east_' int2str(iexp) '(end+1)=sum(tendlibmassbffl_east(1:85));'])
				eval(['results_penin_' int2str(iexp) '(end+1)=sum(tendlibmassbffl_penin(1:85));'])
			end
		end %end of model
	end %end of isexp

	results=[mean(results_west_1), mean(results_west_2),mean(results_west_3), mean(results_west_4),mean(results_west_5), mean(results_west_6),...
		mean(results_west_7), mean(results_west_8),mean(results_west_9), mean(results_west_10),mean(results_west_11), mean(results_west_12);...
		mean(results_east_1),mean(results_east_2),mean(results_east_3),mean(results_east_4),mean(results_east_5),mean(results_east_6),...
		mean(results_east_7), mean(results_east_8),mean(results_east_9), mean(results_east_10),mean(results_east_11), mean(results_east_12);...
		mean(results_penin_1),mean(results_penin_2),mean(results_penin_3),mean(results_penin_4),mean(results_penin_5),mean(results_penin_6),...
		mean(results_penin_7), mean(results_penin_8),mean(results_penin_9), mean(results_penin_10),mean(results_penin_11), mean(results_penin_12)];
	results_std=[std(results_west_1), std(results_west_2),std(results_west_3), std(results_west_4),std(results_west_5), std(results_west_6),...
		std(results_west_7), std(results_west_8),std(results_west_9), std(results_west_10),std(results_west_11), std(results_west_12);...
		std(results_east_1),std(results_east_2),std(results_east_3),std(results_east_4),std(results_east_5),std(results_east_6),...
		std(results_east_7), std(results_east_8),std(results_east_9), std(results_east_10),std(results_east_11), std(results_east_12);...
		std(results_penin_1),std(results_penin_2),std(results_penin_3),std(results_penin_4),std(results_penin_5),std(results_penin_6),...
		std(results_penin_7), std(results_penin_8),std(results_penin_9), std(results_penin_10),std(results_penin_11), std(results_penin_12)];
	%Plot results
	figure(1); set(gcf,'color','w'); set(gcf,'Position',[400 500 800 450]);
	hb=bar(-results); hold on 

	for ib = 1:numel(hb)
		%XData property is the tick labels/group centers; XOffset is the offset of each distinct group
		xData = hb(ib).XData+hb(ib).XOffset;
		errorbar(xData(1),-results(1,ib),results_std(1,ib)*(-1),'color','k'); 
		errorbar(xData(2),-results(2,ib),results_std(2,ib)*(-1),'color','k'); 
		errorbar(xData(3),-results(3,ib),results_std(3,ib)*(-1),'color','k'); 
	end
	set(gca, 'XTickLabel', '');
	h=text(1,-0.63*10^5,'WAIS','VerticalAlignment','middle','HorizontalAlignment','center','color','k','fontweight','b');
	h=text(2,-0.63*10^5,'EAIS','VerticalAlignment','middle','HorizontalAlignment','center','color','k','fontweight','b');
	h=text(3,-0.63*10^5,'Peninsula','VerticalAlignment','middle','HorizontalAlignment','center','color','k','fontweight','b');
	for i=1:6,
		hb(2*i-1).FaceColor = colors(i,:);
		hb(2*i).FaceColor = colors(i,:);
	end
	ylabel('Total ocean melt (Gt)');
	xlim([0.5 3.5]);
	ylim([-5 30]*10^4);
	text(0.35,-0.65*10^5,'a','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',18,'fontweight','b');
	legend({'NorESM open','NorESM standard','MIROC open','MIROC standard','CCSM open','CCSM standard',...
		'HadGEM2 open','HadGEM2 standard','CSIRO open','CSIRO standard','IPSL open','IPSL standard'},'location','EastOutside')

	h = gcf;
	set(h,'Units','Inches');
	pos = get(h,'Position');
	set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
	print(gcf, '-dpdf', '-painters', 'FiguresPaperFinal/Figure11a.pdf');

end %}}}
if step==18 % {{{Figure 11b

	colors = get(gca,'colororder');
	experiments_list={'exp01','exp05','exp02','exp06','exp04','exp08','expA1','expA5','expA2','expA6','expA3','expA7'}; %open and standard RCP8.5
 
	for iexp=1:length(experiments_list),
		eval(['results_west_' int2str(iexp) '=[];'])
		eval(['results_east_' int2str(iexp) '=[];'])
		eval(['results_penin_' int2str(iexp) '=[];'])
	end

	for imodel=1:length(model_list),
		modelname=model_list{imodel};
		for iexp=1:length(experiments_list),
			expename=experiments_list{iexp};
			specifics;

			eval(['isexp=is' expename ';'])
			if isexp,
				explimnsw_file=['' scalarpath '/' group '/' simul '/' expename '/computed_ivaf_minus_ctrl_proj_AIS_' group '_' simul '_' expename '.nc'];
				time_model=ncread(explimnsw_file,'time');
				limnsw_west=ncread(explimnsw_file,'ivaf_region_1')*ice_density/(10^9*1000); %from m^3 to Gt
				limnsw_east=ncread(explimnsw_file,'ivaf_region_2')*ice_density/(10^9*1000); %from m^3 to Gt
				limnsw_penin=ncread(explimnsw_file,'ivaf_region_3')*ice_density/(10^9*1000); %from m^3 to Gt
				eval(['results_west_' int2str(iexp) '(end+1)=limnsw_west(85);'])
				eval(['results_east_' int2str(iexp) '(end+1)=limnsw_east(85);'])
				eval(['results_penin_' int2str(iexp) '(end+1)=limnsw_penin(85);'])
			end
		end %end of model
	end %end of isexp

	results=[mean(results_west_1), mean(results_west_2),mean(results_west_3), mean(results_west_4),mean(results_west_5), mean(results_west_6),...
		mean(results_west_7), mean(results_west_8),mean(results_west_9), mean(results_west_10),mean(results_west_11), mean(results_west_12);...
		mean(results_east_1),mean(results_east_2),mean(results_east_3),mean(results_east_4),mean(results_east_5),mean(results_east_6),...
		mean(results_east_7), mean(results_east_8),mean(results_east_9), mean(results_east_10),mean(results_east_11), mean(results_east_12);...
		mean(results_penin_1),mean(results_penin_2),mean(results_penin_3),mean(results_penin_4),mean(results_penin_5),mean(results_penin_6),...
		mean(results_penin_7), mean(results_penin_8),mean(results_penin_9), mean(results_penin_10),mean(results_penin_11), mean(results_penin_12)];
	results_std=[std(results_west_1), std(results_west_2),std(results_west_3), std(results_west_4),std(results_west_5), std(results_west_6),...
		std(results_west_7), std(results_west_8),std(results_west_9), std(results_west_10),std(results_west_11), std(results_west_12);...
		std(results_east_1),std(results_east_2),std(results_east_3),std(results_east_4),std(results_east_5),std(results_east_6),...
		std(results_east_7), std(results_east_8),std(results_east_9), std(results_east_10),std(results_east_11), std(results_east_12);...
		std(results_penin_1),std(results_penin_2),std(results_penin_3),std(results_penin_4),std(results_penin_5),std(results_penin_6),...
		std(results_penin_7), std(results_penin_8),std(results_penin_9), std(results_penin_10),std(results_penin_11), std(results_penin_12)];
	%Plot results
	figure(1); set(gcf,'color','w'); set(gcf,'Position',[400 500 800 450]);
	hb=bar(results*(-1/361.8)); hold on 

	for ib = 1:numel(hb)
		%XData property is the tick labels/group centers; XOffset is the offset of each distinct group
		xData = hb(ib).XData+hb(ib).XOffset;
		errorbar(xData(1),results(1,ib)*(-1/361.8),results_std(1,ib)*(-1/361.8),'color','k'); 
		errorbar(xData(2),results(2,ib)*(-1/361.8),results_std(2,ib)*(-1/361.8),'color','k'); 
		errorbar(xData(3),results(3,ib)*(-1/361.8),results_std(3,ib)*(-1/361.8),'color','k'); 
	end
	set(gca, 'XTickLabel', '');
	h=text(1,-88,'WAIS','VerticalAlignment','middle','HorizontalAlignment','center','color','k','fontweight','b');
	h=text(2,-88,'EAIS','VerticalAlignment','middle','HorizontalAlignment','center','color','k','fontweight','b');
	h=text(3,-88,'Peninsula','VerticalAlignment','middle','HorizontalAlignment','center','color','k','fontweight','b');
	for i=1:6,
		hb(2*i-1).FaceColor = colors(i,:);
		hb(2*i).FaceColor = colors(i,:);
	end
	ylabel('Additional SLE (mm)');
	xlim([0.5 3.5]);
	ylim([-80 200]);
	text(0.35,-85,'b','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',18,'fontweight','b');
	legend({'NorESM open','NorESM standard','MIROC open','MIROC standard','CCSM4 open','CCSM4 standard',...
		'HadGEM2 open','HadGEM2 standard','CSIRO open','CSIRO standard','IPSL open','IPSL standard'},'location','EastOutside')

	h = gcf;
	set(h,'Units','Inches');
	pos = get(h,'Position');
	set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
	print(gcf, '-dpdf', '-painters', 'FiguresPaperFinal/Figure11b.pdf');
end %}}}

%Figure 12: Melt parameterization
if step==19 % {{{Figure 12a

	%Remove UTAS to include only models that did both the 3 model uncertainty NorESM
	model_list={'AWI_PISM1','DOE_MALI','ILTS_PIK_SICOPOLIS','IMAU_IMAUICE1','IMAU_IMAUICE2','JPL1_ISSM','LSCE_GRISLI','NCAR_CISM','PIK_PISM1','PIK_PISM2','UCIJPL_ISSM','ULB_FETISH32','ULB_FETISH16','VUB_AISMPALEO','VUW_PISM'};
	experiments_list={'exp05','exp09','exp10'}; %standard NorESM 

	mean_basalmelt=zeros(86,length(experiments_list));
	min_basalmelt=9999*ones(86,length(experiments_list));
	max_basalmelt=-9999*ones(86,length(experiments_list));
	num_models=zeros(length(experiments_list),1);
	time=[2015:2100];

	for iexp=1:length(experiments_list),
		expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			modelname=model_list{imodel};
			specifics;

			eval(['isexp=is' expename ';'])
			if isexp,
				num_models(iexp)=num_models(iexp)+1;
				exptendlibmassbffl_file=['' scalarpath '/' group '/' simul '/' expename '/computed_bmbfl_minus_ctrl_proj_AIS_' group '_' simul '_' expename '.nc'];
				time_model=ncread(exptendlibmassbffl_file,'time');
				tendlibmassbffl_model=ncread(exptendlibmassbffl_file,'bmbfl')*yearday_model*3600*24/(10^9*1000); %from kg/s to Gt/yr
				mean_basalmelt(:,iexp)=mean_basalmelt(:,iexp)+[0;tendlibmassbffl_model(1:85)];
				min_basalmelt(:,iexp)=min(min_basalmelt(:,iexp),[0;tendlibmassbffl_model(1:85)]);
				max_basalmelt(:,iexp)=max(max_basalmelt(:,iexp),[0;tendlibmassbffl_model(1:85)]);
			else
				%Experiment not done by this model, do nothing
			end

		end %end of model
	end %end of isexp

	mean_basalmelt=-(mean_basalmelt./num_models'); %positive ocean melt
	
	%Figure NorESM
	figure(1); set(gcf,'color','w'); set(gcf,'Position',[400 500 800 450]);
	plot(time,mean_basalmelt(:,1),'color','r','linewidth',2); hold on
	plot(time,mean_basalmelt(:,2),'color',[0.5172 0.5172 1],'linewidth',2); hold on
	plot(time,mean_basalmelt(:,3),'color',[0 0.5172 0.5862],'linewidth',2); hold on
	patch([time';flipud(time')],[-min_basalmelt(:,1);flipud(-max_basalmelt(:,1))],[1 0 0]+(1-[1 0 0])*0.25,'FaceAlpha',.3,'EdgeColor','None')
	patch([time';flipud(time')],[-min_basalmelt(:,2);flipud(-max_basalmelt(:,2))],[0.5172 0.5172 1]+(1-[0.5172 0.5172 1])*0.25,'FaceAlpha',.3,'EdgeColor','None')
	patch([time';flipud(time')],[-min_basalmelt(:,3);flipud(-max_basalmelt(:,3))],[0 0.5172 0.5862]+(1-[0 0.5172 0.5862])*0.25,'FaceAlpha',.3,'EdgeColor','None')
	legend({'Median Melt (exp05)','95% Melt (exp09)','5% Melt (exp10)'},'location','NorthWest')
	legend boxoff
	xlim([2015 2100])
	xlabel('Time (yr)');
	ylabel('Ocean Basal Melt (Gt/yr)');
	set(gcf,'Position',[400 500 800 450]);
	h = gcf;
	set(h,'Units','Inches');
	pos = get(h,'Position');
	set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
	text(2010,-1600,'a','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',18,'fontweight','b');
	print(gcf, '-dpdf', '-painters', 'FiguresPaperFinal/Figure12a.pdf');

end %}}}
if step==20 % {{{Figure 12b

	%All standard models did these runs, so no need to remove anyone

	experiments_list={'exp05','exp13'}; %standard NorESM 

	mean_basalmelt=zeros(86,length(experiments_list));
	min_basalmelt=9999*ones(86,length(experiments_list));
	max_basalmelt=-9999*ones(86,length(experiments_list));
	num_models=zeros(length(experiments_list),1);
	time=[2015:2100];

	for iexp=1:length(experiments_list),
		expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			modelname=model_list{imodel};
			specifics;

			eval(['isexp=is' expename ';'])
			if isexp,
				num_models(iexp)=num_models(iexp)+1;
				exptendlibmassbffl_file=['' scalarpath '/' group '/' simul '/' expename '/computed_bmbfl_minus_ctrl_proj_AIS_' group '_' simul '_' expename '.nc'];
				time_model=ncread(exptendlibmassbffl_file,'time');
				tendlibmassbffl_model=ncread(exptendlibmassbffl_file,'bmbfl')*yearday_model*3600*24/(10^9*1000); %from kg/s to Gt/yr
				mean_basalmelt(:,iexp)=mean_basalmelt(:,iexp)+[0;tendlibmassbffl_model(1:85)];
				min_basalmelt(:,iexp)=min(min_basalmelt(:,iexp),[0;tendlibmassbffl_model(1:85)]);
				max_basalmelt(:,iexp)=max(max_basalmelt(:,iexp),[0;tendlibmassbffl_model(1:85)]);
			else
				%Experiment not done by this model, do nothing
			end

		end %end of model
	end %end of isexp

	mean_basalmelt=-(mean_basalmelt./num_models'); %positive ocean melt

	%Figure NorESM
	figure(1); set(gcf,'color','w'); set(gcf,'Position',[400 500 800 450]);
	plot(time,mean_basalmelt(:,1),'color','r','linewidth',2); hold on
	plot(time,mean_basalmelt(:,2),'color',[0.8276 0.069 1],'linewidth',2); hold on
	patch([time';flipud(time')],[-min_basalmelt(:,1);flipud(-max_basalmelt(:,1))],[1 0 0]+(1-[1 0 0])*0.25,'FaceAlpha',.3,'EdgeColor','None')
	patch([time';flipud(time')],[-min_basalmelt(:,2);flipud(-max_basalmelt(:,2))],[0.8276 0.069 1]+(1-[0.8276 0.069 1])*0.25,'FaceAlpha',.3,'EdgeColor','None')
	legend({'MeanAnt (exp05)','PIGL (exp13)'},'location','NorthWest')
	text(2010,-7000,'b','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',18,'fontweight','b');
	legend boxoff
	xlim([2015 2100])
	xlabel('Time (yr)');
	ylabel('Ocean Basal Melt (Gt/yr)');
	set(gcf,'Position',[400 500 800 450]);
	h = gcf;
	set(h,'Units','Inches');
	pos = get(h,'Position');
	set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
	print(gcf, '-dpdf', '-painters', 'FiguresPaperFinal/Figure12b.pdf');

end %}}}
if step==21 % {{{Figure 12c

	%Remove UTAS to have models that did both the 3 model uncertainty NorESM
	model_list={'AWI_PISM1','DOE_MALI','ILTS_PIK_SICOPOLIS','IMAU_IMAUICE1','IMAU_IMAUICE2','JPL1_ISSM','LSCE_GRISLI','NCAR_CISM','PIK_PISM1','PIK_PISM2','UCIJPL_ISSM','ULB_FETISH32','ULB_FETISH16','VUB_AISMPALEO','VUW_PISM'};
	experiments_list={'exp05','exp09','exp10'}; %standard NorESM 

	mean_vaf=zeros(86,length(experiments_list));
	min_vaf=9999*ones(86,length(experiments_list));
	max_vaf=-9999*ones(86,length(experiments_list));
	num_models=zeros(length(experiments_list),1);
	time=[2015:2100];

	for iexp=1:length(experiments_list),
		expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			modelname=model_list{imodel};
			specifics;

			eval(['isexp=is' expename ';'])
			if isexp,
				explimnsw_file=['' scalarpath '/' group '/' simul '/' expename '/computed_ivaf_minus_ctrl_proj_AIS_' group '_' simul '_' expename '.nc'];
				time_model=ncread(explimnsw_file,'time');
				limnsw_model=ncread(explimnsw_file,'ivaf')*ice_density/(10^9*1000); %from m^3 to Gt
				mean_vaf(:,iexp)=mean_vaf(:,iexp)+[0;limnsw_model(1:85)];
				min_vaf(:,iexp)=min(min_vaf(:,iexp),[0;limnsw_model(1:85)]);
				max_vaf(:,iexp)=max(max_vaf(:,iexp),[0;limnsw_model(1:85)]);
				num_models(iexp)=num_models(iexp)+1;
			else
				%Experiment not done by this model, do nothing
			end

		end %end of model
	end %end of isexp

	mean_vaf=-(mean_vaf./num_models')/362.5; %vaf in SLE mm

	%Figure NorESM
	figure(1); set(gcf,'color','w'); set(gcf,'Position',[400 500 800 450]);
	plot(time,mean_vaf(:,1),'color','r','linewidth',2); hold on
	plot(time,mean_vaf(:,2),'color',[0.5172 0.5172 1],'linewidth',2); hold on
	plot(time,mean_vaf(:,3),'color',[0 0.5172 0.5862],'linewidth',2); hold on
	patch([time';flipud(time')],[min_vaf(:,1)/-362.5;flipud(max_vaf(:,1)/-362.5)],[1 0 0]+(1-[1 0 0])*0.25,'FaceAlpha',.3,'EdgeColor','None')
	patch([time';flipud(time')],[min_vaf(:,2)/-362.5;flipud(max_vaf(:,2)/-362.5)],[0.5172 0.5172 1]+(1-[0.5172 0.5172 1])*0.25,'FaceAlpha',.3,'EdgeColor','None')
	patch([time';flipud(time')],[min_vaf(:,3)/-362.5;flipud(max_vaf(:,3)/-362.5)],[0 0.5172 0.5862]+(1-[0 0.5172 0.5862])*0.25,'FaceAlpha',.3,'EdgeColor','None')
	legend({'Median Melt (exp05)','95% Melt (exp09)','5% Melt (exp10)'},'location','NorthWest')
	text(2010,-37,'c','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',18,'fontweight','b');
	legend boxoff
	xlim([2015 2100])
	xlabel('Time (yr)');
	ylabel('Sea Level Contribution (mm SLE)');
	set(gcf,'Position',[400 500 800 450]);

	h = gcf;
	set(h,'Units','Inches');
	pos = get(h,'Position');
	set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
	print(gcf, '-dpdf', '-painters', 'FiguresPaperFinal/Figure12c.pdf');

end %}}}
if step==22 % {{{Figure 12d

	%Evebody did it, no need to remove anyone
	experiments_list={'exp05','exp13'}; %standard NorESM AllMelt and PIGL calibration

	mean_vaf=zeros(86,length(experiments_list));
	min_vaf=9999*ones(86,length(experiments_list));
	max_vaf=-9999*ones(86,length(experiments_list));
	num_models=zeros(length(experiments_list),1);
	time=[2015:2100];

	for iexp=1:length(experiments_list),
		expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			modelname=model_list{imodel};
			specifics;

			eval(['isexp=is' expename ';'])
			if isexp,
				explimnsw_file=['' scalarpath '/' group '/' simul '/' expename '/computed_ivaf_minus_ctrl_proj_AIS_' group '_' simul '_' expename '.nc'];
				limnsw_model=ncread(explimnsw_file,'ivaf')*ice_density/(10^9*1000); %from m^3 to Gt
				mean_vaf(:,iexp)=mean_vaf(:,iexp)+[0;limnsw_model(1:85)];
				min_vaf(:,iexp)=min(min_vaf(:,iexp),[0;limnsw_model(1:85)]);
				max_vaf(:,iexp)=max(max_vaf(:,iexp),[0;limnsw_model(1:85)]);
				num_models(iexp)=num_models(iexp)+1;
			else
				%Experiment not done by this model, do nothing
			end

		end %end of model
	end %end of isexp

	mean_vaf=-(mean_vaf./num_models')/362.5; %vaf in SLE mm

	%Figure NorESM
	figure(1); set(gcf,'color','w'); set(gcf,'Position',[400 500 800 450]);
	plot(time,mean_vaf(:,1),'color','r','linewidth',2); hold on
	plot(time,mean_vaf(:,2),'color',[0.8276 0.069 1],'linewidth',2); hold on
	patch([time';flipud(time')],[min_vaf(:,1)/-362.5;flipud(max_vaf(:,1)/-362.5)],[1 0 0]+(1-[1 0 0])*0.25,'FaceAlpha',.3,'EdgeColor','None')
	patch([time';flipud(time')],[min_vaf(:,2)/-362.5;flipud(max_vaf(:,2)/-362.5)],[0.8276 0.069 1]+(1-[0.8276 0.069 1])*0.25,'FaceAlpha',.3,'EdgeColor','None')
	legend({'MeanAnt (exp05)','PIGL (exp13)'},'location','NorthWest')
	text(2010,-75,'d','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',18,'fontweight','b');
	legend boxoff
	xlim([2015 2100])
	xlabel('Time (yr)');
	ylabel('Sea Level Contribution (mm SLE)');
	set(gcf,'Position',[400 500 800 450]);

	h = gcf;
	set(h,'Units','Inches');
	pos = get(h,'Position');
	set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
	print(gcf, '-dpdf', '-painters', 'FiguresPaperFinal/Figure12d.pdf');

end %}}}

%Figure 13: Ice shelf collapse
if step==23 % {{{Figure 13a

	%Only models that did ice shelf collapse experiment
	model_list={'AWI_PISM1','DOE_MALI','ILTS_PIK_SICOPOLIS','IMAU_IMAUICE1','IMAU_IMAUICE2','JPL1_ISSM','LSCE_GRISLI','UCIJPL_ISSM','ULB_FETISH32','ULB_FETISH16'};
	experiments_list={'exp04','exp11','exp08','exp12'}; %open and standard NorESM with and without shelf collapse

	mean_floatingarea=zeros(86,length(experiments_list));
	std_total=zeros(86,length(experiments_list));
	num_models=zeros(length(experiments_list),1);
	time=[2015:2100];

	for iexp=1:length(experiments_list),
		expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			modelname=model_list{imodel};
			specifics;

			eval(['isexp=is' expename ';'])
			if isexp,
				num_models(iexp)=num_models(iexp)+1;
				expiareafl_file=['' scalarpath '/' group '/' simul '/' expename '/computed_iareafl_minus_ctrl_proj_AIS_' group '_' simul '_' expename '.nc'];
				time_model=ncread(expiareafl_file,'time');
				iareafl_model=ncread(expiareafl_file,'iareafl')/1000^2; % m^2 to km^2
				mean_floatingarea(:,iexp)=mean_floatingarea(:,iexp)+[0;iareafl_model(1:85)];
			else
				%Experiment not done by this model, do nothing
			end

		end %end of model
	end %end of isexp
	mean_exp=zeros(86,2);
	mean_exp(:,1)=(mean_floatingarea(:,1)+mean_floatingarea(:,3))/(num_models(1)+num_models(3)); 
	mean_exp(:,2)=(mean_floatingarea(:,2)+mean_floatingarea(:,4))/(num_models(2)+num_models(4)); 
	mean_exp_extended=repmat(mean_exp,1,2); %copy twice experiments 4 and 8 and exp 11 and 12

	for iexp=1:length(experiments_list),
		expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			modelname=model_list{imodel};
			specifics;

			eval(['isexp=is' expename ';'])
			if isexp,
				expiareafl_file=['' scalarpath '/' group '/' simul '/' expename '/computed_iareafl_minus_ctrl_proj_AIS_' group '_' simul '_' expename '.nc'];
				time_model=ncread(expiareafl_file,'time');
				iareafl_model=ncread(expiareafl_file,'iareafl')/1000^2; % m^2 to km^2
				std_total(:,iexp)=std_total(:,iexp)+([0;iareafl_model(1:85)]-mean_exp_extended(:,iexp)).^2;
			else
				%Experiment not done by this model, do nothing
			end

		end %end of model
	end %end of isexp

	std_exp=zeros(86,4);
	std_exp(:,1)=sqrt((std_total(:,1)+std_total(:,3))/(num_models(1)+num_models(3)));
	std_exp(:,2)=sqrt((std_total(:,2)+std_total(:,4))/(num_models(2)+num_models(4)));

	%Figure NorESM
	figure(1); set(gcf,'color','w'); set(gcf,'Position',[400 500 800 450]);
	plot(time,mean_exp(:,1),'color','r','linewidth',2); hold on
	plot(time,mean_exp(:,2),'color','c','linewidth',2); hold on
	patch([time';flipud(time')],[mean_exp(:,1)-std_exp(:,1);flipud(mean_exp(:,1)+std_exp(:,1))],[1 0 0]+(1-[1 0 0])*0.25,'FaceAlpha',.3,'EdgeColor','None')
	patch([time';flipud(time')],[mean_exp(:,2)-std_exp(:,2);flipud(mean_exp(:,2)+std_exp(:,2))],[0 1 1]+(1-[0 1 1])*0.25,'FaceAlpha',.3,'EdgeColor','None')
	text(2010,-21*10^4,'a','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',18,'fontweight','b');
	legend({'No shelf collapse (exp04 & exp08)','Ice shelf collapse (exp11 & exp12)'},'location','NorthWest')
	legend boxoff
	xlim([2015 2100])
	xlabel('Time (yr)');
	ylabel('Floating ice area (km^2)');

	set(gcf,'Position',[400 500 800 450]);
	h = gcf;
	set(h,'Units','Inches');
	pos = get(h,'Position');
	set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
	print(gcf, '-dpdf', '-painters', 'FiguresPaperFinal/Figure13a.pdf');

end %}}}
if step==24 % {{{Figure 13b

	%Only models that did ice shelf collapse experiment
	model_list={'AWI_PISM1','DOE_MALI','ILTS_PIK_SICOPOLIS','IMAU_IMAUICE1','IMAU_IMAUICE2','JPL1_ISSM','LSCE_GRISLI','UCIJPL_ISSM','ULB_FETISH32','ULB_FETISH16'};

	experiments_list={'exp04','exp11','exp08','exp12'}; %open and standard NorESM with and without shelf collapse

	mean_vaf=zeros(86,length(experiments_list));
	std_total=zeros(86,length(experiments_list));
	num_models=zeros(length(experiments_list),1);
	time=[2015:2100];

	for iexp=1:length(experiments_list),
		expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			modelname=model_list{imodel};
			specifics;

			eval(['isexp=is' expename ';'])
			if isexp,
				num_models(iexp)=num_models(iexp)+1;
				explimnsw_file=['' scalarpath '/' group '/' simul '/' expename '/computed_ivaf_minus_ctrl_proj_AIS_' group '_' simul '_' expename '.nc'];
				limnsw_model=ncread(explimnsw_file,'ivaf')*ice_density/(10^9*1000); %from m^3 to Gt
				time_model=ncread(explimnsw_file,'time');
				mean_vaf(:,iexp)=mean_vaf(:,iexp)+[0;limnsw_model(1:85)];
			else
				%Experiment not done by this model, do nothing
			end

		end %end of model
	end %end of isexp
	mean_exp=zeros(86,2);
	mean_exp(:,1)=(mean_vaf(:,1)+mean_vaf(:,3))/(num_models(1)+num_models(3)); 
	mean_exp(:,2)=(mean_vaf(:,2)+mean_vaf(:,4))/(num_models(2)+num_models(4)); 
	mean_exp_extended=repmat(mean_exp,1,2);

	for iexp=1:length(experiments_list),
		expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			modelname=model_list{imodel};
			specifics;

			eval(['isexp=is' expename ';'])
			if isexp,
				explimnsw_file=['' scalarpath '/' group '/' simul '/' expename '/computed_ivaf_minus_ctrl_proj_AIS_' group '_' simul '_' expename '.nc'];
				time_model=ncread(explimnsw_file,'time');
				limnsw_model=ncread(explimnsw_file,'ivaf')*ice_density/(10^9*1000); %from m^3 to Gt
				std_total(:,iexp)=std_total(:,iexp)+([0;limnsw_model(1:85)]-mean_exp_extended(:,iexp)).^2;
			else
				%Experiment not done by this model, do nothing
			end

		end %end of model
	end %end of isexp

	std_exp=zeros(86,3);
	std_exp(:,1)=sqrt((std_total(:,1)+std_total(:,3))/(num_models(1)+num_models(3)));
	std_exp(:,2)=sqrt((std_total(:,2)+std_total(:,4))/(num_models(2)+num_models(4)));

	mean_exp=-mean_exp/362.5; %vaf in SLE mm
	std_exp=std_exp/362.5; %in SLE mm

	%Figure NorESM
	figure(1); set(gcf,'color','w'); set(gcf,'Position',[400 500 800 450]);
	plot(time,mean_exp(:,1),'color','r','linewidth',2); hold on
	plot(time,mean_exp(:,2),'color','c','linewidth',2); hold on
	patch([time';flipud(time')],[mean_exp(:,1)-std_exp(:,1);flipud(mean_exp(:,1)+std_exp(:,1))],[1 0 0]+(1-[1 0 0])*0.25,'FaceAlpha',.3,'EdgeColor','None')
	patch([time';flipud(time')],[mean_exp(:,2)-std_exp(:,2);flipud(mean_exp(:,2)+std_exp(:,2))],[0 1 1]+(1-[0 1 1])*0.25,'FaceAlpha',.3,'EdgeColor','None')
	text(2010,-85,'b','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',18,'fontweight','b');
	legend({'No shelf collapse (exp04 & exp08)','Ice shelf collapse (exp11 & exp12)'},'location','NorthWest')
	legend boxoff
	xlim([2015 2100])
	xlabel('Time (yr)');
	ylabel('Sea Level Contribution (mm SLE)');
	set(gcf,'Position',[400 500 800 450]);
	h = gcf;
	set(h,'Units','Inches');
	pos = get(h,'Position');
	set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
	print(gcf, '-dpdf', '-painters', 'FiguresPaperFinal/Figure13b.pdf');

end %}}}

%Figure 14: Ice shelf collapse spatial changes
if step==25 % {{{Figure 14a

	experiments_list={'exp11','exp12'}; %ice shelf collapse experiments
	nummodels=zeros(761,761);
	totalthicknesschange=zeros(761,761);
	thicknessstd=zeros(761,761);

	for iexp=1:length(experiments_list),
		expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			modelname=model_list{imodel};
			specifics;
			eval(['exp_name=' expename '_regrid;'])
			exp_directory=['' gridpath '/' group '/' simul '/' exp_name '/'];
			if strcmp(expename,'exp11'),
				eval(['ctrl_proj_name=exp04_regrid;'])
				ctrl_directory=['' gridpath '/' group '/' simul '/' ctrl_proj_name '/'];
			elseif strcmp(expename,'exp12'),
				eval(['ctrl_proj_name=exp08_regrid;'])
				ctrl_directory=['' gridpath '/' group '/' simul '/' ctrl_proj_name '/'];
			end

			eval(['isexp=is' expename ';'])
			if isexp,
				if exist(exp_directory,'dir')==0,
					error(['directory ' exp_directory ' not found']);
				end

				if is_lithk==1,
					field='lithk';
					exp_file=[exp_directory '/' field '_AIS_' group '_' simul '_' expename '.nc'];
					mask_file=[exp_directory '/sftgif_AIS_' group '_' simul '_' expename '.nc'];
					if strcmp(expename,'exp11'),
						ctrl_file=[ctrl_directory '/' field '_AIS_' group '_' simul '_exp04.nc'];
						maskctrl_file=[ctrl_directory '/sftgif_AIS_' group '_' simul '_exp04.nc'];
					elseif strcmp(expename,'exp12'),
						ctrl_file=[ctrl_directory '/' field '_AIS_' group '_' simul '_exp08.nc'];
						maskctrl_file=[ctrl_directory '/sftgif_AIS_' group '_' simul '_exp08.nc'];
					end

					data = rot90(double(ncread(exp_file,field)));
					datac = rot90(double(ncread(ctrl_file,field)));
					mask = rot90(double(ncread(mask_file,'sftgif')));
					maskctrl= rot90(double(ncread(maskctrl_file,'sftgif')));
					if size(data)~=[761,761,21];
						error(['warming: file ' exp_file ' has the wrong size']);
					end
					data_init=data(:,:,end-85).*mask(:,:,end-85);
					data_end=data(:,:,end).*mask(:,:,end);
					datac_init=datac(:,:,end-85).*maskctrl(:,:,end-85);
					datac_end=datac(:,:,end).*maskctrl(:,:,1);
					pos=find(data_init~=0 & ~isnan(data_init) & data_end~=0 & ~isnan(data_end) & datac_init~=0 & ~isnan(datac_init) & datac_end~=0 & ~isnan(datac_end));
					nummodels(pos)=nummodels(pos)+1;
					totalthicknesschange(pos)=totalthicknesschange(pos)+data_end(pos)-data_init(pos)-(datac_end(pos)-datac_init(pos));

				end
			end
		end
	end
	mean_thicknesschange=totalthicknesschange./nummodels;
	pos=find(nummodels<5);
	mean_thicknesschange(pos)=NaN;
	collapsed_shelf=rot90(ncread('/u/astrid-r1b/seroussi/issmjpl/proj-seroussi/ISMIP6Projections/ForcingsPrep/OutputFiles/ice_shelf_collapse_mask_CCSM4_1995-2100_08km.nc','mask')); %ice shelf mask provided in ISMIP6 forcings at 8 km resolution
	collapsed_shelf_end=collapsed_shelf(:,:,end);
	pos=find(collapsed_shelf_end>0.5);
	mean_thicknesschange(pos)=NaN;

	%plot results
	close; set(gcf,'color','w'); set(gcf,'Position',[400 400 700 500]);
	[pos_nani pos_nanj]=find(isnan(mean_thicknesschange));
	data_min=-100; data_max=100;
	colorm = jet(100);
	%colormap used in the paper can be found here: https://www.mathworks.com/matlabcentral/fileexchange/17555-light-bartlein-color-maps
	image_rgb = ind2rgb(uint16((max(data_min,min(data_max,mean_thicknesschange)) - data_min)*(size(colorm,1)/(data_max-data_min))),colorm);
	image_rgb(sub2ind(size(image_rgb),repmat(pos_nani,1,3),repmat(pos_nanj,1,3),repmat(1:3,size(pos_nani,1),1))) = repmat([1 1 1],size(pos_nani,1),1);
	imagesc(-3040:6080/size(mean_thicknesschange,2):3040,-3040:6080/size(mean_thicknesschange,1):3040,image_rgb);
	colormap(jet); colorbar; caxis([data_min data_max]); set(gca,'fontsize',14); hcb=colorbar; title(hcb,'m');
	hold on
	x=repmat(-3040:6080/size(mean_thicknesschange,2):3040,length(-3040:6080/size(mean_thicknesschange,2):3040),1);
	y=x';
	x=(x(1:end-1,1:end-1)+x(2:end,2:end))/2;
	y=(y(1:end-1,1:end-1)+y(2:end,2:end))/2;
	[c1, h1]=contourf(x,y,collapsed_shelf_end,[0.5 0.5]);
	set(h1,'linestyle','none','Tag','HatchingRegion');
	ax1 = gca;
	ax1.XLim=[-3040 3032];
	ax1.YLim=[-3040 3032];
	hp = findobj(ax1,'Tag','HatchingRegion');
	hh = hatchfill2(hp,'single','LineWidth',1,'Fill','off','HatchDensity',150,'HatchColor','k');
	colorbar; axis('equal','off'); xlim([-3040 3040]); ylim([-3040 3040]);
	text(-2400,2400,'a','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',18,'fontweight','b');

	h = gcf;
	set(h,'Units','Inches');
	pos = get(h,'Position');
	set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
	print(gcf, '-dpdf', '-painters', 'FiguresPaperFinal/Figure14a.pdf');

end %}}}
if step==26 % {{{Figure 14b

	experiments_list={'exp11','exp12'}; %ice shelf collapse experiments
	nummodels=zeros(761,761);
	totalvelocitychange=zeros(761,761);
	velocitystd=zeros(761,761);

	for iexp=1:length(experiments_list),
		expename=experiments_list{iexp};

		for imodel=1:length(model_list),
			modelname=model_list{imodel};
			specifics;
			eval(['exp_name=' expename '_regrid;'])
			exp_directory=['' gridpath '/' group '/' simul '/' exp_name '/'];
			if strcmp(expename,'exp11'),
				eval(['ctrl_proj_name=exp04_regrid;'])
				ctrl_directory=['' gridpath '/' group '/' simul '/' ctrl_proj_name '/'];
			elseif strcmp(expename,'exp12'),
				eval(['ctrl_proj_name=exp08_regrid;'])
				ctrl_directory=['' gridpath '/' group '/' simul '/' ctrl_proj_name '/'];
			end

			eval(['isexp=is' expename ';'])
			if isexp,
				if exist(exp_directory,'dir')==0,
					error(['directory ' exp_directory ' not found']);
				end

			if is_xvelmean==1, 
				field='xvelmean'; field2='yvelmean';
			else error('should have some velocity fields');
			end
			exp_file=[exp_directory '/' field '_AIS_' group '_' simul '_' expename '.nc'];
			exp_file2=[exp_directory '/' field2 '_AIS_' group '_' simul '_' expename '.nc'];
			if strcmp(expename,'exp11'),
				ctrl_file=[ctrl_directory '/' field '_AIS_' group '_' simul '_exp04.nc'];
				ctrl_file2=[ctrl_directory '/' field2 '_AIS_' group '_' simul '_exp04.nc'];
			elseif strcmp(expename,'exp12'),
				ctrl_file=[ctrl_directory '/' field '_AIS_' group '_' simul '_exp08.nc'];
				ctrl_file2=[ctrl_directory '/' field2 '_AIS_' group '_' simul '_exp08.nc'];
			end

			datau = rot90(double(ncread(exp_file,field)));
			datav = rot90(double(ncread(exp_file2,field2)));
			data=sqrt(datau.^2+datav.^2)*31556926;
			datacu = rot90(double(ncread(ctrl_file,field)));
			datacv = rot90(double(ncread(ctrl_file2,field2)));
			datac=sqrt(datacu.^2+datacv.^2)*31556926;
			if size(data)~=[761,761,21];
				error(['warming: file ' abmb_file ' has the wrong size']);
			end
			data_init=data(:,:,1);
			data_end=data(:,:,end);
			datac_init=datac(:,:,1);
			datac_end=datac(:,:,end);
			pos=find(datac_init~=0 & ~isnan(datac_init) & datac_end~=0 & ~isnan(datac_end) & datac_init~=0 & ~isnan(datac_init) & datac_end~=0 & ~isnan(datac_end)); 
			nummodels(pos)=nummodels(pos)+1;
			totalvelocitychange(pos)=totalvelocitychange(pos)+(data_end(pos)-data_init(pos)) - (datac_end(pos)-datac_init(pos));

		end
	end
	end
	mean_velocitychange=totalvelocitychange./nummodels;
	pos=find(nummodels<5);
	mean_velocitychange(pos)=NaN;
	collapsed_shelf=rot90(ncread('/u/astrid-r1b/seroussi/issmjpl/proj-seroussi/ISMIP6Projections/ForcingsPrep/OutputFiles/ice_shelf_collapse_mask_CCSM4_1995-2100_08km.nc','mask')); %ice shelf mask provided in ISMIP6 forcings at 8 km resolution
	collapsed_shelf_end=collapsed_shelf(:,:,end);
	pos=find(collapsed_shelf_end>0.5);
	mean_velocitychange(pos)=NaN;

	%plot results
	set(gcf,'color','w'); set(gcf,'Position',[400 400 700 500]);
	[pos_nani pos_nanj]=find(isnan(mean_velocitychange));
	data_min=-100; data_max=100;
	colorm = jet(100);
	%colormap used in the paper can be found here: https://www.mathworks.com/matlabcentral/fileexchange/17555-light-bartlein-color-maps
	image_rgb = ind2rgb(uint16((max(data_min,min(data_max,mean_velocitychange)) - data_min)*(size(colorm,1)/(data_max-data_min))),colorm);
	image_rgb(sub2ind(size(image_rgb),repmat(pos_nani,1,3),repmat(pos_nanj,1,3),repmat(1:3,size(pos_nani,1),1))) = repmat([1 1 1],size(pos_nani,1),1);
	imagesc(-3040:6080/size(mean_velocitychange,2):3040,-3040:6080/size(mean_velocitychange,1):3040,image_rgb); 
	hold on
	x=repmat(-3040:6080/size(mean_velocitychange,2):3040,length(-3040:6080/size(mean_velocitychange,2):3040),1);
	y=x';
	x=(x(1:end-1,1:end-1)+x(2:end,2:end))/2;
	y=(y(1:end-1,1:end-1)+y(2:end,2:end))/2;
	[c1, h1]=contourf(x,y,collapsed_shelf_end,[0.5 0.5]);
	set(h1,'linestyle','none','Tag','HatchingRegion');
	ax1 = gca;
	ax1.XLim=[-3040 3032];
	ax1.YLim=[-3040 3032];
	hp = findobj(ax1,'Tag','HatchingRegion');
	hh = hatchfill2(hp,'single','LineWidth',1,'Fill','off','HatchDensity',150,'HatchColor','k');
	axis('equal','off'); xlim([-3040 3040]); ylim([-3040 3040]);
	colormap(jet); colorbar; caxis([data_min data_max]); set(gca,'fontsize',14); hcb=colorbar; title(hcb,'m/yr');
	text(-2400,2400,'b','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',18,'fontweight','b');

	h = gcf;
	set(h,'Units','Inches');
	pos = get(h,'Position');
	set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
	print(gcf, '-dpdf', '-painters', 'FiguresPaperFinal/Figure14b.pdf');
end %}}}

%Figure 15: regional sensitivity
if step==27 % {{{Figure 15

	colors = distinguishable_colors(18); %number of basins
	experiments_list={'exp01','exp02','exp04','exp05','exp06','exp08','expA1','expA2','expA3','expA5','expA6','expA7'}; %open and standard RCP8.5

	for isector=1:18,
		eval(['sensitivity_results_limnsw_sector' num2str(isector) '=[];'])
		eval(['sensitivity_results_tendacabfgr_sector' num2str(isector) '=[];'])
		eval(['sensitivity_results_tendlibmassbffl_sector' num2str(isector) '=[];'])
	end

	for imodel=1:length(model_list),
		modelname=model_list{imodel};
		for iexp=1:length(experiments_list),
			expename=experiments_list{iexp};
			specifics;

			eval(['isexp=is' expename ';'])
			if isexp,
				explimnsw_file=['' scalarpath '/' group '/' simul '/' expename '/computed_ivaf_minus_ctrl_proj_AIS_' group '_' simul '_' expename '.nc'];
				exptendacabfgr_file=['' scalarpath '/' group '/' simul '/' expename '/computed_smbgr_minus_ctrl_proj_AIS_' group '_' simul '_' expename '.nc'];
				exptendlibmassbffl_file=['' scalarpath '/' group '/' simul '/' expename '/computed_bmbfl_minus_ctrl_proj_AIS_' group '_' simul '_' expename '.nc'];

				for isector=1:18,
					eval(['results_tendlibmassbffl_sector' num2str(isector) '=ncread(exptendlibmassbffl_file,''bmbfl_sector_' num2str(isector) ''') *yearday_model*3600*24/(10^9*1000);']) %from kg/s to Gt/yr
					eval(['sensitivity_results_tendlibmassbffl_sector' num2str(isector) '(end+1)=sum(results_tendlibmassbffl_sector' num2str(isector) '(1:85));'])

					eval(['results_tendacabfgr_sector' num2str(isector) '=ncread(exptendacabfgr_file,''smbgr_sector_' num2str(isector) ''')*yearday_model*3600*24/(10^9*1000);']) %from kg/s to Gt/yr
					eval(['sensitivity_results_tendacabfgr_sector' num2str(isector) '(end+1)=sum(results_tendacabfgr_sector' num2str(isector) '(1:85));'])

					eval(['results_limnsw_sector' num2str(isector) '=ncread(explimnsw_file,''ivaf_sector_' num2str(isector) ''') *ice_density/(10^9*1000);']) %from m^3 to Gt
					eval(['sensitivity_results_limnsw_sector' num2str(isector) '(end+1)=results_limnsw_sector' num2str(isector) '(85);'])
				end
			end

		end %end of model
	end %end of isexp

	%Plot results per region
	figure(1); set(gcf,'color','w'); set(gcf,'Position',[400 500 800 500]);
	for ibasin=1:18,
		eval(['plot(-sensitivity_results_tendlibmassbffl_sector' num2str(ibasin) ',-(sensitivity_results_limnsw_sector' num2str(ibasin) '-sensitivity_results_tendacabfgr_sector' num2str(ibasin) ')/362.5,''.'',''color'',colors(ibasin,:));'])
		hold on
	end
	ylabel('Dynamic mass loss (mm SLE)');
	xlabel('Ice shelf melt (Gt)');
	xlim([-0.1 1.6]*10^5)
	ylim([-20 140])
	text(-1.2*10^4,-28,'b','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',18,'fontweight','b');
	h = gcf;
	set(h,'Units','Inches');
	pos = get(h,'Position');
	set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
	print(gcf, '-dpdf', '-painters', 'FiguresPaperFinal/Figure15.pdf');

	%Now create inset map to show the color of each region
	colorm = distinguishable_colors(18);
	sectors_8km=flipud(rot90(double(ncread('Data/sectors_8km_iceonly.nc','sectors')))); %18 sectors with ice regions only
	set(gcf,'color','w');
	[pos_nani pos_nanj]=find(isnan(sectors_8km) | sectors_8km==0);
	data_min=1; data_max=19;
	image_rgb = ind2rgb(uint16((max(data_min,min(data_max,sectors_8km)) - data_min)*(size(colorm,1)/(data_max-data_min))),colorm);
	image_rgb(sub2ind(size(image_rgb),repmat(pos_nani,1,3),repmat(pos_nanj,1,3),repmat(1:3,size(pos_nani,1),1))) = repmat([1 1 1],size(pos_nani,1),1);
	imagesc(-3040:6080/size(sectors_8km,2):3040,-3040:6080/size(sectors_8km,1):3040,image_rgb); colorbar('off');
	set(gca,'fontsize',14);
	colorbar('off'); axis('equal','off'); xlim([-3040 3040]); ylim([-3040 3040]);
	h = gcf;
	set(h,'Units','Inches');
	pos = get(h,'Position');
	set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
	print(gcf, '-dpdf', '-painters', 'FiguresPaperFinal/Figure15inset.pdf');

end %}}}
