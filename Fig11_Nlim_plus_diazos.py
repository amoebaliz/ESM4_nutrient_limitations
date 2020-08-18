import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

#% mapping routines and color maps
#addpath /home/cas/matlab/m_map
#addpath /home/cas/matlab/kakearney-cptcmap-pkg-845bf83/cptcmap/
#addpath /home/cas/matlab/kakearney-cptcmap-pkg-845bf83/cptcmap/cptfiles/

arch_dir = '/archive/oar.gfdl.cmip6/ESM4/DECK/ESM4_historical_D1/gfdl.ncrc4-intel16-prod-openmp/pp/ocean_cobalt_omip_2d/av/annual_5yr/'

#NOTE: are these two files spatially different?
mod_file = arch_dir + 'ocean_cobalt_omip_2d.1995-1999.ann.nc'
grid_file = '/archive/oar.gfdl.cmip6/ESM4/DECK/ESM4_historical_D1/gfdl.ncrc4-intel16-prod-openmp/pp/ocean_annual/ocean_annual.static.nc';
fid_mod = nc.Dataset(mod_file)
fid_grd = nc.Dataset(grid_file)

#diminfo = nc_getdiminfo(mod_file,'yh'); nmod_lat = diminfo.Length;
#diminfo = nc_getdiminfo(mod_file,'xh'); nmod_lon = diminfo.Length;

#lon_mod = nc_varget(grid_file,'geolon',[0 0],[nmod_lat nmod_lon]);
#lat_mod = nc_varget(grid_file,'geolat',[0 0],[nmod_lat nmod_lon]);

lon_mod = fid_grd.variables['geolon'][:]
lat_mod = fid_grd.variables['geolat'][:]
area    = fid_grd.variables['areacello'][:]

nmod_lat, nmod_lon = lat_mod.shape

#intpp_mod = squeeze(nc_varget(mod_file,'intpp',[0 0 0],[1 nmod_lat nmod_lon]));
#intpp_mod = squeeze(intpp_mod)*86400*1000*12;

intpp_mod = 86400*1000*12*fid_mod.variables['intpp'][0,:].squeeze()

#aa = find(isfinite(intpp_mod));
aa = np.isfinite(intpp)

# calculate the total primary production (Pg C yr-1)
#intpp_tot = sum(intpp_mod(aa).*area(aa))/1e15;
intpp_tot = np.dot(intpp_mod[aa],area[aa])/1e15

num_days = [31 28 31 30 31 30 31 31 30 31 30 31]

nlim_weighted_month  = np.zeros((12,nmod_lat,nmod_lon))
plim_weighted_month  = np.zeros((12,nmod_lat,nmod_lon))
felim_weighted_month = np.zeros(12,nmod_lat,nmod_lon))

pp_month = np.zeros((12,nmod_lat,nmod_lon))

#for m = 1:12
#    if m < 10
#        mod_file1 = ['/archive/oar.gfdl.cmip6/ESM4/DECK/ESM4_historical_D1/gfdl.ncrc4-intel16-prod-openmp/pp/ocean_cobalt_omip_2d/av/monthly_5yr/ocean_cobalt_omip_2d.1995-1999.0',num2str(m),'.nc'];
#        mod_file2 = ['/archive/oar.gfdl.cmip6/ESM4/DECK/ESM4_historical_D1/gfdl.ncrc4-intel16-prod-openmp/pp/ocean_cobalt_omip_2d/av/monthly_5yr/ocean_cobalt_omip_2d.2000-2004.0',num2str(m),'.nc'];
#        mod_file3 = ['/archive/oar.gfdl.cmip6/ESM4/DECK/ESM4_historical_D1/gfdl.ncrc4-intel16-prod-openmp/pp/ocean_cobalt_omip_2d/av/monthly_5yr/ocean_cobalt_omip_2d.2005-2009.0',num2str(m),'.nc'];
#        mod_file4 = ['/archive/oar.gfdl.cmip6/ESM4/DECK/ESM4_historical_D1/gfdl.ncrc4-intel16-prod-openmp/pp/ocean_cobalt_omip_2d/av/monthly_5yr/ocean_cobalt_omip_2d.2010-2014.0',num2str(m),'.nc'];
#    else
#        mod_file1 = ['/archive/oar.gfdl.cmip6/ESM4/DECK/ESM4_historical_D1/gfdl.ncrc4-intel16-prod-openmp/pp/ocean_cobalt_omip_2d/av/monthly_5yr/ocean_cobalt_omip_2d.1995-1999.',num2str(m),'.nc'];
#        mod_file2 = ['/archive/oar.gfdl.cmip6/ESM4/DECK/ESM4_historical_D1/gfdl.ncrc4-intel16-prod-openmp/pp/ocean_cobalt_omip_2d/av/monthly_5yr/ocean_cobalt_omip_2d.2000-2004.',num2str(m),'.nc'];
#        mod_file3 = ['/archive/oar.gfdl.cmip6/ESM4/DECK/ESM4_historical_D1/gfdl.ncrc4-intel16-prod-openmp/pp/ocean_cobalt_omip_2d/av/monthly_5yr/ocean_cobalt_omip_2d.2005-2009.',num2str(m),'.nc'];
#        mod_file4 = ['/archive/oar.gfdl.cmip6/ESM4/DECK/ESM4_historical_D1/gfdl.ncrc4-intel16-prod-openmp/pp/ocean_cobalt_omip_2d/av/monthly_5yr/ocean_cobalt_omip_2d.2010-2014.',num2str(m),'.nc'];
#    end

for m in range(12):
    # create list of files
    y1 = 1995
    nc_fils = []
    for fils in range(4):
        nc_fil = arch_dir + 'ocean_cobalt_omip_2d.' + str(y1) + '-' + str(y1+4) + '.' + str(m).zfill(2) + '.nc'
        nc_fils.append(nc_fil)
        y1+=5
      
#    % get constituent productivity terms
#    intppdiat1 = squeeze(nc_varget(mod_file1,'intppdiat',[0 0 0],[1 nmod_lat nmod_lon]));
#    intppdiat2 = squeeze(nc_varget(mod_file2,'intppdiat',[0 0 0],[1 nmod_lat nmod_lon]));
#    intppdiat3 = squeeze(nc_varget(mod_file3,'intppdiat',[0 0 0],[1 nmod_lat nmod_lon]));
#    intppdiat4 = squeeze(nc_varget(mod_file4,'intppdiat',[0 0 0],[1 nmod_lat nmod_lon]));
#    intppdiat = (intppdiat1 + intppdiat2 + intppdiat3 + intppdiat4)/4;

#    intppmisc1 = squeeze(nc_varget(mod_file1,'intppmisc',[0 0 0],[1 nmod_lat nmod_lon]));
#    intppmisc2 = squeeze(nc_varget(mod_file2,'intppmisc',[0 0 0],[1 nmod_lat nmod_lon]));
#    intppmisc3 = squeeze(nc_varget(mod_file3,'intppmisc',[0 0 0],[1 nmod_lat nmod_lon]));
#    intppmisc4 = squeeze(nc_varget(mod_file4,'intppmisc',[0 0 0],[1 nmod_lat nmod_lon]));
#    intppmisc = (intppmisc1 + intppmisc2 + intppmisc3 + intppmisc4)/4;

#    intpppico1 = squeeze(nc_varget(mod_file1,'intpppico',[0 0 0],[1 nmod_lat nmod_lon]));
#    intpppico2 = squeeze(nc_varget(mod_file2,'intpppico',[0 0 0],[1 nmod_lat nmod_lon]));
#    intpppico3 = squeeze(nc_varget(mod_file3,'intpppico',[0 0 0],[1 nmod_lat nmod_lon]));
#    intpppico4 = squeeze(nc_varget(mod_file4,'intpppico',[0 0 0],[1 nmod_lat nmod_lon]));
#    intpppico = (intpppico1 + intpppico2 + intpppico3 + intpppico4)/4;

# Limitation terms (picophytoplankton)
#    limnpico1 = squeeze(nc_varget(mod_file1,'limnpico',[0 0 0],[1 nmod_lat nmod_lon]));
#    limnpico2 = squeeze(nc_varget(mod_file2,'limnpico',[0 0 0],[1 nmod_lat nmod_lon]));
#    limnpico3 = squeeze(nc_varget(mod_file3,'limnpico',[0 0 0],[1 nmod_lat nmod_lon]));
#    limnpico4 = squeeze(nc_varget(mod_file4,'limnpico',[0 0 0],[1 nmod_lat nmod_lon]));
#    limnpico = (limnpico1 + limnpico2 + limnpico3 + limnpico4)/4;
#    limppico1 = squeeze(nc_varget(mod_file1,'limppico',[0 0 0],[1 nmod_lat nmod_lon]));
#    limppico2 = squeeze(nc_varget(mod_file2,'limppico',[0 0 0],[1 nmod_lat nmod_lon]));
#    limppico3 = squeeze(nc_varget(mod_file3,'limppico',[0 0 0],[1 nmod_lat nmod_lon]));
#    limppico4 = squeeze(nc_varget(mod_file4,'limppico',[0 0 0],[1 nmod_lat nmod_lon]));
#    limppico = (limppico1 + limppico2 + limppico3 + limppico4)/4;
#    limfepico1 = squeeze(nc_varget(mod_file1,'limfepico',[0 0 0],[1 nmod_lat nmod_lon]));
#    limfepico2 = squeeze(nc_varget(mod_file2,'limfepico',[0 0 0],[1 nmod_lat nmod_lon]));
#    limfepico3 = squeeze(nc_varget(mod_file3,'limfepico',[0 0 0],[1 nmod_lat nmod_lon]));
#    limfepico4 = squeeze(nc_varget(mod_file4,'limfepico',[0 0 0],[1 nmod_lat nmod_lon]));
#    limfepico = (limfepico1 + limfepico2 + limfepico3 + limfepico4)/4;

#    % diatoms
#    limndiat1 = squeeze(nc_varget(mod_file1,'limndiat',[0 0 0],[1 nmod_lat nmod_lon]));
#    limndiat2 = squeeze(nc_varget(mod_file2,'limndiat',[0 0 0],[1 nmod_lat nmod_lon]));
#    limndiat3 = squeeze(nc_varget(mod_file3,'limndiat',[0 0 0],[1 nmod_lat nmod_lon]));
#    limndiat4 = squeeze(nc_varget(mod_file4,'limndiat',[0 0 0],[1 nmod_lat nmod_lon]));
#    limndiat = (limndiat1 + limndiat2 + limndiat3 + limndiat4)/4;
#    limpdiat1 = squeeze(nc_varget(mod_file1,'limpdiat',[0 0 0],[1 nmod_lat nmod_lon]));
#    limpdiat2 = squeeze(nc_varget(mod_file2,'limpdiat',[0 0 0],[1 nmod_lat nmod_lon]));
#    limpdiat3 = squeeze(nc_varget(mod_file3,'limpdiat',[0 0 0],[1 nmod_lat nmod_lon]));
#    limpdiat4 = squeeze(nc_varget(mod_file4,'limpdiat',[0 0 0],[1 nmod_lat nmod_lon]));
#    limpdiat = (limpdiat1 + limpdiat2 + limpdiat3 + limpdiat4)/4;
#    limfediat1 = squeeze(nc_varget(mod_file1,'limfediat',[0 0 0],[1 nmod_lat nmod_lon]));
#    limfediat2 = squeeze(nc_varget(mod_file2,'limfediat',[0 0 0],[1 nmod_lat nmod_lon]));
#    limfediat3 = squeeze(nc_varget(mod_file3,'limfediat',[0 0 0],[1 nmod_lat nmod_lon]));
#    limfediat4 = squeeze(nc_varget(mod_file4,'limfediat',[0 0 0],[1 nmod_lat nmod_lon]));
#    limfediat = (limfediat1 + limfediat2 + limfediat3 + limfediat4)/4;

    # Miscellaneous
#    limnmisc1 = squeeze(nc_varget(mod_file1,'limnmisc',[0 0 0],[1 nmod_lat nmod_lon]));
#    limnmisc2 = squeeze(nc_varget(mod_file2,'limnmisc',[0 0 0],[1 nmod_lat nmod_lon]));
#    limnmisc3 = squeeze(nc_varget(mod_file3,'limnmisc',[0 0 0],[1 nmod_lat nmod_lon]));
#    limnmisc4 = squeeze(nc_varget(mod_file4,'limnmisc',[0 0 0],[1 nmod_lat nmod_lon]));
#    limnmisc = (limnmisc1 + limnmisc2 + limnmisc3 + limnmisc4)/4;
#    limpmisc1 = squeeze(nc_varget(mod_file1,'limpmisc',[0 0 0],[1 nmod_lat nmod_lon]));
#    limpmisc2 = squeeze(nc_varget(mod_file2,'limpmisc',[0 0 0],[1 nmod_lat nmod_lon]));
#    limpmisc3 = squeeze(nc_varget(mod_file3,'limpmisc',[0 0 0],[1 nmod_lat nmod_lon]));
#    limpmisc4 = squeeze(nc_varget(mod_file4,'limpmisc',[0 0 0],[1 nmod_lat nmod_lon]));
#    limpmisc = (limpmisc1 + limpmisc2 + limpmisc3 + limpmisc4)/4;
#    limfemisc1 = squeeze(nc_varget(mod_file1,'limfemisc',[0 0 0],[1 nmod_lat nmod_lon]));
#    limfemisc2 = squeeze(nc_varget(mod_file2,'limfemisc',[0 0 0],[1 nmod_lat nmod_lon]));
#    limfemisc3 = squeeze(nc_varget(mod_file3,'limfemisc',[0 0 0],[1 nmod_lat nmod_lon]));
#    limfemisc4 = squeeze(nc_varget(mod_file4,'limfemisc',[0 0 0],[1 nmod_lat nmod_lon]));
#    limfemisc = (limfemisc1 + limfemisc2 + limfemisc3 + limfemisc4)/4;

#    nlim_weighted_month(m,:,:) = (intppdiat.*limndiat+intppmisc.*limnmisc+intpppico.*limnpico)./ ...
#        (intppdiat + intppmisc + intpppico);

#    plim_weighted_month(m,:,:) = (intppdiat.*limpdiat+intppmisc.*limpmisc+intpppico.*limppico)./ ...
#        (intppdiat + intppmisc + intpppico);

#    felim_weighted_month(m,:,:) = (intppdiat.*limfediat+intppmisc.*limfemisc+intpppico.*limfepico)./ ...
#        (intppdiat + intppmisc + intpppico);

#    pp_month(m,:,:) = (intppdiat + intpppico + intppmisc)*86400*num_days(m)*6.625*12;
# end

    fid_mod = nc.MFDataset(nc_fils) 

    # Constituent productivity terms
    intppdiat = np.mean(fid_mod.variables['intppdiat'][:].squeeze(),axis=0)    
    intppmisc = np.mean(fid_mod.variables['intppmisc'][:].squeeze(),axis=0)
    intpppico = np.mean(fid_mod.variables['intpppico'][:].squeeze(),axis=0)   

    # Limitation terms (picophytoplankton)
    limnpico  = np.mean(fid_mod.variables['limnpico'][:].squeeze(),axis=0)
    limppico  = np.mean(fid_mod.variables['limppico'][:].squeeze(),axis=0)
    limfepico = np.mean(fid_mod.variables['limfepico'][:].squeeze(),axis=0)

    # Diatoms
    limndiat  = np.mean(fid_mod.variables['limndiat'][:].squeeze(),axis=0)
    limpdiat  = np.mean(fid_mod.variables['limpdiat'][:].squeeze(),axis=0)
    limfediat = np.mean(fid_mod.variables['limfediat'][:].squeeze(),axis=0)

    # Miscellaneous
    limnmisc  = np.mean(fid_mod.variables['limnmisc'][:].squeeze(),axis=0)
    limpmisc  = np.mean(fid_mod.variables['limpmisc'][:].squeeze(),axis=0)
    limfemisc = np.mean(fid_mod.variables['limfemisc'][:].squeeze(),axis=0)

    nlim_weighted_month[m,:]  = (intppdiat*limndiat + intppmisc*limnmisc + intpppico*limnpico)/\
                                (intppdiat + intppmisc + intpppico)

    plim_weighted_month[m,:]  = (intppdiat*limpdiat + intppmisc*limpmisc + intpppico*limppico)/\
                                (intppdiat + intppmisc + intpppico);

    felim_weighted_month[m,:] = (intppdiat*limfediat + intppmisc*limfemisc + intpppico*limfepico)/\
                                (intppdiat + intppmisc + intpppico);
    
    pp_month[m,:] = (intppdiat + intpppico + intppmisc)*86400*num_days[m]*6.625*12

# Monthly weighted mean
sum_pp_month = np.sum(pp_month,axis=0)

nlim_weighted  = np.sum(nlim_weighted_month*pp_month,axis=0)/sum_pp_month
plim_weighted  = np.sum(plim_weighted_month*pp_month,axis=0)/sum_pp_month
felim_weighted = np.sum(felim_weighted_month*pp_month,axis=0)/sum_pp_month
    
#felim   = np.zeros((nmod_lat,nmod_lon))
#nlim    = np.zeros((nmod_lat,nmod_lon))
#plim    = np.zeros((nmod_lat,nmod_lon))
#lim_ind = np.zeros((nmod_lat,nmod_lon))

#for i = 1:nmod_lat
#    for j = 1:nmod_lon
#        if isfinite(felim_weighted(i,j));
#            felim(i,j) = felim_weighted(i,j) - min(nlim_weighted(i,j),plim_weighted(i,j));
#            nlim(i,j) = nlim_weighted(i,j) - min(felim_weighted(i,j),plim_weighted(i,j));
#            plim(i,j) = plim_weighted(i,j) - min(nlim_weighted(i,j),felim_weighted(i,j));
#        else
#            felim(i,j) = NaN;
#            nlim(i,j) = NaN;
#            plim(i,j) = NaN;
#        end
#    end
#end

#liebig_macro = min(nlim_weighted,plim_weighted);
#liebig = min(liebig_macro,felim_weighted);

aa = ~isfinite(felim_weighted)

felim = felim_weighted - np.min(np.stack((nlim_weighted,plim_weighted)),axis=0)
felim[aa] = np.nan 
nlim  = nlim_weighted  - np.min(np.stack((felim_weighted,plim_weighted)),axis=0)
nlim[aa] = np.nan
plim  = plim_weighted  - np.min(np.stack((nlim_weighted,felim_weighted)),axis=0)
plim[aa] = np.nan

liebig_macro = np.min(np.stack((nlim_weighted,plim_weighted)),axis=0)
liebig = np.min(np.stack((liebig_macro,felim_weighted)),axis=0)

lim_ind[nlim <= 0] = 1;
lim_ind[(plim <= 0 and plim > -0.25)] = 1.5
lim_ind[plim <= -0.25] = 2
lim_ind[(felim <= 0 and felim > -0.25)] = 4.5
lim_ind[felim <= -0.25] = 5

# Include if you'd like to see areas where nutrients are generally not
# limiting
# lim_ind(liebig > 0.9) = NaN;

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% Get data for plotting nitrogen fixation                                 %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#mod_file1 = '/archive/oar.gfdl.cmip6/ESM4/DECK/ESM4_historical_D1/gfdl.ncrc4-intel16-prod-openmp/pp/ocean_cobalt_omip_2d/av/annual_5yr/ocean_cobalt_omip_2d.1995-1999.ann.nc';
#mod_file2 = '/archive/oar.gfdl.cmip6/ESM4/DECK/ESM4_historical_D1/gfdl.ncrc4-intel16-prod-openmp/pp/ocean_cobalt_omip_2d/av/annual_5yr/ocean_cobalt_omip_2d.2000-2004.ann.nc';
#mod_file3 = '/archive/oar.gfdl.cmip6/ESM4/DECK/ESM4_historical_D1/gfdl.ncrc4-intel16-prod-openmp/pp/ocean_cobalt_omip_2d/av/annual_5yr/ocean_cobalt_omip_2d.2005-2009.ann.nc';
#mod_file4 = '/archive/oar.gfdl.cmip6/ESM4/DECK/ESM4_historical_D1/gfdl.ncrc4-intel16-prod-openmp/pp/ocean_cobalt_omip_2d/av/annual_5yr/ocean_cobalt_omip_2d.2010-2014.ann.nc';

#intpn2_1 = squeeze(nc_varget(mod_file1,'intpn2',[0 0 0],[1 nmod_lat nmod_lon]));
#intpn2_2 = squeeze(nc_varget(mod_file2,'intpn2',[0 0 0],[1 nmod_lat nmod_lon]));
#intpn2_3 = squeeze(nc_varget(mod_file3,'intpn2',[0 0 0],[1 nmod_lat nmod_lon]));
#intpn2_4 = squeeze(nc_varget(mod_file4,'intpn2',[0 0 0],[1 nmod_lat nmod_lon]));
#intpn2 = squeeze(intpn2_1+intpn2_2+intpn2_3+intpn2_4)/4;
#intpn2 = intpn2*86400*1000; % mmoles day-1
#aa = find(isfinite(intpn2));
# calculate the total primary production (Tg N yr-1)
#intpn2_tot = sum(intpn2(aa).*area(aa)*365*14/1000/1e12);

y1 = 1995
nc_fils = []
for fils in range(4):
    nc_fil = arch_dir + 'ocean_cobalt_omip_2d.' + str(y1) + '-' + str(y1+4) + '.ann.nc'
        nc_fils.append(nc_fil)
        y1+=5

fid_mod = nc.MFDataset(nc_fils)
intpn2 = 86400*1000*np.mean(fid_mod.variables['intpn2'][:].squeeze(),axis=0)
aa = np.isfinite(intpn2)

# calculate the total primary production (Tg N yr-1)
intpn2_tot = np.dot((intpn2[aa],area[aa]))*365*14/(1000*1e12)

width = 0.8
height = 0.4
lfs = 10
fs = 10

###################
# Recreate figure #
###################
cmap1 = cptcmap('GMT_no_green');
cmap = interp1([1:4:61],cmap1,1:61);

figure(1);
clf

axes('position',[0.05 0.52 width height]);
hold on
title('Phytoplankton Limiting Nutrient');
colormap(cmap);
m_proj('Mollweide','long',[min(lon_mod(:)) max(lon_mod(:))],'lat',[min(lat_mod(:)) max(lat_mod(:))]);
m_pcolor(lon_mod,lat_mod,lim_ind);
shading flat
caxis([1 5]);
m_coast('patch',[0.7 0.7 0.7],'edgecolor','k');
m_grid('fontsize',fs,'xticklabels',[],'yticklabels',[]);
axis off

axes('position',[0.825 0.52 0.15 height]);
hold on
set(gca,'Xlim',[0 1],'Ylim',[0 1]);
patch([0.1 0.3 0.3 0.1],[0.05 0.05 0.15 0.15],5);
text(0.35,0.1,'Fe','FontSize',fs);
patch([0.1 0.3 0.3 0.1],[0.25 0.25 0.35 0.35],4.5);
text(0.35,0.3,'weakly Fe','FontSize',fs);
patch([0.1 0.3 0.3 0.1],[0.45 0.45 0.55 0.55],2);
text(0.35,0.5,'P','FontSize',fs)
patch([0.1 0.3 0.3 0.1],[0.65 0.65 0.75 0.75],1.5);
text(0.35,0.7,'weakly P','FontSize',fs)
patch([0.1 0.3 0.3 0.1],[0.85 0.85 0.95 0.95],1);
text(0.35,0.9,'N','FontSize',fs)
caxis([1 5]);
axis off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot nitrogen fixation                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes('position',[0.05 0.05 width height]);
hold on
cmin = 0; cmax = 0.3; cint = 0.0025;

colormap(cmap);
m_proj('Mollweide','long',[min(lon_mod(:)) max(lon_mod(:))],'lat',[min(lat_mod(:)) max(lat_mod(:))]);
m_pcolor(lon_mod,lat_mod,intpn2);
title('Nitrogen Fixation, mmoles m^{-2} day^{-1}');
shading flat
caxis([cmin cmax]);
m_coast('line','Color',[0 0 0]);
m_grid('fontsize',fs,'xticklabels',[],'yticklabels',[]);
axis off;

axes('position',[0.85 0.1 0.03 0.3]);
hold on
caxis([cmin cmax]);
h = colorbar('East');
p = get(h,'position');
set(h,'YAxisLocation','right');
set(h,'position',[0.85 0.1 0.03 0.3]);
axis off;

set(gcf,'PaperPosition',[1.25 3 6 5]);
fname = 'Fig11_NPP_limitation_Nfix_300dpi';
print('-r300','-dpng',fname);
print('-r300','-dtiff',fname);

