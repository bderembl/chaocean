%ncks -d MT,0,20 uvts.nc uvts_small.nc
% usefull trick
% i=find(a==0);a(i)=nan;
clear;

ieee='b';
accuracy='real*4';

% output dir
dir_o = '/tank/groups/climode/chaocean/init_cond97_12/';

% ocean input dir
dir_1_init = '/tank/chaocean/initial_data/';
dir_1_obcs = '/tank/chaocean/boundary_data/';

% atm input dir
dir_era = '/tank/groups/climode/chaocean/era/';

build_init = 1; % 1: generate initial conditions, 0: no
build_obcs = 1; % 1: generate open boundaries files, 0: no

build_atmosphere = 1; % 1: generate atmospheric forcing files, 0: no

flag_cut = 1;  % 1: put eastern atlantic in south america, 0: no
flag_interp = 0; % 1: matlab interpolation (slower), 0: mygriddata
flag_plot = 0; % 1: plot a couple of diags, 0: no

% ================== NEW GRID =====================================

% horizontal resolution
dla = 1/12.;
dlo = 1/11.62;

%dla = 1/4.;
%dlo = 1/4.;
  
la1 = -20+dla;
la2 = 55.0;
lo1 = -98;
lo2 = 14-dlo;

yy = [la1:dla:la2];
xx = [lo1:dlo:lo2];
[yla,xlo]=meshgrid(yy,xx);
[si_x_mit,si_y_mit] = size(xlo);

dla_km = dla*2*pi*6370/360;
dlo_km = dlo*2*pi*6370*cos(2*pi*yy/360)/360;

dx_fi = dlo*ones(si_x_mit,1);
dy_fi = dla*ones(si_y_mit,1);

% vertical grid
%dz_fi = [5,5,5,5,5,5,10,10,10,10,10,10,10,10,15,15,15,15,15,15,20,20,20,20,20,20,30,30,30,30,30,30,40,40,40,40,50,50,50,60,60,60,80,80,80,100,100,100,100,120,120,120,120,140,140,140,140,160,160,160,160,180,180,180,180,200,200,200,200,200,200,200,200];

% 46 layers
dz_fi=[3.0468, 6.9099,7.5347,8.3106,9.2726,10.4633,11.9339,13.7455,15.9696,18.6894,21.9988,26.0017,30.8085,36.5308,43.2729,51.1207,60.1275,70.2985,81.5768,93.8328,106.8625,120.3962,134.1177, 147.6931,  160.8020,  173.1690,  184.5790,  194.8930,  204.0490,  212.0400,  218.9170,  224.7620,  229.6790,  233.7760,  237.1680,  239.9560,  242.2380,  244.0970,  245.6070,  246.8310,  247.8180,  248.6160,  249.2580,  249.7750,  250.1890,  250.0000];


%dz_fi = [500,1000];

si_z_mit = size(dz_fi,2);

lay_de = zeros(si_z_mit,1);

lay_de(1) = dz_fi(1)/2;
for nz = 2:si_z_mit
  lay_de(nz) = lay_de(nz-1)+dz_fi(nz-1)/2+dz_fi(nz)/2;
end

fid=fopen([dir_o,'dx.box'],'w',ieee);  fwrite(fid,dx_fi,accuracy); fclose(fid);
fid=fopen([dir_o,'dy.box'],'w',ieee);  fwrite(fid,dy_fi,accuracy); fclose(fid);
fid=fopen([dir_o,'dz.box'],'w',ieee);  fwrite(fid,dz_fi,accuracy); fclose(fid);

fprintf('MIT max depth %f m \n',sum(dz_fi));
fprintf('si_x: %i \n',si_x_mit);
fprintf('si_y: %i \n',si_y_mit);
fprintf('si_z: %i \n',si_z_mit);

% ybc : variable for north-south BC cut
ybc = 1;
if flag_cut
  % which lat to cut
  lat_cut = -12.;
  [argvalue, x_cut] = min(abs(xx-lat_cut));

  if x_cut > 0.5*si_x_mit
    x_c1 = x_cut;
    x_c2 = si_x_mit;
    % cut west
    flag_cut = 2;
  else
    x_c1 = 1;
    x_c2 = x_cut;
    % cut east
    flag_cut = 1;
  end    
  
  si_x_mit2 = si_x_mit-(x_c2-x_c1);
  
  % if it is a regular grid: doesn't matter which dx you pick
  fid=fopen([dir_o,'dx.box'],'w',ieee);  fwrite(fid,dx_fi(1:si_x_mit2),accuracy); fclose(fid);  
  fprintf('*******************\n');
  fprintf('new si_x: %i \n',si_x_mit2);
  fprintf('lon cut: %6.2f \n',xx(x_cut+1));
  fprintf('*******************\n');
end

fprintf('Building oceanic files \n');
  
% ================ VERTICAL GRID =========================================

fileT = '1987/CI-MJM88_y1987m01d05_gridT.nc';


dep = ncread([dir_1_init,fileT],'deptht');
[si_z,naux] = size(dep);

%find weights for vertical interpolation
w_z = zeros(si_z_mit,3);
for nz = 1:si_z_mit
  nz2 = 1;
  while dep(nz2)<lay_de(nz)
    nz2 = nz2 + 1;
    if nz2>si_z; break; end;
  end
  if nz2>1 && nz2 <= si_z
    w_z(nz,1) = nz2-1; % no of the layer
    w_z(nz,2) = (dep(nz2)-lay_de(nz))/(dep(nz2)-dep(nz2-1)); % weight for the upper layer
    w_z(nz,3) = (lay_de(nz)-dep(nz2-1))/(dep(nz2)-dep(nz2-1)); % weight for the lower layer
  elseif nz2 == 1
    w_z(nz,1) = 1; % no of the layer
    w_z(nz,2) = 1; % weight for the upper layer
    w_z(nz,3) = 0; % weight for the lower layer
  elseif nz2 > si_z
    w_z(nz,1) = si_z-1; % no of the layer
    w_z(nz,2) = 0; % weight for the upper layer
    w_z(nz,3) = 1; % weight for the lower layer
  end
end

% ================ TOPO =========================================
fprintf('Topography \n');

filetopo = 'grid/CI-MJM88_bathy_meter.nc';
  
bathy   = double(ncread([dir_1_init,filetopo],'Bathymetry'));
nav_lat = double(ncread([dir_1_init,filetopo],'nav_lat'));
nav_lon = double(ncread([dir_1_init,filetopo],'nav_lon'));

% do some checks
fprintf('check grid: \n');
fprintf('min lat (old, new): %f %f \n', min(min(nav_lat)), min(min(yla)));
fprintf('max lat (old, new): %f %f \n', max(max(nav_lat)), max(max(yla)));

fprintf('min lon (old, new): %f %f \n', min(min(nav_lon)), min(min(xlo)));
fprintf('max lon (old, new): %f %f \n', max(max(nav_lon)), max(max(xlo)));


xlo_t = xlo;
yla_t = yla;
h = griddata(nav_lon,nav_lat,bathy,xlo_t,yla_t); 
h(isnan(h)) = 0;

% sign for mitgcm
h = -h;

% modif topo
%-------- WARNING --------
% may not work at all resolutions 
%-------- WARNING --------


mask_topo = h;
mask_topo(mask_topo ~= 0) = 1;
label = bwlabel(mask_topo);

% pacific
lab_pacific = label(1,1);
mask_topo(label == lab_pacific) = 0;

% mediterranean
%gibraltar
% lat first because the grid is not regular
[argvalue, ila_gib] = min(abs(yy-35.9));
[argvalue, ilo_gib] = min(abs(xx+5.7));
djgib = floor(70.0/dla_km)+1;

% straight line across gibraltar 
mask_topo(ilo_gib,ila_gib-djgib:ila_gib+djgib) = -0;

label = bwlabel(mask_topo);
% hopefully this is a sea point
lab_med = label(ilo_gib+2,ila_gib);
mask_topo(label == lab_med) = 0;

% fill up the rest (North of Gibralar)
%-- MAY BE DANGEROUS --
for ny = ila_gib:si_y_mit
  if mask_topo(si_x_mit,ny) ~= 0
    lab_o = label(si_x_mit,ny);
    mask_topo(label == lab_o) = 0;
  end
end

% northern canada
[argvalue, ila_can] = min(abs(yy-51.0));
[argvalue, ilo_can] = min(abs(xx+71.3));
mask_topo(1:ilo_can,ila_can:si_y_mit) = 0;

% apply mask
h = h.*mask_topo;

mask_mit = 0*h;
if flag_cut
  % gulf of mexico
  for nx = x_c1:x_c2
    for ny = 1:si_y_mit
      if  h(nx,ny) ~=0
        mask_mit(nx,ny) = 1;
        % buffer zone
        if flag_cut == 1
          mask_mit(max(nx-3,1):min(nx+3,x_cut),max(ny-3,1):min(ny+3,si_y_mit)) = 1;
        else
          mask_mit(max(nx-3,x_cut):min(nx+3,si_x_mit),max(ny-3,1):min(ny+3,si_y_mit)) = 1;
        end
      end
    end
  end
end


mask_check = cut_gulf(mask_topo,mask_mit,x_cut,flag_cut,ybc);

h2 = cut_gulf(h,mask_mit,x_cut,flag_cut,ybc);


fid=fopen([dir_o,'topo.box'],'w',ieee); fwrite(fid,h2,accuracy); fclose(fid);

% make a plot with the tile subdivision
if (flag_plot == 1)
  figure;
  contourf(xlo(1:si_x_mit2,:),yla(1:si_x_mit2,:),h2)
  hold on;
  nxp = 10;
  nyp = 9;
  si_xp = floor(si_x_mit2/nxp);
  si_yp = floor(si_y_mit /nyp);
  for nx = 1:nxp
    plot([xlo(1+(nx-1)*si_xp,1),xlo(1+(nx-1)*si_xp,1)],[yla(1+(nx-1)*si_xp,1),yla(1+(nx-1)*si_xp,end)],'r')
  end
  
  for ny = 1:nyp
    plot([xlo(1,1),xlo(si_x_mit2,1)],[yla(1,1+(ny-1)*si_yp),yla(1,1+(ny-1)*si_yp)],'r')
  end
  saveas(gcf,'topo_tiles.png');
  close;
end

fprintf('end topo \n');

% =============== some global data
array_v = {'u';'v';'t';'s';'e'};
def_val = [0.0, 0.0, 20.0, 30.0, 0.0];
array_v2 = {'vozocrtx';
            'vomecrty';
            'votemper';
            'vosaline';
            'sossheig'};


% ================ INITIAL CONDITIONS =========================================
if build_init
  
fprintf('Building Initial conditions \n');

array_f = {'1987/CI-MJM88_y1987m01d05_gridU.nc';
           '1987/CI-MJM88_y1987m01d05_gridV.nc';
           '1987/CI-MJM88_y1987m01d05_gridT.nc';
           '1987/CI-MJM88_y1987m01d05_gridT.nc';
           '1987/CI-MJM88_y1987m01d05_SSHaa.nc'};


for nv = 1:4
  fprintf(array_v{nv});

  fileT = array_f{nv};

  T_i = ncread([dir_1_init,fileT],array_v2{nv});
  nav_lat = double(ncread([dir_1_init,fileT],'nav_lat'));
  nav_lon = double(ncread([dir_1_init,fileT],'nav_lon'));

  [si_x,si_y,si_z] = size(T_i);
  t = zeros(si_x,si_y,si_z_mit);
  
  % Vertical interpolation
  tlast = def_val(nv);
  for nx = 1:si_x
    for ny = 1:si_y
      for nz = 1:si_z_mit
        t(nx,ny,nz) = w_z(nz,3)*T_i(nx,ny,w_z(nz,1)+1) + w_z(nz,2)*T_i(nx,ny,w_z(nz,1));
        if (isnan(t(nx,ny,nz)))
          if (isnan(T_i(nx,ny,w_z(nz,1))))
            t(nx,ny,nz) = tlast;
          else
            t(nx,ny,nz) = T_i(nx,ny,w_z(nz,1));
          end
        else
          tlast = t(nx,ny,nz);
        end
      end
    end
  end
  
  % interpolate data on the NEW grid
  %if ntime == 1; mod2 = 'w'; end
  mod2 = 'w';
  for nz = 1:si_z_mit
    if flag_interp
    t_mit = griddata(nav_lon,nav_lat,squeeze(t(:,:,nz)),xlo_t,yla_t); 
    elseif nz == 1
    fprintf('Creating new interpolant\n')
    [t_mit,tri,wei] = my_griddata1(nav_lon,nav_lat,squeeze(t(:,:,nz)),xlo_t,yla_t,{'QJ'}); %
    t_mit = my_griddata2(nav_lon,nav_lat,squeeze(t(:,:,nz)),xlo_t,yla_t,tri,wei);
    else
    t_mit = my_griddata2(nav_lon,nav_lat,squeeze(t(:,:,nz)),xlo_t,yla_t,tri,wei);
    end
    t_mit(isnan(t_mit)) = def_val(nv);
    t_mit = cut_gulf(t_mit,mask_mit,x_cut,flag_cut,ybc);
    fid=fopen([dir_o, array_v{nv}, '_init.box'],mod2,ieee); fwrite(fid,t_mit,accuracy); fclose(fid);
    mod2 = 'a';
  end
  
  clear t_mit t T_i;
end

% SSH:
nv = 5;
fileT = array_f{nv};
T_i = ncread([dir_1_init,fileT],array_v2{nv});
nav_lat = double(ncread([dir_1_init,fileT],'nav_lat'));
nav_lon = double(ncread([dir_1_init,fileT],'nav_lon'));

T_i(isnan(T_i)) = 0.0;
% eta is on a different grid
t_mit = griddata(nav_lon,nav_lat,squeeze(T_i),xlo_t,yla_t); 

t_mit = cut_gulf(t_mit,mask_mit,x_cut,flag_cut,ybc);
fid=fopen([dir_o, array_v{nv}, '_init.box'],'w',ieee); fwrite(fid,t_mit,accuracy); fclose(fid);

fprintf('\n');
end % build_init

% ================ OBCS =========================================
if build_obcs
  
fprintf('Building OBCS \n');

bc_loc = ['NORTH';
          'SOUTH'];
bc_loc = cellstr(bc_loc);
[si_b, naux] = size(bc_loc);

for nb =1:si_b
  
  fprintf([bc_loc{nb},'\n'])
  
  dir_2 = [bc_loc{nb} 'BDY-MJM88/1987/'];
  
  %gridu = dir([dir_1_obcs, dir_2,'*m01d05*gridU*']);
  %gridv = dir([dir_1_obcs, dir_2,'*m01d05*gridV*']);
  %gridt = dir([dir_1_obcs, dir_2,'*m01d05*gridT*']);

  gridu = dir([dir_1_obcs, dir_2,'*gridU*']);
  gridv = dir([dir_1_obcs, dir_2,'*gridV*']);
  gridt = dir([dir_1_obcs, dir_2,'*gridT*']);

  [si_f, naux] = size(gridu);
  
  for nv = 1:4
    
    ntime1 = 0;
    for nf = 1:si_f
      ntime1 = ntime1 + 1;
      
      if nv == 1
        fileT  = gridu(nf).name;
      elseif nv == 2
        fileT  = gridv(nf).name;
      else
        fileT  = gridt(nf).name;
      end
      
      T_i = ncread([dir_1_obcs,dir_2,fileT],array_v2{nv});
      nav_lat = double(ncread([dir_1_obcs,dir_2,fileT],'nav_lat'));
      nav_lon = double(ncread([dir_1_obcs,dir_2,fileT],'nav_lon'));
      [si_x,si_y,si_z] = size(T_i);

      if (nb == 1)  % North
        xlo_t = xlo(:,si_y_mit);
        yla_t = yla(:,si_y_mit);
      elseif (nb == 1)  % South
        xlo_t = xlo(:,1);
        yla_t = yla(:,1);
      end
      
      t = zeros(si_x,si_y,si_z_mit);
      % Vertical interpolation
      tlast = def_val(nv);
      for nx = 1:si_x
        for ny = 1:si_y
          for nz = 1:si_z_mit
            t(nx,ny,nz) = w_z(nz,3)*T_i(nx,ny,w_z(nz,1)+1) + w_z(nz,2)*T_i(nx,ny,w_z(nz,1));
            if (isnan(t(nx,ny,nz)))
              if (isnan(T_i(nx,ny,w_z(nz,1))))
                t(nx,ny,nz) = tlast;
              else
                t(nx,ny,nz) = T_i(nx,ny,w_z(nz,1));
              end
            else
              tlast = t(nx,ny,nz);
            end
          end
        end
      end
      
      % interpolate data on the NEW grid
      if nf == 1; 
        mod2 = 'w'; 
      else
        mod2 = 'a'; 
      end
      for nz = 1:si_z_mit
        if (nb == 1)  % North
          if flag_interp
            ybc = si_y_mit;
            t_mit = griddata(nav_lon,nav_lat,squeeze(t(:,:,nz)),xlo_t,yla_t); 
          elseif (nf == 1 && nz == 1)
          fprintf('Creating new interpolant\n')
          [t_mit,tri,wei] = my_griddata1(nav_lon,nav_lat,squeeze(t(:,:,nz)),xlo_t,yla_t,{'QJ'}); %
          t_mit = my_griddata2(nav_lon,nav_lat,squeeze(t(:,:,nz)),xlo_t,yla_t,tri,wei);
          else
          t_mit = my_griddata2(nav_lon,nav_lat,squeeze(t(:,:,nz)),xlo_t,yla_t,tri,wei);             
          end
        elseif (nb == 2)  % South
          ybc = 1;
          t_mit = interp1(nav_lon,squeeze(t(:,:,nz)),xlo_t); 
        end
        t_mit(isnan(t_mit)) = def_val(nv);
        %TODO check cut_gulf for BC.
        t_mit = cut_gulf(t_mit,mask_mit,x_cut,flag_cut,ybc);
        fid=fopen([dir_o, array_v{nv}, '_', bc_loc{nb}, '.box'],mod2,ieee); fwrite(fid,t_mit,accuracy); fclose(fid);
        mod2 = 'a'; 
      end

      clear t_mit t T_i;

    end % nv
    
  end % nf
  
end % BC_loc

fprintf('OBCS: %i pts\n',ntime1);
fprintf('period OBCS 86400 : %i\n',86400*5*(ntime1));

end % build OBCS


% % ================ ATMOSPHERIC FILES =========================================

  
if (build_atmosphere)
fprintf('Building atmospheric files \n');

dir_fc = dir_era;

extf = '.nc';

lat_era = ncread([dir_era,'t2.198701', extf],'lat');
lon_era = ncread([dir_era,'t2.198701', extf],'lon');
lat_era = flipdim(lat_era,1);

si_x_era = size(lon_era,1);
si_y_era = size(lat_era,1);
for nx =1:si_x_era
  if (lon_era(nx)> 180); lon_era(nx) = lon_era(nx) - 360; end
end

[y_era,x_era]=meshgrid(lat_era,lon_era);

t_air = zeros(si_x_era,si_y_era);
q_air = zeros(si_x_era,si_y_era);
u_air = zeros(si_x_era,si_y_era);
v_air = zeros(si_x_era,si_y_era);
s_air = zeros(si_x_era,si_y_era);
b_air = zeros(si_x_era,si_y_era);
c_air = zeros(si_x_era,si_y_era);
d_air = zeros(si_x_era,si_y_era);

if flag_interp == 0
fprintf('Creating new interpolant\n')
[t_aux,tri,wei] = my_griddata1(x_era,y_era,t_air,xlo,yla,{'QJ'}); %
end  

%date = {'198701','198702','198703','198704','198705','198706','198707','198708','198709','198710','198711','198712','200701','200702','200703','200704','200705','200706','200707','200708','200709','200710','200711','200712'};
date = {'198701','198702','198703','198704','198705','198706','198707','198708','198709','198710','198711','198712'};
%date = {'199701','199702','199703','199704','199705','199706','199707','199708','199709','199710','199711','199712'};
%date = {'198701'};


ntime2 = 0;
for nfi=1:length(date)
  fprintf([date{nfi},' ']);
time = ncread([dir_era,'t2.',date{nfi},extf],'time');
si_t_era = size(time,1);

for nt = 1:si_t_era
  ntime2 = ntime2 + 1;
  
  t2m = ncread([dir_era,'t2.',date{nfi},extf],'t2',[1 1 nt],[Inf Inf 1]);
  t2m = t2m - 273.16;
  d2m = ncread([dir_era,'d2.',date{nfi},extf],'d2',[1 1 nt],[Inf Inf 1]);
  d2m = d2m - 273.16;
  u10 = ncread([dir_era,'u10.',date{nfi},extf],'u10',[1 1 nt],[Inf Inf 1]);
  v10 = ncread([dir_era,'v10.',date{nfi},extf],'v10',[1 1 nt],[Inf Inf 1]);
  ntsol = nt + floor((nt-0.1)/4);
  sol = ncread([dir_fc,'ssrd.',date{nfi},extf],'ssrd',[1 1 ntsol],[Inf Inf 2]);
  % in W.m-2
  sol = squeeze(sol(:,:,2)-sol(:,:,1))/6./60./60.;
  %  blh = ncread([dir_era,'blh_clim.',date{nfi}(5:6),extf],'blh',[ntsol -1 -1],[ntsol -1 -1]);
  %tcc = ncread([dir_era,'tcc_clim.',date{nfi}(5:6),extf],'tcc',[nt -1 -1],[nt -1 -1]);

  strd = ncread([dir_fc,'strd.',date{nfi},extf],'strd',[1 1 ntsol],[Inf Inf 2]);
  strd = squeeze(strd(:,:,2)-strd(:,:,1))/6./60./60.;
  
  % flipdim
  for nx = 1:si_x_era
    for ny = 1:si_y_era
      t_air(nx,ny) = t2m(nx,si_y_era-ny+1);
      q_air(nx,ny) = d2m(nx,si_y_era-ny+1);
      q_air(nx,ny) = 6.112*exp(17.67*q_air(nx,ny)/(243.5+q_air(nx,ny)))*0.622/1000;
      u_air(nx,ny) = u10(nx,si_y_era-ny+1);
      v_air(nx,ny) = v10(nx,si_y_era-ny+1);
      s_air(nx,ny) = 0.94*sol(nx,si_y_era-ny+1);% 0.94:albedo(ssrd)
      d_air(nx,ny) = strd(nx,si_y_era-ny+1);
    end
  end
  
  % cubic for u and v for cheapaml precip field
  uair_mit = griddata(y_era,x_era,u_air,yla,xlo,'cubic');
  vair_mit = griddata(y_era,x_era,v_air,yla,xlo,'cubic');
  if flag_interp
  tair_mit = griddata(x_era,y_era,t_air,xlo,yla);
  qair_mit = griddata(x_era,y_era,q_air,xlo,yla);
  sair_mit = griddata(x_era,y_era,s_air,xlo,yla);
  dair_mit = griddata(x_era,y_era,d_air,xlo,yla);
  else
  tair_mit = my_griddata2(x_era,y_era,t_air,xlo,yla,tri,wei);
  qair_mit = my_griddata2(x_era,y_era,q_air,xlo,yla,tri,wei);
  sair_mit = my_griddata2(x_era,y_era,s_air,xlo,yla,tri,wei);
  dair_mit = my_griddata2(x_era,y_era,d_air,xlo,yla,tri,wei);
  end
  
  
  if (ntime2 == 1)
    mod1 = 'w';
    
    tair_mit_sauv = tair_mit;
    qair_mit_sauv = qair_mit;
    uair_mit_sauv = uair_mit;
    vair_mit_sauv = vair_mit;
    sair_mit_sauv = sair_mit;
    dair_mit_sauv = dair_mit;

    tair_mit_me = 0.0*tair_mit;
    qair_mit_me = 0.0*qair_mit;
    uair_mit_me = 0.0*uair_mit;
    vair_mit_me = 0.0*vair_mit;
    sair_mit_me = 0.0*sair_mit;

    % ERA-I solar forcast starts at 12h: put 0 for the first two solar
    sair_mit_te = cut_gulf(sair_mit_sauv,mask_mit,x_cut,flag_cut,ybc);
    fid=fopen([dir_o, 'solar.box'],'w',ieee); fwrite(fid,sair_mit_te,accuracy); fclose(fid);
    fid=fopen([dir_o, 'solar.box'],'a',ieee); fwrite(fid,sair_mit_te,accuracy); fclose(fid);

    dair_mit_te = cut_gulf(dair_mit_sauv,mask_mit,x_cut,flag_cut,ybc);
    fid=fopen([dir_o 'longwave.box'],'w',ieee); fwrite(fid,dair_mit_te,accuracy); fclose(fid);
    fid=fopen([dir_o, 'longwave.box'],'a',ieee); fwrite(fid,dair_mit_te,accuracy); fclose(fid);
    
  else
    mod1 = 'a';    
  end
  
  tair_mit_me = ((ntime2-1)*tair_mit_me + tair_mit)/ntime2;
  qair_mit_me = ((ntime2-1)*qair_mit_me + qair_mit)/ntime2;
  uair_mit_me = ((ntime2-1)*uair_mit_me + uair_mit)/ntime2;
  vair_mit_me = ((ntime2-1)*vair_mit_me + vair_mit)/ntime2;
  sair_mit_me = ((ntime2-1)*sair_mit_me + sair_mit)/ntime2;
  
  tair_mit = cut_gulf(tair_mit,mask_mit,x_cut,flag_cut,ybc);
  qair_mit = cut_gulf(qair_mit,mask_mit,x_cut,flag_cut,ybc);
  uair_mit = cut_gulf(uair_mit,mask_mit,x_cut,flag_cut,ybc);
  vair_mit = cut_gulf(vair_mit,mask_mit,x_cut,flag_cut,ybc);
  sair_mit = cut_gulf(sair_mit,mask_mit,x_cut,flag_cut,ybc);
  dair_mit = cut_gulf(dair_mit,mask_mit,x_cut,flag_cut,ybc);

  
  fid=fopen([dir_o,'tair.box'] ,mod1,ieee); fwrite(fid,tair_mit,accuracy); fclose(fid);
  fid=fopen([dir_o,'qair.box'] ,mod1,ieee); fwrite(fid,qair_mit,accuracy); fclose(fid);
  fid=fopen([dir_o,'windx.box'],mod1,ieee); fwrite(fid,uair_mit,accuracy); fclose(fid);
  fid=fopen([dir_o,'windy.box'],mod1,ieee); fwrite(fid,vair_mit,accuracy); fclose(fid);
  % solar already started
  fid=fopen([dir_o,'solar.box'],'a',ieee); fwrite(fid,sair_mit,accuracy); fclose(fid);
  fid=fopen([dir_o,'longwave.box'],'a',ieee); fwrite(fid,dair_mit,accuracy); fclose(fid);

end
end

tair_mit_sauv = cut_gulf(tair_mit_sauv,mask_mit,x_cut,flag_cut,ybc);
qair_mit_sauv = cut_gulf(qair_mit_sauv,mask_mit,x_cut,flag_cut,ybc);
uair_mit_sauv = cut_gulf(uair_mit_sauv,mask_mit,x_cut,flag_cut,ybc);
vair_mit_sauv = cut_gulf(vair_mit_sauv,mask_mit,x_cut,flag_cut,ybc);
sair_mit_sauv = cut_gulf(sair_mit_sauv,mask_mit,x_cut,flag_cut,ybc);
dair_mit_sauv = cut_gulf(dair_mit_sauv,mask_mit,x_cut,flag_cut,ybc);

fid=fopen([dir_o,'tair.box'] ,mod1,ieee); fwrite(fid,tair_mit_sauv,accuracy); fclose(fid);
fid=fopen([dir_o,'qair.box'] ,mod1,ieee); fwrite(fid,qair_mit_sauv,accuracy); fclose(fid);
fid=fopen([dir_o,'windx.box'],mod1,ieee); fwrite(fid,uair_mit_sauv,accuracy); fclose(fid);
fid=fopen([dir_o,'windy.box'],mod1,ieee); fwrite(fid,vair_mit_sauv,accuracy); fclose(fid);
% in principle, no need to add the last solar
fid=fopen([dir_o,'solar.box'],mod1,ieee); fwrite(fid,sair_mit_sauv,accuracy); fclose(fid);
fid=fopen([dir_o,'longwave.box'],mod1,ieee); fwrite(fid,dair_mit_sauv,accuracy); fclose(fid);


fprintf('atm pts: %i !!!\n',ntime2+1);
fprintf('period CAML 21600 : %i\n',86400/4*(ntime2+1));

end

