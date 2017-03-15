%Domaine initial chaocean
lat_min=-20;
lat_max=54.9;
lon_min=-98;
lon_max=14;
% longitude de coupure pour la reduction du domaine
lon_cut=-11.35;
%bathymrety (je dois changer ça en utilisant les fichiers Grenoble)
[elev1,lon1,lat1]=m_tbase([lon_min lon_max lat_min lat_max]);
% puting all elevation on land to zero
elev1(find(elev1>0))=0;
% puting Mediterranean sea to zero
[lon1,lat1,elev1]=fill_mediterranean(lon1,lat1,elev1);
% puting Pacific Ocean to zero
[lon1,lat1,elev1]=fill_pacific(lon1,lat1,elev1);
% puting Great Lakes to zero
[lon1,lat1,elev1]=fill_greatlakes(lon1,lat1,elev1);
%index coupure
[dummy,closest_lon]=min(abs(lon1(1,:)-lon_cut));
index_lon_cut=closest_lon-1; clear dummy closest_lon
% Making the new grid
elev1_save=elev1;
elev1_save(:,1:end-index_lon_cut+1)=elev1(:,index_lon_cut:end);
mask=elev1 | elev1_save; mask2=elev1 & elev1_save;
elev_tempo=(elev1.*mask + elev1_save.*mask)./(mask2+1);
% Domaine final
lon2=lon1(:,1:index_lon_cut);
lat2=lat1(:,1:index_lon_cut);
elev2=elev_tempo(:,1:index_lon_cut);