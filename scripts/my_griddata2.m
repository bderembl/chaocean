function [xi] = my_griddata2(x,y,z,xi,yi,tri,w)
%GRIDDATA Data gridding and surface fitting.
%   ZI = GRIDDATA(X,Y,Z,XI,YI) fits a surface of the form Z = F(X,Y) to the
%   data in the (usually) nonuniformly-spaced vectors (X,Y,Z). GRIDDATA
%   interpolates this surface at the points specified by (XI,YI) to produce
%   ZI.  The surface always goes through the data points. XI and YI are
%   usually a uniform grid (as produced by MESHGRID) and is where GRIDDATA
%   gets its name.
%
%   XI can be a row vector, in which case it specifies a matrix with
%   constant columns. Similarly, YI can be a column vector and it specifies
%   a matrix with constant rows.
%
%   [XI,YI,ZI] = GRIDDATA(X,Y,Z,XI,YI) also returns the XI and YI formed
%   this way (the results of [XI,YI] = MESHGRID(XI,YI)).
%
%   [...] = GRIDDATA(X,Y,Z,XI,YI,METHOD) where METHOD is one of
%       'linear'    - Triangle-based linear interpolation (default)
%       'cubic'     - Triangle-based cubic interpolation
%       'nearest'   - Nearest neighbor interpolation
%       'v4'        - MATLAB 4 griddata method
%   defines the type of surface fit to the data. The 'cubic' and 'v4'
%   methods produce smooth surfaces while 'linear' and 'nearest' have
%   discontinuities in the first and zero-th derivative respectively.  All
%   the methods except 'v4' are based on a Delaunay triangulation of the
%   data.
%   If METHOD is [], then the default 'linear' method will be used.
%
%   [...] = GRIDDATA(X,Y,Z,XI,YI,METHOD,OPTIONS) specifies a cell array 
%   of strings OPTIONS to be used as options in Qhull via DELAUNAYN. 
%   If OPTIONS is [], the default DELAUNAYN options will be used.
%   If OPTIONS is {''}, no options will be used, not even the default.
%
%   Example:
%      rand('seed',0)
%      x = rand(100,1)*4-2; y = rand(100,1)*4-2; z = x.*exp(-x.^2-y.^2);
%      ti = -2:.25:2; 
%      [xi,yi] = meshgrid(ti,ti);
%      zi = griddata(x,y,z,xi,yi);
%      mesh(xi,yi,zi), hold on, plot3(x,y,z,'o'), hold off
%
%   See also GRIDDATA3, GRIDDATAN, DELAUNAY, INTERP2, MESHGRID, DELAUNAYN.

%   Copyright 1984-2006 The MathWorks, Inc. 
%   $Revision: 5.33.4.6 $  $Date: 2006/06/20 20:11:50 $

%error(nargchk(5,7,nargin))
narginchk(5,7)

[msg,x,y,z,xi,yi] = xyzchk(x,y,z,xi,yi);
if ~isempty(msg), error(msg); end

if ( issparse(x) || issparse(y) || issparse(z) || issparse(xi) || issparse(yi) )
    error('MATLAB:griddata:InvalidDataSparse',...
        'Input data cannot be sparse.');
end

if ( ~isreal(x) || ~isreal(y) || ~isreal(xi) || ~isreal(yi) )
    error('MATLAB:griddata:InvalidDataComplex',...
        'Input data cannot be complex.');
end

method = 'linear';
opt = [];

% Sort x and y so duplicate points can be averaged before passing to delaunay

%Need x,y and z to be column vectors
sz = numel(x);
x = reshape(x,sz,1);
y = reshape(y,sz,1);
z = reshape(z,sz,1);
sxyz = sortrows([x y z],[2 1]);
x = sxyz(:,1);
y = sxyz(:,2);
z = sxyz(:,3);
myepsx = eps(0.5 * (max(x) - min(x)))^(1/3);
myepsy = eps(0.5 * (max(y) - min(y)))^(1/3);
ind = [0; ((abs(diff(y)) < myepsy) & (abs(diff(x)) < myepsx)); 0];

if sum(ind) > 0
  warning('MATLAB:griddata:DuplicateDataPoints',['Duplicate x-y data points ' ...
            'detected: using average of the z values.']);
  fs = find(ind(1:end-1) == 0 & ind(2:end) == 1);
  fe = find(ind(1:end-1) == 1 & ind(2:end) == 0);
  for i = 1 : length(fs)
    % averaging z values
    z(fe(i)) = mean(z(fs(i):fe(i)));
  end
  x = x(~ind(2:end));
  y = y(~ind(2:end));
  z = z(~ind(2:end));
end

    zi = linear(x,y,z,xi,yi,opt,tri,w);
  
if nargout<=1, xi = zi; end


%------------------------------------------------------------
function zi = linear(x,y,z,xi,yi,opt,tri,w)
%LINEAR Triangle-based linear interpolation

%   Reference: David F. Watson, "Contouring: A guide
%   to the analysis and display of spacial data", Pergamon, 1994.

siz = size(xi);
xi = xi(:); yi = yi(:); % Treat these as columns
x = x(:); y = y(:); % Treat these as columns


z = z(:).'; % Treat z as a row so that code below involving
            % z(tri) works even when tri is 1-by-3.
zi = sum(z(tri) .* w,2);

zi = reshape(zi,siz);

%------------------------------------------------------------
