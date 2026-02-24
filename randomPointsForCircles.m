function pts = randomPointsForCircles(lat0, lon0, radius_km, N)
% randomPointsForCircles  Generate N random points within multiple
%   lat/lon circles on Earth.
%
%   pts = randomPointsForCircles(lat0, lon0, radius_km, N)
%
% Inputs:
%   lat0        - Mx1 vector of center latitudes (deg)
%   lon0        - Mx1 vector of center longitudes (deg)
%   radius_km   - Mx1 vector of circle radii (km)
%   N           - number of random points per circle
%
% Output:
%   pts         - Mx1 cell array; each pts{i} is Nx2 [lat lon] (deg)

% Validate vector sizes
assert(isvector(lat0) && isvector(lon0) && isvector(radius_km), ...
    'lat0, lon0, radius_km must be vectors of equal length.');
assert(numel(lat0)==numel(lon0) && numel(lat0)==numel(radius_km), ...
    'lat0, lon0, and radius_km must be same length.');

M = numel(lat0);
pts = cell(M,1);

% Earth radius (km)
R = 6371;

% Loop through each circle
for i = 1:M
    % center in radians
    lat_c = deg2rad(lat0(i));
    lon_c = deg2rad(lon0(i));
    r = radius_km(i);

    % random radial distances and bearings
    u = rand(N,1);
    v = rand(N,1);
    r_km = r * sqrt(u);            % uniform area
    delta = r_km / R;              % angular distance
    theta = 2*pi*v;

    % spherical computations
    lat = asin( sin(lat_c)*cos(delta) + cos(lat_c).*sin(delta).*cos(theta) );
    lon = lon_c + atan2( sin(theta).*sin(delta).*cos(lat_c), ...
                         cos(delta) - sin(lat_c).*sin(lat) );

    % store in degrees
    pts{i} = [rad2deg(lat), rad2deg(lon)];
end

end