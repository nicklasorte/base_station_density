function [lat1D, lon1D, id] = randCircumference_scircle1(lat0, lon0, radius_km, N, npts)
% randCircumference_scircle1  Random points on circumference using scircle1
% (Output reshaped to 1-D column vectors)
%
%   [lat1D, lon1D, id] = randCircumference_scircle1(lat0, lon0, radius_km, N)
%
% Inputs:
%   lat0, lon0   - Mx1 centers (degrees)
%   radius_km    - Mx1 radii (km)
%   N            - number of random points per circle
%   npts         - (optional) vertices used to construct circle (default 2000)
%
% Outputs:
%   lat1D        - (M*N)x1 latitude vector (deg)
%   lon1D        - (M*N)x1 longitude vector (deg)
%   id           - (M*N)x1 circle ID (1..M)
%
% Requires Mapping Toolbox (scircle1).

if nargin < 5 || isempty(npts)
    npts = 2000;
end

lat0 = lat0(:);
lon0 = lon0(:);
radius_km = radius_km(:);

M = numel(lat0);
assert(numel(lon0)==M && numel(radius_km)==M, ...
    'lat0, lon0, radius_km must have the same length.');

% WGS84 ellipsoid in km
wgs84 = wgs84Ellipsoid("km");

% Generate full small circles
[latC, lonC] = scircle1(lat0, lon0, radius_km, [], ...
                        wgs84, "degrees", npts);

% Random sampling indices
idx = randi(npts, M, N);
col = repmat(1:M, N, 1).';
lin = sub2ind([npts, M], idx, col);

% Extract sampled points
latSample = reshape(latC(lin), M, N);
lonSample = reshape(lonC(lin), M, N);

% Wrap longitude
lonSample = mod(lonSample + 180, 360) - 180;

% ---- Reshape to 1-D column vectors ----
lat1D = latSample(:);
lon1D = lonSample(:);
id    = repelem((1:M).', N);

end

% % function [latP, lonP, id, stacked] = randCircumference_scircle1(lat0, lon0, radius_km, N, npts)
% % % randCircumference_scircle1  Random points on circumference using scircle1
% % %
% % %   [latP, lonP] = randCircumference_scircle1(lat0, lon0, radius_km, N)
% % %   returns MxN arrays of random points along each small-circle perimeter,
% % %   where M = numel(lat0). Uses WGS84 ellipsoid in kilometers.
% % %
% % % Inputs:
% % %   lat0, lon0    - Mx1 (or 1xM) centers (degrees)
% % %   radius_km     - Mx1 radii (km)
% % %   N             - number of random points per circle (scalar)
% % %   npts          - (optional) number of vertices used to represent circle (default 2000)
% % %
% % % Outputs:
% % %   latP, lonP    - MxN random perimeter points (degrees)
% % %   id            - MxN circle id labels (1..M)
% % %   stacked       - (M*N)x3 [lat lon id]
% % %
% % % Notes:
% % %   - scircle1(...,[],ellipsoid,units,npts) interprets r as linear distance
% % %     in ellipsoid units (here km). :contentReference[oaicite:1]{index=1}
% % %   - Randomness is uniform in the circle parameterization (index along vertices),
% % %     which is what you usually want for “random around the circumference”.
% % 
% % if nargin < 5 || isempty(npts)
% %     npts = 2000;  % dense enough for sampling; increase if you need smoother
% % end
% % 
% % % Columnize
% % lat0 = lat0(:);
% % lon0 = lon0(:);
% % radius_km = radius_km(:);
% % 
% % M = numel(lat0);
% % assert(numel(lon0)==M && numel(radius_km)==M, ...
% %     'lat0, lon0, radius_km must have the same length.');
% % 
% % % WGS84 ellipsoid with km units (Mapping Toolbox)
% % wgs84 = wgs84Ellipsoid("km");
% % 
% % % Full circle: set az = [] (per scircle1 docs) :contentReference[oaicite:2]{index=2}
% % az = [];
% % 
% % % scircle1 can return circles for multiple centers/radii:
% % % latC, lonC are npts-by-M, one column per circle. :contentReference[oaicite:3]{index=3}
% % [latC, lonC] = scircle1(lat0, lon0, radius_km, az, wgs84, "degrees", npts);
% % 
% % % Sample random vertices along each column
% % idx = randi(npts, M, N);                   % MxN indices into rows
% % col = repmat(1:M, N, 1).';                 % MxN column indices
% % 
% % lin = sub2ind([npts, M], idx, col);
% % 
% % latP = reshape(latC(lin), M, N);
% % lonP = reshape(lonC(lin), M, N);
% % 
% % % Wrap lon to [-180,180)
% % lonP = mod(lonP + 180, 360) - 180;
% % 
% % % Labels
% % id = repmat((1:M).', 1, N);
% % stacked = [latP(:), lonP(:), id(:)];
% % end