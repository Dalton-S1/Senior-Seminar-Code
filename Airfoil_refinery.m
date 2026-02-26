%% Airfoil densifier: airfoil.dat -> airfoil_refined.dat
% Reads an airfoil coordinate file (plain XFOIL format):
%   Line 1: name / comment
%   Lines 2..N: x y
% Increases point count via interpolation and writes a new .dat file.

clear; clc;

%% USER SETTINGS
inFile      = '2412airfoil.dat';          % original airfoil file
outFile     = '2412airfoil_refined.dat';  % output file
nPtsTarget  = 400;                    % desired total number of points

%% READ ORIGINAL AIRFOIL FILE
fid = fopen(inFile,'r');
if fid < 0
    error('Could not open %s', inFile);
end

% First line is name/comment
nameLine = strtrim(fgetl(fid));

% Remaining lines: coordinates
data = textscan(fid,'%f %f','CollectOutput',true);
fclose(fid);

coords = data{1};            % [x y] size N x 2
x      = coords(:,1);
y      = coords(:,2);
N      = numel(x);

if N < 5
    error('Airfoil file has too few points (%d).', N);
end

%% PARAMETRIZE BY CUMULATIVE ARC LENGTH
% Arc-length parameter s from 0 to 1, following existing point order
dx = diff(x);
dy = diff(y);
ds = sqrt(dx.^2 + dy.^2);
s  = [0; cumsum(ds)];
s  = s / s(end);              % normalize 0..1

% Query parameter values for refined distribution
sq = linspace(0,1,nPtsTarget).';

%% INTERPOLATE X AND Y ALONG ARC LENGTH
xq = interp1(s, x, sq, 'pchip');
yq = interp1(s, y, sq, 'pchip');

%% WRITE REFINED AIRFOIL FILE
fid = fopen(outFile,'w');
if fid < 0
    error('Could not open %s for writing.', outFile);
end

% Write name/comment line
fprintf(fid,'%s (refined to %d pts)\n', nameLine, numel(xq));

% Write coordinates
for i = 1:numel(xq)
    fprintf(fid,'%.6f %.6f\n', xq(i), yq(i));
end
fclose(fid);

fprintf('Wrote %s with %d points (input had %d points).\n', outFile, numel(xq), N);
