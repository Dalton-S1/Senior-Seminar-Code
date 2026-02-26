%% XROTOR -> SolidWorks XYZ exporter
% Requires:
%   - 'testrotor.txt' (XROTOR output)
%   - 'airfoil.dat' with columns: x_local (0-1), y_local (-t/2..+t/2)
% Output:
%   - One file per radial station: section_##.txt with X Y Z columns

clear; clc;

%% USER INPUTS
xrotorFile  = 'testrotor.txt';   % your XROTOR output file
airfoilFile = ['2412airfoil_refined.dat'];     % base airfoil coordinates
nPtsChord   = 200;               % points per airfoil (if resampling)
scaleChord  = 1.0;               % additional chord scale factor (usually 1)
rAxisSign   = 1;                 % +1 if radius grows along +X, -1 if -X

% Choose global axes:
%   Z = prop axis, X-Y = rotor plane
% We'll place the quarter-chord on a circle at radius r.

%% LOAD AIRFOIL
afData = load(airfoilFile);                    % [x_local, y_local]
x_af   = afData(:,1);                          % chordwise (0 to 1)
y_af   = afData(:,2);                          % thickness / camber line

% Use raw airfoil coordinates directly
xq = x_af;
yq = y_af;

%% PARSE XROTOR FILE
fid = fopen(xrotorFile,'r');
if fid < 0
    error('Could not open %s', xrotorFile);
end

R = [];    % rotor radius
data = []; % [r/R, c/R, Beta0deg]

while ~feof(fid)
    line = fgetl(fid);
    if ~ischar(line), break; end

    % Get radius (Rad) line
    if contains(line,'Rad') && isempty(R)
        line2 = fgetl(fid);
        vals  = sscanf(line2,'%f');
        R     = vals(1);              % first value is Rad
    end

    % Get radial distribution table start
    if contains(line,'r/R') && contains(line,'C/R') && contains(line,'Beta0deg')
        % Next lines until a line starting with '!' or EOF
        while true
            pos = ftell(fid);
            l2  = fgetl(fid);
            if ~ischar(l2), break; end
            if startsWith(strtrim(l2),'!')
                fseek(fid,pos,'bof');
                break;
            end
            vals = sscanf(l2,'%f');
            if numel(vals) >= 3
                data = [data; vals(1:3).']; %#ok<AGROW>
            end
        end
    end
end
fclose(fid);

if isempty(R) || isempty(data)
    error('Failed to read radius or blade table from XROTOR file.');
end

r_R   = data(:,1);   % non-dimensional radius
c_R   = data(:,2);   % non-dimensional chord
beta0 = data(:,3);   % pitch angle (deg)

r     = r_R * R;                 % [m] radius at each station
c     = c_R * R * scaleChord;    % [m] chord at each station

%% LOOP OVER SECTIONS AND WRITE XYZ FILES
nSec = numel(r);
for i = 1:nSec
    ri    = r(i);
    ci    = c(i);
    beta  = beta0(i);              % deg

    % Local airfoil coordinates, scaled by chord
    x_local = xq * ci;
    y_local = yq * ci;

    % Place quarter-chord at radius ri in rotor plane.
    % Define local coordinates so x_local = 0 at LE, x_local = ci at TE.
    % Quarter-chord point at x_qc = 0.25 * ci.
    x_qc = 0.25 * ci;

    % Local frame: 
    %   s-axis along tangential direction (in rotor plane),
    %   r-axis along radial direction.
    % We rotate the chord about the s-axis by beta.

    % Shift to quarter-chord centering in r-direction
    r_local = (x_local - x_qc);
    s_local = y_local;

    % Convert pitch angle beta about s-axis:
    %   r' =  r_local * cos(beta)
    %   z' =  r_local * sin(beta)
    %   s' =  s_local
    beta_rad = deg2rad(beta);
    r_prime  = r_local * cos(beta_rad);
    z_prime  = r_local * sin(beta_rad);
    s_prime  = s_local;

    % Global rotor-plane coordinates:
    % radial direction = ri + r'
    % tangential direction = s'
    r_global = ri + r_prime;
    s_global = s_prime;

    % Map rotor-plane (r_global, s_global) into 3D (X,Y,Z)
    % Here: radius along X, tangential along Y, prop axis along Z.
    X = rAxisSign * r_global;
    Y = s_global;
    Z = z_prime;

    % Write XYZ file
    outName = sprintf('section_%02d.txt', i);
    outData = [X Y Z];

    fid_out = fopen(outName,'w');
    if fid_out < 0
        error('Could not open %s for writing.', outName);
    end
    fprintf(fid_out,'%.6f %.6f %.6f\n', outData.');
    fclose(fid_out);

    fprintf('Wrote %s (r/R = %.3f, c/R = %.3f, beta = %.2f deg)\n', ...
        outName, r_R(i), c_R(i), beta);
end

disp('Done. Import each section_XX.txt in SolidWorks via Insert -> Curve -> Curve Through XYZ Points.');
