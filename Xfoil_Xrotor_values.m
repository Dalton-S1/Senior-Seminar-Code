function A = parse_xfoil_to_xrotor()

filename = 'C:\XFOIL6.99\3117.txt';

txt = fileread(filename);

% X-foil Header (Re, Mach)
A = struct();
A.Re = NaN;
A.Mach = NaN;

tok = regexp(txt,'Re\s*=\s*([0-9\.Ee\+\- ]+)','tokens','once');
if ~isempty(tok)
    s = regexprep(tok{1},'\s*e\s*','e');
    A.Re = str2double(s);
end

tok = regexp(txt,'Mach\s*=\s*([0-9\.Ee\+\- ]+)','tokens','once');
if ~isempty(tok)
    A.Mach = str2double(tok{1});
end

% Read table
table_start = regexp(txt,'\n\s*alpha\s+CL\s+CD\s+CDp\s+CM','once');
if isempty(table_start), error("Table header not found."); end

table_txt = txt(table_start:end);
lines = regexp(table_txt,'\n\s*([-0-9\.]+[^\n]*)','tokens');

alpha = []; CL = []; CD = []; CM = [];

for i = 1:numel(lines)
    row = strtrim(lines{i}{1});
    nums = regexp(row,'([-+]?\d*\.?\d+(?:[eE][\+\-]?\d+)?)','match');
    if numel(nums) >= 5
        alpha(end+1,1) = str2double(nums{1});
        CL(end+1,1)    = str2double(nums{2});
        CD(end+1,1)    = str2double(nums{3});
        CM(end+1,1)    = str2double(nums{5});
    end
end

if numel(alpha) < 3
    error("Too few polar rows parsed (%d).", numel(alpha));
end



% "-0"
alpha(abs(alpha) < 1e-12) = 0;

% Remove non-finite rows
good = isfinite(alpha) & isfinite(CL) & isfinite(CD) & isfinite(CM);
alpha = alpha(good); CL = CL(good); CD = CD(good); CM = CM(good);

% Sort by alpha
[alpha, idx] = sort(alpha);
CL = CL(idx); CD = CD(idx); CM = CM(idx);

% Remove duplicate alpha rows 
[alpha_u, ~, g] = unique(alpha,'stable');
if numel(alpha_u) < numel(alpha)
    CL = accumarray(g, CL, [], @mean);
    CD = accumarray(g, CD, [], @mean);
    CM = accumarray(g, CM, [], @mean);
    alpha = alpha_u;
end

%Step grid
alpha_grid = (-6:0.1:16)';

CL = interp1(alpha, CL, alpha_grid, 'linear','extrap');
CD = interp1(alpha, CD, alpha_grid, 'linear','extrap');
CM = interp1(alpha, CM, alpha_grid, 'linear','extrap');

alpha = alpha_grid;


% Zero-lift angle

good2 = isfinite(CL) & isfinite(alpha);
CL2 = CL(good2);
alpha2 = alpha(good2);

[CLu, ~, gg] = unique(CL2, 'stable');
if numel(CLu) < numel(CL2)
    alpha_u2 = accumarray(gg, alpha2, [], @mean);
else
    alpha_u2 = alpha2;
end

A.zero_lift_alpha_deg = interp1(CLu, alpha_u2, 0, 'linear','extrap');


% dCL/dalpha (linear region)

mask = (alpha > 2.5) & (alpha < 5);
if nnz(mask) < 2
    mask = (alpha >= 0) & (alpha <= 6);
end
p = polyfit(alpha(mask), CL(mask), 1);
A.dCl_dalpha_per_rad = p(1) * 180/pi;

% dCL/dalpha near stall
A.dCl_dalpha_stall_per_rad = 0.10;

% Max CL
[A.Cl_max, idxMax] = max(CL);
A.alpha_at_Cl_max = alpha(idxMax);

% Min CL
[A.Cl_min, idxMin] = min(CL);
A.alpha_at_Cl_min = alpha(idxMin);

% Increment from linear prediction to stall 
A.deltaCl_last_linear_to_stall = 0.10;

% Cd(0), Cl(0), Cm(0)
i0 = find(alpha == 0, 1);
A.Cd_alpha0 = CD(i0);
A.Cl_alpha0 = CL(i0);
A.Cm_alpha0 = CM(i0);

% Slope of Cd vs Cl^2
Cl2 = CL.^2;
k = polyfit(Cl2, CD, 1);
A.dCd_dCl2 = k(1);

% Re-scaling exponent
if ~isnan(A.Re) && A.Re < 2e5
    A.Re_scaling_exponent = -0.8;
else
    A.Re_scaling_exponent = -0.2;
end


% SAVE/APPEND RESULTS TO CSV FILE
outputFile = 'airfoil_database.csv';

[~, airfoilName, ~] = fileparts(filename);
A.airfoil_name = airfoilName;

T = struct2table(A);

if isfile(outputFile)
    writetable(T, outputFile, 'WriteMode','append');
else
    writetable(T, outputFile);
end

fprintf('Saved results to %s (airfoil: %s)\n', outputFile, airfoilName);
disp(A)

end
