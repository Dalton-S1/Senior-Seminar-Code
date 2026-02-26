function Tsorted = rank_NACA_from_polars_only()
% Rank NACA*.txt XFOIL polars:
%   Cd_alpha0, dCd_dCl2, Cl_max, dCl_dalpha_per_rad, Cm_alpha0
%
% 
%
% Incomplete:
%   - <16 numeric rows OR missing header => Incomplete
%
% Ranking:
%   Cl_star = 0.70
%   Score = Cl_star / (Cd_alpha0 + dCd_dCl2 * Cl_star^2)
%
% Outputs:
%   NACA_ranked.csv
%   NACA_ranked.xlsx


folderPath = 'C:\XFOIL6.99\polars\Airfoil Data';
pattern    = 'NACA*.txt';

outCSV  = fullfile(folderPath, 'NACA_ranked.csv');
outXLSX = fullfile(folderPath, 'NACA_ranked.xlsx');

Cl_star = 0.70;            
min_rows_required = 16;    % <-- incomplete threshold

% Optional constraints
min_Clmax   = -Inf;  % e.g. 1.20
max_abs_Cm0 = Inf;   % e.g. 0.20

alpha_zero_tol = 1e-12;    % fixes "-0"

files = dir(fullfile(folderPath, pattern));
if isempty(files)
    error('No files found: %s', fullfile(folderPath, pattern));
end

Aall = [];

for k = 1:numel(files)
    fname = fullfile(files(k).folder, files(k).name);
    [~, airfoilName, ~] = fileparts(fname);

    % Read only polar table
    [alpha, CL, CD, CM, nRows, hasHeader] = read_xfoil_polar_table(fname);

    % Snap "-0"
    if ~isempty(alpha)
        alpha(abs(alpha) < alpha_zero_tol) = 0;
    end

    % Base record
    A = struct();
    A.airfoil_name  = airfoilName;
    A.source_file   = fname;
    A.raw_row_count = nRows;
    A.Incomplete    = (~hasHeader) || (nRows < min_rows_required);

    % Metrics
    A.Cd_alpha0 = NaN;
    A.dCd_dCl2  = NaN;
    A.Cl_max    = NaN;
    A.alpha_at_Cl_max = NaN;
    A.dCl_dalpha_per_rad = NaN;
    A.Cm_alpha0 = NaN;

    if A.Incomplete
        fprintf('INCOMPLETE: %s (rows=%d)\n', files(k).name, nRows);
        Aall = [Aall; A]; 
        continue
    end

    % Sort by alpha
    [alpha, idx] = sort(alpha);
    CL = CL(idx); CD = CD(idx); CM = CM(idx);

    % Snap again after sorting
    alpha(abs(alpha) < alpha_zero_tol) = 0;

    % Remove duplicate alpha rows 
    [alpha_u, ~, g] = unique(alpha, 'stable');
    if numel(alpha_u) < numel(alpha)
        CL = accumarray(g, CL, [], @mean);
        CD = accumarray(g, CD, [], @mean);
        CM = accumarray(g, CM, [], @mean);
        alpha = alpha_u;
    end

    % If duplicates collapse, mark incomplete
    if numel(alpha) < min_rows_required
        A.Incomplete = true;
        A.raw_row_count = numel(alpha);
        fprintf('INCOMPLETE after de-dup: %s (unique rows=%d)\n', files(k).name, numel(alpha));
        Aall = [Aall; A]; 
        continue
    end

    % ---- Metric: Cl_max ----
    [A.Cl_max, iMax] = max(CL);
    A.alpha_at_Cl_max = alpha(iMax);

    % ---- Metrics: Cd_alpha0, Cm_alpha0 at alpha=0 (prefer exact 0) ----
    i0 = find(alpha == 0, 1, 'first');
    if ~isempty(i0)
        A.Cd_alpha0 = CD(i0);
        A.Cm_alpha0 = CM(i0);
    else
        A.Cd_alpha0 = interp1(alpha, CD, 0, 'linear', NaN);
        A.Cm_alpha0 = interp1(alpha, CM, 0, 'linear', NaN);
    end

    % ---- Metric: dCd_dCl2 from Cd vs Cl^2 ----
    good = isfinite(CL) & isfinite(CD);
    if nnz(good) >= 3
        Cl2 = (CL(good)).^2;
        p = polyfit(Cl2, CD(good), 1);  
        A.dCd_dCl2 = p(1);
    end

    % ---- Metric: dCl/dalpha near low alpha ----
    mask = (alpha >= 0) & (alpha <= 5) & isfinite(CL);
    if nnz(mask) < 2
        mask = (alpha >= 1) & (alpha <= 6) & isfinite(CL);
    end
    if nnz(mask) >= 2
        p = polyfit(alpha(mask), CL(mask), 1);     % per-deg slope = p(1)
        A.dCl_dalpha_per_rad = p(1) * (180/pi);    % per-deg -> per-rad
    end

    fprintf('OK: %s (rows=%d, unique=%d)\n', files(k).name, nRows, numel(alpha));
    Aall = [Aall; A]; %#ok<AGROW>
end

T = struct2table(Aall);

% ---------- Ranking ----------
T.Cl_star  = repmat(Cl_star, height(T), 1);
T.Cd_at_Cl = T.Cd_alpha0 + T.dCd_dCl2 .* (Cl_star.^2);
T.Score_LD = Cl_star ./ T.Cd_at_Cl;

% Incomplete or NaNs -> bottom
badScore = T.Incomplete | ~isfinite(T.Score_LD) | isnan(T.Score_LD);
T.Score_LD(badScore) = -Inf;

% Optional constraints (only for complete)
T.Pass = true(height(T),1);
idx = ~T.Incomplete;

if ~isinf(min_Clmax)
    T.Pass(idx) = T.Pass(idx) & (T.Cl_max(idx) >= min_Clmax);
end
if ~isinf(max_abs_Cm0)
    T.Pass(idx) = T.Pass(idx) & (abs(T.Cm_alpha0(idx)) <= max_abs_Cm0);
end

T.Score_LD(~T.Pass & ~T.Incomplete) = -Inf;

% Sort and rank
Tsorted = sortrows(T, 'Score_LD', 'descend');
Tsorted.Rank = (1:height(Tsorted))';

% Save outputs
writetable(Tsorted, outCSV);
writetable(Tsorted, outXLSX);

fprintf('\nSaved:\n  %s\n  %s\n', outCSV, outXLSX);

% Read rows 
    function [alpha, CL, CD, CM, nRows, hasHeader] = read_xfoil_polar_table(file)
        txt = fileread(file);

        table_start = regexp(txt,'\n\s*alpha\s+CL\s+CD\s+CDp\s+CM','once');
        hasHeader = ~isempty(table_start);

        alpha = []; CL = []; CD = []; CM = [];
        nRows = 0;

        if ~hasHeader
            return
        end

        table_txt = txt(table_start:end);
        lines = regexp(table_txt,'\n\s*([-0-9\.]+[^\n]*)','tokens');

        for ii = 1:numel(lines)
            row = strtrim(lines{ii}{1});
            nums = regexp(row,'([-+]?\d*\.?\d+(?:[eE][\+\-]?\d+)?)','match');
            if numel(nums) >= 5
                nRows = nRows + 1;
                alpha(end+1,1) = str2double(nums{1});
                CL(end+1,1)    = str2double(nums{2});
                CD(end+1,1)    = str2double(nums{3});
                CM(end+1,1)    = str2double(nums{5});
            end
        end
    end


end
