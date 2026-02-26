function [nRows, hasHeader] = count_xfoil_polar_rows(filename)

txt = fileread(filename);

table_start = regexp(txt,'\n\s*alpha\s+CL\s+CD\s+CDp\s+CM','once');
hasHeader = ~isempty(table_start);

if ~hasHeader
    nRows = 0;
    return
end

table_txt = txt(table_start:end);
lines = regexp(table_txt,'\n\s*([-0-9\.]+[^\n]*)','tokens');

nRows = 0;
for i = 1:numel(lines)
    row = strtrim(lines{i}{1});
    nums = regexp(row,'([-+]?\d*\.?\d+(?:[eE][\+\-]?\d+)?)','match');
    if numel(nums) >= 5
        nRows = nRows + 1;
    end
end

end
