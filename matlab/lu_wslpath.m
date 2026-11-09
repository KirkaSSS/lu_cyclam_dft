function p = lu_wslpath(winPath)
%LU_WSLPATH  Pretvara Windows putanju (C:\Users\...) u WSL putanju (/mnt/c/Users/...).
p = strrep(char(winPath), '\', '/');
if numel(p) >= 2 && p(2) == ':'
    p = ['/mnt/' lower(p(1)) p(3:end)];
end
end
