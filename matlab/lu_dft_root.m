function root = lu_dft_root()
%LU_DFT_ROOT  Koreni folder paketa (folder iznad matlab/).
root = fileparts(fileparts(mfilename('fullpath')));
end
