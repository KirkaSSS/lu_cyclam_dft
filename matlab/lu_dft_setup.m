function lu_dft_setup(varargin)
%LU_DFT_SETUP  Jednom: napravi Python okruzenje i instalira pakete (ase, rdkit, tblite, pyscf, dftd3, bse).
%   lu_dft_setup                      % Linux/macOS: venv u <paket>/venv
%   lu_dft_setup('Python','python3.12')
%   Na Windowsu sve ide kroz WSL2: okruzenje se pravi u ~/lu-venv unutar Ubuntu-a.
p = inputParser;
addParameter(p, 'Python', 'python3');
parse(p, varargin{:}); o = p.Results;
root = lu_dft_root();
if ispc
    setenv('WSL_UTF8', '1');
    [s0, out0] = system('wsl -l -q'); out0 = strrep(out0, char(0), '');
    if s0 ~= 0 || isempty(strtrim(out0))
        error('lu_dft:wsl', 'WSL nije instaliran. U PowerShell-u (administrator): wsl --install -d Ubuntu, restart, pa u Ubuntu-u: sudo apt update && sudo apt install -y python3-venv python3-pip. Poruka: %s', strtrim(out0));
    end
    req = lu_wslpath(fullfile(root, 'requirements.txt'));
    cmd = sprintf(['wsl bash -lc "python3 -m venv ~/lu-venv && ~/lu-venv/bin/pip install -U pip ' ...
                   '&& ~/lu-venv/bin/pip install -r ''%s''"'], req);
else
    venv = fullfile(root, 'venv');
    cmd = sprintf('unset LD_LIBRARY_PATH DYLD_LIBRARY_PATH; "%s" -m venv "%s" && "%s" install -U pip && "%s" install -r "%s"', ...
        o.Python, venv, fullfile(venv, 'bin', 'pip'), fullfile(venv, 'bin', 'pip'), fullfile(root, 'requirements.txt'));
end
fprintf('%s\n', cmd);
status = system(cmd, '-echo');
if status ~= 0
    error('lu_dft:setup', 'Instalacija nije uspela (kod %d). Pogledajte poruke iznad.', status);
end
fprintf('Okruzenje je spremno.\n');
end
