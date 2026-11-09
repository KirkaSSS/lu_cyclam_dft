function T = lu_dft_run(varargin)
%LU_DFT_RUN  Pokrece ceo racun (GFN2-xTB + PBE0-D3 DFT) i vraca tabelu rezultata.
%   T = lu_dft_run('SkipXtb', true, 'Threads', 4)   % samo DFT na geometrijama iz analize (~40-60 min na 8 jezgara)
%   T = lu_dft_run('Threads', 4, 'Jobs', 4)          % ceo tok, ukljucujuci novu pretragu konformera
%   Opcije: Threads (DFT jezgra, 2), Jobs (paralelni xTB poslovi, 2), Memory (MB po DFT racunu, 5500),
%           Seeds (konformera po raspodeli vode, 16), SkipXtb (false), Python (putanja do pythona; podrazumevano venv iz lu_dft_setup).
%   Na Linuxu/macOS se MATLAB-ove biblioteke (LD_LIBRARY_PATH) iskljucuju za Python, da ne bi bilo sukoba.
%   Racun traje dugo; MATLAB je za to vreme zauzet. Izlaz se ispisuje uzivo u Command Window.
p = inputParser;
addParameter(p, 'Threads', 2); addParameter(p, 'Jobs', 2); addParameter(p, 'Memory', 5500);
addParameter(p, 'Seeds', 16); addParameter(p, 'SkipXtb', false); addParameter(p, 'Python', '');
parse(p, varargin{:}); o = p.Results;
root = lu_dft_root();
lu_wsl_check();
args = sprintf('--threads %d --jobs %d --memory %d --seeds %d', o.Threads, o.Jobs, o.Memory, o.Seeds);
if o.SkipXtb, args = [args ' --skip-xtb']; end
if ispc
    py = o.Python; if isempty(py), py = '~/lu-venv/bin/python'; end
    cmd = sprintf('wsl bash -lc "cd ''%s'' && %s run_all.py %s"', lu_wslpath(fullfile(root, 'scripts')), py, args);
else
    py = o.Python; if isempty(py), py = fullfile(root, 'venv', 'bin', 'python'); end
    cmd = sprintf('unset LD_LIBRARY_PATH DYLD_LIBRARY_PATH; cd "%s" && "%s" run_all.py %s', fullfile(root, 'scripts'), py, args);
end
fprintf('%s\n', cmd);
t0 = tic;
status = system(cmd, '-echo');
fprintf('Trajanje: %.1f min\n', toc(t0) / 60);
if status ~= 0
    error('lu_dft:run', 'Racun je prekinut (kod %d). Pogledajte poruke iznad.', status);
end
T = lu_dft_results(fullfile(root, 'work', 'dft.jsonl'));
end
