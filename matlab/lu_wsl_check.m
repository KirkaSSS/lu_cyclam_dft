function lu_wsl_check()
%LU_WSL_CHECK  Windows: proverava WSL i Python okruzenje pre racuna i daje citljivu poruku o gresci.
if ~ispc, return; end
setenv('WSL_UTF8', '1');   % bez ovoga MATLAB od WSL poruka prikaze samo prvo slovo (UTF-16)
[s, out] = system('wsl -l -q');
out = strrep(out, char(0), '');
if s ~= 0 || isempty(strtrim(out))
    error('lu_dft:wsl', ['WSL nije instaliran ili nema Linux distribuciju.\n' ...
        'U PowerShell-u (kao administrator) pokrenite:  wsl --install -d Ubuntu\n' ...
        'pa restartujte racunar, otvorite Ubuntu jednom (korisnik i lozinka) i u njemu:\n' ...
        '  sudo apt update && sudo apt install -y python3-venv python3-pip\n' ...
        'Poruka WSL-a: %s'], strtrim(out));
end
[s, out] = system('wsl bash -lc "test -x ~/lu-venv/bin/python && echo ok"');
if s ~= 0 || ~contains(out, 'ok')
    error('lu_dft:venv', ['Python okruzenje u WSL-u ne postoji. Pokrenite prvo:  lu_dft_setup\n' ...
        '(WSL distribucije: %s)'], strtrim(strrep(evalc('system(''wsl -l -q'');'), char(0), '')));
end
end
