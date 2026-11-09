function T = lu_dft_results(file, doPlot)
%LU_DFT_RESULTS  Ucitava dft.jsonl i vraca tabelu relativnih energija (kcal/mol), naboja i dipola.
%   T = lu_dft_results                                  % originalni rezultati iz analize
%   T = lu_dft_results(fullfile(lu_dft_root,'work','dft.jsonl'))   % vasi novi rezultati
%   T = lu_dft_results(file, false)                     % bez grafika
if nargin < 1 || isempty(file), file = fullfile(lu_dft_root(), 'results', 'dft', 'dft.jsonl'); end
if nargin < 2, doPlot = true; end
lines = strtrim(strsplit(fileread(file), {'\r\n', '\n'}));
lines = lines(~cellfun('isempty', lines));
n = numel(lines); H2KCAL = 627.509;
Struktura = strings(n, 1); Sistem = strings(n, 1);
E = zeros(n, 1); Ed3 = zeros(n, 1); q = zeros(n, 1); muCoC = zeros(n, 1); muLu = zeros(n, 1); Konv = false(n, 1);
for i = 1:n
    r = jsondecode(lines{i});
    [~, nm] = fileparts(r.f); Struktura(i) = string(nm);
    if r.nao < 400, Sistem(i) = "Lu-akva kalibracija"; else, Sistem(i) = "Lu-ciklam kompleks"; end
    E(i) = r.E; Ed3(i) = r.Ed3; q(i) = r.qLuMull; muCoC(i) = r.dipCoC; muLu(i) = r.dipLu; Konv(i) = r.conv;
end
dE = nan(n, 1); dE_bezD3 = nan(n, 1);
for s = unique(Sistem)'
    k = Sistem == s;
    dE(k) = (E(k) + Ed3(k) - min(E(k) + Ed3(k))) * H2KCAL;
    dE_bezD3(k) = (E(k) - min(E(k))) * H2KCAL;
end
T = table(Sistem, Struktura, dE, dE_bezD3, q, muCoC, muLu, Konv, ...
    'VariableNames', {'Sistem', 'Struktura', 'dE_kcal', 'dE_bezD3_kcal', 'qLu_Mulliken', 'mu_CoC_D', 'mu_Lu_D', 'Konvergirao'});
T = sortrows(T, {'Sistem', 'dE_kcal'});
disp(T);
if doPlot
    K = T(T.Sistem == "Lu-ciklam kompleks", :);
    if ~isempty(K)
        figure('Name', 'Lu-ciklam: relativne energije');
        barh(categorical(K.Struktura, flipud(K.Struktura)), flipud(K.dE_kcal), 'FaceColor', [0.6 0.6 0.6]);
        xlabel('\DeltaE PBE0-D3 (kcal/mol)'); title('Relativne energije Lu(III)-ciklam struktura'); grid on;
    end
end
end
