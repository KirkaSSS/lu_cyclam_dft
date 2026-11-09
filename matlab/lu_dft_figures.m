function lu_dft_figures(outDir)
%LU_DFT_FIGURES  Pravi Fig. S1 (strukture CN 8 i CN 9) i Fig. S2 (relativne energije) za SI.
%   lu_dft_figures            % PNG (300 dpi) i PDF u <paket>/figures
%   lu_dft_figures('C:\rad\slike')
%   Podaci: results/xtb_structures/*.xyz i results/dft/dft.jsonl iz paketa.
%   Radi u MATLAB R2020a ili novijem (exportgraphics); bez toolbox-ova.
root = lu_dft_root();
if nargin < 1 || isempty(outDir), outDir = fullfile(root, 'figures'); end
if ~exist(outDir, 'dir'), mkdir(outDir); end
xyzDir = fullfile(root, 'results', 'xtb_structures');

% ---- energije iz dft.jsonl (iste kao u Tabeli 1) ----
T = lu_dft_results(fullfile(root, 'results', 'dft', 'dft.jsonl'), false);
K = T(T.Sistem == "Lu-ciklam kompleks", :);
dE = @(name) K.dE_kcal(K.Struktura == name);

% ---- Slika 1: strukture ----
f1 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 25 13]);
ax1 = axes(f1, 'Position', [0.02 0.12 0.46 0.78]);
draw_structure(ax1, fullfile(xyzDir, 's3_2_10.xyz'));
ax2 = axes(f1, 'Position', [0.52 0.12 0.46 0.78]);
draw_structure(ax2, fullfile(xyzDir, 'cn9_s4_1_6.xyz'));
annotation(f1, 'textbox', [0.03 0.90 0.45 0.07], 'String', sprintf('(a) CN 8, %.1f kcal/mol', dE("s3_2_10")), ...
    'EdgeColor', 'none', 'FontSize', 13, 'FontName', 'Helvetica');
annotation(f1, 'textbox', [0.53 0.90 0.45 0.07], 'String', sprintf('(b) CN 9, +%.1f kcal/mol', dE("cn9_s4_1_6")), ...
    'EdgeColor', 'none', 'FontSize', 13, 'FontName', 'Helvetica');
% legenda (nevidljivi markeri u posebnim osama)
axL = axes(f1, 'Position', [0.05 0.01 0.9 0.08], 'Visible', 'off'); hold(axL, 'on');
c = atom_colors();
h = gobjects(7, 1);
h(1) = plot(axL, NaN, NaN, 'o', 'MarkerFaceColor', c.Lu, 'MarkerEdgeColor', 'k', 'MarkerSize', 10, 'LineStyle', 'none');
h(2) = plot(axL, NaN, NaN, 'o', 'MarkerFaceColor', c.N, 'MarkerEdgeColor', 'k', 'MarkerSize', 10, 'LineStyle', 'none');
h(3) = plot(axL, NaN, NaN, 'o', 'MarkerFaceColor', c.O, 'MarkerEdgeColor', 'k', 'MarkerSize', 10, 'LineStyle', 'none');
h(4) = plot(axL, NaN, NaN, 'o', 'MarkerFaceColor', c.C, 'MarkerEdgeColor', 'k', 'MarkerSize', 10, 'LineStyle', 'none');
h(5) = plot(axL, NaN, NaN, 'o', 'MarkerFaceColor', c.H, 'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'LineStyle', 'none');
h(6) = plot(axL, NaN, NaN, '-', 'Color', c.Lu, 'LineWidth', 2.2);
h(7) = plot(axL, NaN, NaN, ':', 'Color', c.bond, 'LineWidth', 1.6);
legend(axL, h, {'Lu', 'N', 'O', 'C', 'H (C-H omitted)', 'Lu-N/O bond', 'H-bond'}, ...
    'Orientation', 'horizontal', 'Location', 'north', 'Box', 'off', 'FontSize', 10);
exportgraphics(f1, fullfile(outDir, 'fig1_structures.png'), 'Resolution', 300);
exportgraphics(f1, fullfile(outDir, 'fig1_structures.pdf'), 'ContentType', 'vector');

% ---- Slika 2: relativne energije ----
names  = ["s3_2_10", "s4_1_6", "cn9_s4_1_6", "cn9_s3_2_10", "s5_0_6"];
labels = {'CN 8 (4N + 4O + H_2O), A', 'CN 8 (4N + 4O + H_2O), B', 'CN 9 (4N + 5O), B', ...
          'CN 9 (4N + 5O), A', 'CN 8 (3N + 5O)'};
v = arrayfun(@(n) dE(n), names);
f2 = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 18 8]);
ax = axes(f2); hold(ax, 'on');
y = numel(v):-1:1;
grey = [0.604 0.631 0.663]; accent = [0.165 0.435 0.859];
for i = 1:numel(v)
    col = grey; if v(i) == min(v), col = accent; end
    barh(ax, y(i), max(v(i), 0.15), 0.58, 'FaceColor', col, 'EdgeColor', 'none');
    if v(i) == min(v), s = sprintf('%.1f (reference)', v(i)); else, s = sprintf('%.1f', v(i)); end
    text(ax, v(i) + 0.35, y(i), s, 'VerticalAlignment', 'middle', 'FontSize', 10);
end
set(ax, 'YTick', fliplr(y), 'YTickLabel', fliplr(labels), 'TickLabelInterpreter', 'tex', ...
    'FontSize', 10, 'TickLength', [0 0], 'Box', 'off', 'XColor', [0.3 0.3 0.3], 'YColor', [0.3 0.3 0.3]);
xlim(ax, [0 22]); ylim(ax, [0.4 numel(v) + 0.6]);
xlabel(ax, '\DeltaE, PBE0-D3(BJ)/def2-SVP, ddCOSMO water (kcal/mol)');
exportgraphics(f2, fullfile(outDir, 'fig2_energies.png'), 'Resolution', 300);
exportgraphics(f2, fullfile(outDir, 'fig2_energies.pdf'), 'ContentType', 'vector');
fprintf('Slike su sacuvane u: %s\n', outDir);
end

% =====================================================================
function draw_structure(ax, file)
[Z, P] = read_xyz(file);
c = atom_colors();
% orijentacija: Lu u koordinatnom pocetku, normala ravni N4 nagore
lu = find(Z == 71, 1); P = P - P(lu, :);
N = P(Z == 7, :); cN = mean(N, 1);
[~, ~, V] = svd(N - cN, 0);
n = V(:, 3)'; if dot(-cN, n) < 0, n = -n; end
e1 = V(:, 1)'; e2 = cross(n, e1);
Q = P * [e1; n; e2]';
t = deg2rad(8); Rx = [1 0 0; 0 cos(t) -sin(t); 0 sin(t) cos(t)];
Q = Q * Rx';
% izbaci H vezane za C
D0 = squareform_dist(Q);
keep = true(numel(Z), 1);
for i = 1:numel(Z)
    if Z(i) == 1 && any(D0(i, Z == 6) < 1.2), keep(i) = false; end
end
Z = Z(keep); Q = Q(keep, :); lu = find(Z == 71, 1);
D = squareform_dist(Q); n = numel(Z);
% veze: [i j tip], tip 1 = kovalentna, 2 = Lu-N/O, 3 = H-veza
B = zeros(0, 3);
for i = 1:n
    for j = i+1:n
        zi = Z(i); zj = Z(j);
        if zi == 71 || zj == 71
            o = j; if zj == 71, o = i; end
            if any(Z(o) == [7 8]) && D(i, j) < 3.0, B(end+1, :) = [i j 2]; end %#ok<AGROW>
        elseif zi == 1 || zj == 1
            if D(i, j) < 1.2, B(end+1, :) = [i j 1]; end %#ok<AGROW>
        elseif D(i, j) < 1.75
            B(end+1, :) = [i j 1]; %#ok<AGROW>
        end
    end
end
dLu = vecnorm(Q - Q(lu, :), 2, 2);
outer = find(Z == 8 & dLu > 3.0)';
for o = outer
    for hh = find(Z == 1)'
        if D(o, hh) > 1.5 && D(o, hh) < 2.3, B(end+1, :) = [o hh 3]; end %#ok<AGROW>
    end
end
% crtanje po dubini (z): dalje prvo
zB = (Q(B(:, 1), 3) + Q(B(:, 2), 3)) / 2;
items = [zB, ones(size(zB)), (1:size(B, 1))'; Q(:, 3) + 0.005, 2 * ones(n, 1), (1:n)'];
items = sortrows(items, 1);
hold(ax, 'on');
for k = 1:size(items, 1)
    if items(k, 2) == 1
        b = B(items(k, 3), :); xy = Q(b(1:2), 1:2);
        switch b(3)
            case 1, plot(ax, xy(:, 1), xy(:, 2), '-', 'Color', c.bond, 'LineWidth', 2.2);
            case 2, plot(ax, xy(:, 1), xy(:, 2), '-', 'Color', c.Lu, 'LineWidth', 2.2);
            case 3, plot(ax, xy(:, 1), xy(:, 2), ':', 'Color', c.bond, 'LineWidth', 1.6);
        end
    else
        i = items(k, 3); [col, r] = atom_style(Z(i), c);
        rectangle(ax, 'Position', [Q(i, 1) - r, Q(i, 2) - r, 2 * r, 2 * r], 'Curvature', [1 1], ...
            'FaceColor', col, 'EdgeColor', [0.13 0.15 0.16], 'LineWidth', 0.7);
    end
end
for o = outer
    text(ax, Q(o, 1) + 0.75, Q(o, 2) + 0.55, 'H_2O (2nd sphere)', 'FontSize', 9, 'Color', [0.29 0.31 0.34]);
end
axis(ax, 'equal'); axis(ax, 'off');
xlim(ax, [min(Q(:, 1)) - 0.8, max(Q(:, 1)) + 2.2]); ylim(ax, [-4.6 4.6]);
end

function [Z, P] = read_xyz(file)
L = splitlines(strtrim(fileread(file)));
nat = str2double(L{1});
sym = {'H', 'C', 'N', 'O', 'Lu'}; num = [1 6 7 8 71];
Z = zeros(nat, 1); P = zeros(nat, 3);
for i = 1:nat
    t = strsplit(strtrim(L{i + 2}));
    Z(i) = num(strcmp(sym, t{1}));
    P(i, :) = str2double(t(2:4));
end
end

function D = squareform_dist(Q)
D = sqrt(max(sum(Q.^2, 2) + sum(Q.^2, 2)' - 2 * (Q * Q'), 0));
end

function c = atom_colors()
c.Lu = [0.165 0.616 0.561]; c.N = [0.231 0.357 0.859]; c.O = [0.878 0.192 0.192];
c.C = [0.525 0.557 0.588]; c.H = [0.945 0.953 0.961]; c.bond = [0.286 0.314 0.341];
end

function [col, r] = atom_style(z, c)
switch z
    case 71, col = c.Lu; r = 0.62;
    case 7,  col = c.N;  r = 0.30;
    case 8,  col = c.O;  r = 0.30;
    case 6,  col = c.C;  r = 0.27;
    otherwise, col = c.H; r = 0.15;
end
end
