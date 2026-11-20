function fig = plotProgression(out, varargin)
%FLPM.PLOTPROGRESSION  Show how the cracks grow, interact and link up.
%
%   r   = flpm.fatigueLife(dSigma, defects, mat, flpm.options('RecordHistory', true));
%   fig = flpm.plotProgression(r)
%   fig = flpm.plotProgression(r, 'Frames', 5, 'Window', [s0 s1 z0 z1])
%   flpm.plotProgression(r, 'Animate', 'progression.gif')
%
%   Left: the fracture plane (normal to the load), position along the
%   surface against depth below it, with the free surface at depth 0. Each
%   crack is the disc of its sqrt(area), radius sqrt(area/pi), seen along
%   the load axis (axial offsets are not shown). Outlines show the cracks at
%   the selected moments, light to dark in time; the last of them is
%   filled, with the crack of largest dK in red. Discs cut by the surface
%   are clipped at depth 0. By default the moments cover the interaction
%   stage: until the largest crack is 4 times its initial size, or up to
%   the first link-up if that is later; the long-crack stage that follows
%   is in the right panel. Right: the largest crack against cycles (log
%   scale), with the shown moments marked and link-ups as dotted lines.
%
%   Name-value pairs
%     'Frames'   number of moments to show, or indices into r.frames.
%                Default: start, every link-up and the end of the
%                interaction stage, plus evenly spaced moments, 6 in all
%     'Window'   [s0 s1 z0 z1] (mm) of the left panel. Default: around the
%                governing crack at the end of the interaction stage and
%                the defects it absorbed (other growing cracks may lie
%                outside; pass a window to see them)
%     'Animate'  GIF file: the moments in turn (at most 'MaxImages', most
%                of them in the interaction stage), the left panel at that
%                moment and the right panel up to its cycles
%     'MaxImages' images in the GIF                                    [60]
%     'Delay'    seconds per GIF image                                [0.15]
%     'Visible'  figure visibility, 'on' or 'off'                       ['on']
%
%   See also flpm.fatigueLife, flpm.options.

ip = inputParser;
ip.FunctionName = 'flpm.plotProgression';
ip.addParameter('Frames', []);
ip.addParameter('Window', []);
ip.addParameter('Animate', '');
ip.addParameter('Delay', 0.15);
ip.addParameter('MaxImages', 60);
ip.addParameter('Visible', 'on');
ip.parse(varargin{:});
o = ip.Results;
if ~isfield(out, 'frames') || isempty(out.frames)
    error('flpm:plotProgression:noFrames', ...
        'Run flpm.fatigueLife with flpm.options(''RecordHistory'', true).');
end
F = out.frames;
n = numel(F);
link = find(arrayfun(@(k) numel(F(k).id) < numel(F(k-1).id), 2:n)) + 1;
am = arrayfun(@(f) max(f.size), F);
stage = find(am <= 4 * am(1), 1, 'last');
if ~isempty(link), stage = max(stage, link(1)); end

% moments to show
if isempty(o.Frames) || isscalar(o.Frames)
    m = 6;
    if isscalar(o.Frames), m = o.Frames; end
    key = unique([1, link(link <= stage), stage]);
    if numel(key) > m
        key = key(unique(round(linspace(1, numel(key), m))));
    end
    extra = round(linspace(1, stage, m));
    k = unique([key, extra(~ismember(extra, key))]);
    while numel(k) > m
        free = find(~ismember(k, key));
        [~, j] = min(diff(k(free)));
        if isempty(j), break, end
        k(free(j + 1)) = [];
    end
    sel = k;
else
    sel = unique(o.Frames(:).');
end
sel = sel(sel >= 1 & sel <= n);

win = o.Window;
if isempty(win), win = defaultWindow(F(1), F(stage)); end

fig = figure('Color', 'w', 'Visible', o.Visible, 'Position', [100 100 1000 420]);
if exist('theme', 'file') == 2 || exist('theme', 'builtin') == 5
    try, theme(fig, 'light'); catch, end
end
t = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax1 = nexttile(t);
ax2 = nexttile(t);
cmap = shades(numel(sel));

drawPlane(ax1, F, sel, cmap, win);
drawCurve(ax2, out, F, sel, link, cmap, numel(F));
title(t, sprintf('\\Delta\\sigma = %.0f MPa: %s', out.stress, outcome(out)));

if ~isempty(o.Animate)
    late = min(10, floor(o.MaxImages / 4));
    early = round(linspace(1, stage, min(stage, o.MaxImages - late)));
    for j = unique([early, round(linspace(stage, n, late + 1))])
        drawPlane(ax1, F, j, [0.35 0.35 0.35], win);
        drawCurve(ax2, out, F, j, link(link <= j), [0.35 0.35 0.35], j);
        drawnow;
        im = frame2im(getframe(fig));
        [A, map] = rgb2ind(im, 128);
        if j == 1
            imwrite(A, map, o.Animate, 'gif', 'LoopCount', Inf, 'DelayTime', o.Delay);
        else
            imwrite(A, map, o.Animate, 'gif', 'WriteMode', 'append', ...
                'DelayTime', o.Delay * (1 + 14 * (j == n)));
        end
    end
    % leave the static figure in place
    drawPlane(ax1, F, sel, cmap, win);
    drawCurve(ax2, out, F, sel, link, cmap, numel(F));
end
end

% -------------------------------------------------------------------------
function drawPlane(ax, F, sel, cmap, win)
cla(ax); hold(ax, 'on');
th = linspace(0, 2*pi, 90);
last = F(sel(end));
% filled last state
[~, kk] = max(last.dK);
for i = 1:numel(last.size)
    [x, y] = disc(last, i, th);
    col = [0.55 0.65 0.8];
    if i == kk, col = [0.85 0.2 0.15]; end
    patch(ax, x, y, col, 'EdgeColor', 'none', 'FaceAlpha', 0.45);
end
% outlines of the selected moments
for j = 1:numel(sel)
    f = F(sel(j));
    for i = 1:numel(f.size)
        [x, y] = disc(f, i, th);
        plot(ax, x, y, '-', 'Color', cmap(j, :), 'LineWidth', 1);
    end
end
yline(ax, 0, 'k-', 'LineWidth', 1.5);
text(ax, win(1), 0, ' free surface', 'VerticalAlignment', 'bottom', 'FontSize', 8);
set(ax, 'YDir', 'reverse');
axis(ax, 'equal');
xlim(ax, win(1:2)); ylim(ax, [min(win(3), -0.04 * diff(win(3:4))) win(4)]);
xlabel(ax, 'position along the surface (mm)');
ylabel(ax, 'depth below the surface (mm)');
box(ax, 'on'); hold(ax, 'off');
end

function [x, y] = disc(f, i, th)
r = f.size(i) / sqrt(pi);
x = f.position(i) + r * cos(th);
y = max(f.depth(i) + r * sin(th), 0);
end

function drawCurve(ax, out, F, sel, link, cmap, upto)
cla(ax); hold(ax, 'on');
Nf = [F.N];
am = arrayfun(@(f) max(f.size), F);
set(ax, 'YScale', 'log');
if isfield(out, 'history') && ~isempty(out.history.N)
    Nh = [0; out.history.N];
    ah = [am(1); out.history.aMax];
    keep = Nh <= Nf(upto);
    plot(ax, Nh(keep), ah(keep) * 1e3, 'k-', 'LineWidth', 1.2);
else
    plot(ax, Nf(1:upto), am(1:upto) * 1e3, 'k-', 'LineWidth', 1.2);
end
for L = link(:).'
    xline(ax, Nf(L), ':', 'Color', [0.4 0.4 0.4]);
end
for j = 1:numel(sel)
    plot(ax, Nf(sel(j)), am(sel(j)) * 1e3, 'o', 'MarkerFaceColor', cmap(j, :), ...
        'MarkerEdgeColor', 'k', 'MarkerSize', 6);
end
xlabel(ax, 'cycles N');
ylabel(ax, 'largest \surdarea (\mum)');
if Nf(end) > 0, xlim(ax, [0 1.02 * Nf(end)]); end
grid(ax, 'on'); box(ax, 'on'); hold(ax, 'off');
end

function win = defaultWindow(a, b)
% around the governing crack (largest dK) of snapshot b and the original
% defects it has absorbed, with a margin of its radius
[~, k] = max(b.dK);
R = b.size(k) / sqrt(pi);
in = hypot(a.position - b.position(k), a.depth - b.depth(k)) <= R;
r = [R; a.size(in) / sqrt(pi)];
s = [b.position(k); a.position(in)];
z = [b.depth(k); a.depth(in)];
lo = [min(s - r), max(min(z - r), 0)];
hi = [max(s + r), max(z + r)];
pad = max(R, 0.25 * max(hi - lo));
win = [lo(1) - pad, hi(1) + pad, max(lo(2) - pad, 0), hi(2) + pad];
end

function c = shades(m)
% light to dark blue
c = [linspace(0.75, 0.05, m)', linspace(0.82, 0.2, m)', linspace(0.95, 0.5, m)'];
end

function s = outcome(out)
switch out.status
    case 'Kc',       s = sprintf('failure (K_{max} = K_{IC}) at N = %.3g', out.N);
    case 'arrest',   s = 'all cracks arrested (run-out)';
    case 'noGrowth', s = 'no crack grows';
    otherwise,       s = 'run-out';
end
end
