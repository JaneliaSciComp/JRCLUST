function inClust = MikeSplit(mrSpkWav, mrFet, nSplits)

if nargin<1, mrSpkWav = rand(100, 40); mrFet = rand(100,3); nSplits = 2; end
if nargin<2, mrFet = []; end
if nargin<3, nSplits = 2; end

h.mrSpkWav = mrSpkWav';
h.mrFet = mrFet;
h.nSplits = nSplits;

inClust = true(size(h.mrSpkWav, 1), 1);

cs.Nclust = str2double(inputdlg({'Number of clusters:'},'Input', 1, {'2'}));
if isempty(cs.Nclust)
    return;
end

cs.fig = figure(54354); 
set(cs.fig, 'Units', 'Normalized', 'Position', [0.05 0.05 0.8 0.8]);

cs.ax(1) = axes;
cs.ax(2) = axes;
cs.ax(3) = axes;

set(cs.ax(1), 'Units', 'Normalized', 'OuterPosition', [0 0 0.7 1]);
hold(cs.ax(1), 'on');

set(cs.ax(2), 'Units', 'Normalized', 'OuterPosition', [0.7 0.65 0.3 0.35]);
set(cs.ax(3), 'Units', 'Normalized', 'OuterPosition', [0.7 0.4 0.3 0.25]);


cs.ClustList = uicontrol('Style', 'listbox', 'Units', 'normalized', 'Position', ...
    [0.7 0.05 0.1 0.25], 'String', num2str((1:cs.Nclust)'), 'Value', 1, 'BackgroundColor', [1 1 1], ...
    'Max', 3, 'Min', 1, 'Callback',  {@CSUpdatePlots,cs.fig, h});


cs.Done = uicontrol('Style', 'pushbutton',  'Units', 'normalized', 'Position', ...
    [0.85 0.2 0.1 0.1], 'String', 'Done','FontWeight','Normal', 'Callback', ...
    {@CSDone,cs.fig, h});


drawnow;
pause(0.1);

cs.toTestix = true(size(h.mrFet, 1), 1);
cs.toTestNum = find(cs.toTestix);


dat = h.mrFet(cs.toTestix, :);
cs.spk.time = dat(:,1)./25000;
for i = 1:size(dat, 2)
    dat(:,i) = (dat(:,i) - mean(dat(:,i)))./std(dat(:,i));
end

% dat = h.mrFet(cs.toTestix, 2:end);


cs.clustnum = clusterdata(dat, 'linkage', 'ward', 'maxclust',cs.Nclust, 'savememory', 'on');
% cs.clustnum = kmeans(dat, cs.Nclust);

guidata(cs.fig, cs);
CSUpdatePlots([], [], cs.fig, h);

 
uiwait(cs.fig);
 
if ~ishandle(3918)
    return;
end

cs = guidata(3918);
 
if ishandle(cs.fig)
    jrclust.utils.tryClose(cs.fig);
end
if ishandle(3918)
    jrclust.utils.tryClose(3918);
end

inClust = false(size(h.mrSpkWav, 1), 1);
inClust(cs.toTestNum(ismember(cs.clustnum, cs.keep))) = 1;




function CSUpdatePlots(~, ~, fig, h)

cs = guidata(fig);

keep = get(cs.ClustList, 'Value');
clrs = get(gca,'ColorOrder');

inClust = false(size(h.mrSpkWav, 1), 1);
inClust(cs.toTestNum(ismember(cs.clustnum, keep))) = 1;

ISI = diff([cs.spk.time(inClust); inf]);
% ISI(diff(h.spk.trial(inClust))~=0) = inf;

if isempty(ISI)
    ISI = 1;
end

ISIedges = -0.02:0.0005:0.02;
Nisi = histc([ISI; -ISI], ISIedges);


mx = max(Nisi);
mx = 1.05.*mx;

hold(cs.ax(2), 'off');
thr = 0.0025;
axes(cs.ax(2));
f = fill([-thr thr thr -thr], [0 0 mx mx], [1 0.8 0.8]);
set(f, 'Linestyle', 'none');
hold(cs.ax(2), 'on');
b = bar(cs.ax(2), ISIedges+mean(diff(ISIedges))./2, Nisi);
set(b, 'LineStyle', 'none', 'BarWidth', 1, 'FaceColor', 'r');

xlim(cs.ax(2), [-0.02 0.02]);



hold(cs.ax(1), 'off');
dat = h.mrFet(cs.toTestNum, :);
ds = ceil(numel(inClust)/1e4);

ISI = diff([cs.spk.time(inClust); inf]);
viol = false(size(ISI));
viol(1:end-1) = ISI(1:end-1)<thr;
viol(2:end)   = viol(2:end)|ISI(1:end-1)<thr;

violAll = false(size(inClust));
violAll(inClust) = viol;

str = cell(cs.Nclust, 1);
for i = 1:cs.Nclust
    if ismember(i, keep)
        sz = 10;
    else
        sz = 5;
    end
    
    ix = cs.clustnum==i;
    xdat = dat(ix, 1);
    ydat = dat(ix, 2);
    
    plot(cs.ax(1), xdat(1:ds:end), ydat(1:ds:end), '.', 'Color', clrs(mod(i-1, size(clrs,1))+1, :), 'MarkerSize', sz);
%     plot3(cs.ax(1), 1:sum(ix), dat(ix,1), dat(ix, 2), '.', 'Color', clrs(mod(i-1, size(clrs,1))+1, :), 'MarkerSize', sz);
% 
    hold(cs.ax(1), 'on');
    
    str{i} = ['Cluster ' num2str(i) ' (' num2str(sum(ix)) ' spks)'];
end

for i = 1:cs.Nclust

    
    ix = cs.clustnum==i;
    xdat = dat(ix&violAll, 1);
    ydat = dat(ix&violAll, 2);
    
    plot(cs.ax(1), xdat, ydat, '.', 'Color', clrs(mod(i-1, size(clrs,1))+1, :), 'MarkerSize', sz);
%     plot3(cs.ax(1), 1:sum(ix), dat(ix,1), dat(ix, 2), '.', 'Color', clrs(mod(i-1, size(clrs,1))+1, :), 'MarkerSize', sz);

    plot(cs.ax(1), xdat, ydat, 'ko', 'LineWidth', 2);
    hold(cs.ax(1), 'on');
    
    str{i} = ['Cluster ' num2str(i) ' (' num2str(sum(ix)) ' spks)'];
end



legend(cs.ax(1), str, 'Location', 'Best');
set(cs.ax(1), 'YDir', 'Reverse');
view(cs.ax(1), 0, 90);
% xlim(cs.ax(1), prctile(dat(:, 1), [0.001 99.999]));
% ylim(cs.ax(1), prctile(dat(:, 2), [0.001 99.999]));

xlim(cs.ax(1), [min(dat(:,1)) max(dat(:,1))]);
ylim(cs.ax(1), [min(dat(:,2)) max(dat(:,2))]);

hold(cs.ax(3), 'off');
Nplot = 50;
for i = 1:cs.Nclust
    
    if ismember(i, keep)
        c = clrs(mod(i-1, size(clrs,1))+1, :)+(1-clrs(mod(i-1, size(clrs,1))+1, :))./2;
    else
        c = [0.75 0.75 0.75];
    end
    
    ix = cs.toTestNum(cs.clustnum==i);
    ix = randsample(ix, min(Nplot, numel(ix)));
    plot(cs.ax(3), h.mrSpkWav(ix, :)', 'Color', c', 'LineWidth', 0.5);

    hold(cs.ax(3), 'on');
end

for i = 1:cs.Nclust
    
    if ismember(i, keep)
        c = clrs(mod(i-1, size(clrs,1))+1, :);
    else
        c = [0.5 0.5 0.5];
    end
    ix = cs.toTestNum(cs.clustnum==i);
    spks = h.mrSpkWav(ix, :);
    
    plot(cs.ax(3), mean(spks, 1), 'Color', c, 'LineWidth', 2);
    
end




function CSDone(~, ~, fig, ~)

cs = guidata(fig);
cs.keep = get(cs.ClustList, 'Value');

if ishandle(fig)
    jrclust.utils.tryClose(fig);
end
fig = figure(3918);

guidata(fig, cs);

