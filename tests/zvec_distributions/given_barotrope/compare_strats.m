function compare_strats(degree)

close all
if nargin == 0 || isempty(degree), degree = 2; end
if degree == 2, junod = 0.1087e-6; end
if degree == 4, junod = 0.0944e-6; end
if degree == 6, junod = 0.0695e-6; end

files = uigetfile('*.mat','multi','on');
if ~iscell(files), files = {files}; end

figure;
set(gcf, 'defaultLegendInterpreter', 'none')

% first a straight plot
ah1 = axes;
ah1.Box = 'on';
hold(ah1, 'on')
set(ah1,'XTick',2.^(7:16))
set(ah1,'XScale', 'log','XMinorTick','on')
styles = {'--+', '--o', '--*','--x','--s','--d','--p','--h'};
for k=1:length(files)
    data = load(files{k});
    x = data.x;
    if degree == 2, y = data.y2; end
    if degree == 4, y = data.y4; end
    if degree == 6, y = data.y6; end
    plot(x, y, styles{k},...
        'displayname', data.name,...
        'MarkerSize',8,'LineWidth',2)
end
try
    hline(trueval, 'k--');
catch
end
ylabel(['J_',int2str(degree),'(N)']);
xlabel('N levels');
legend show
legend(gca,'location','southeast')

% then a convergence plot
figure
set(gcf, 'defaultLegendInterpreter', 'none')
ah2 = axes;
ah2.Box = 'on';
hold(ah2, 'on')
set(ah2,'XTick',2.^(7:16))
set(ah2,'XScale', 'log','XMinorTick','on')
set(ah2,'YScale', 'log','YMinorTick','on')
styles = {'--+', '--o', '--*','--x','--s','--d','--p','--h'};
for k=1:length(files)
    data = load(files{k});
    x = data.x;
    if degree == 2, y = data.y2; end
    if degree == 4, y = data.y4; end
    if degree == 6, y = data.y6; end
    plot(x, abs(y - y(end)), styles{k},...
        'displayname', data.name,...
        'MarkerSize',8,'LineWidth',2)
end
hline(junod, 'k--');
ylabel(['\Delta J_',int2str(degree),'(N)']);
xlabel('N levels');
legend show
