function compare_strats()

files = uigetfile('*.mat','multi','on');
if ~iscell(files), files = {files}; end
close all

% NMoI plot
figure
set(gcf, 'defaultLegendInterpreter', 'none')
ah1 = axes;
ah1.Box = 'on';
hold(ah1, 'on')
set(ah1,'XTick',[128 256 512 1024 2048 4096 8192 16384])
set(ah1,'XScale', 'log','XMinorTick','on')
styles = {'--+', '--o', '--*','--x','--s','--d','--p','--h'};
for k=1:length(files)
    data = load(files{k});
    x = data.x;
    y = data.L;
    plot(x, y, styles{k},...
        'displayname', data.name,...
        'MarkerSize',8,'LineWidth',2)
end
ylabel('NMoI');
xlabel('N levels');
legend show

% betanorm plot
figure
set(gcf, 'defaultLegendInterpreter', 'none')
ah1 = axes;
ah1.Box = 'on';
hold(ah1, 'on')
set(ah1,'XTick',[128 256 512 1024 2048 4096 8192 16384])
set(ah1,'XScale', 'log','XMinorTick','on')
styles = {'--+', '--o', '--*','--x','--s','--d','--p','--h'};
for k=1:length(files)
    data = load(files{k});
    x = data.x;
    y = data.beta;
    plot(x, y, styles{k},...
        'displayname', data.name,...
        'MarkerSize',8,'LineWidth',2)
end
ylabel('betanorm');
xlabel('N levels');
legend show
