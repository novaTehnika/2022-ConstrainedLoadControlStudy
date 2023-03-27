
files = ls;
nfiles = size(files,1);
for j = 1:nfiles
display(['file ',num2str(j),' of ',num2str(nfiles)])
    if strfind(files(j,:),"data_test_MPLSparameters")
        load(files(j,:))
%         PP_array(iVar) = model_OLloadcontrol(tMPC(1),y0,Tpto,tMPC(end),par,1);
        PP_array(iVar) = PP;
        dur_array(iVar) = dur;
        par.dTdt_max
    end

end



%%

files = ls;
nfiles = size(files,1);
notDone = 1:100;
for j = 1:nfiles

    if strfind(files(j,:),"data_test_MPLSparameters")
        load(files(j,:))
        [r,c,val] = find(notDone==iVar);
        notDone = [notDone(1:c-1), notDone(c+1:end)];
        
    end

end
jobArrayStr = num2str(notDone(1));
for j = 2:length(notDone)
    jobArrayStr = append(jobArrayStr,[',',num2str(notDone(j))]);
        PP_array(notDone(j)) = nan;
        dur_array(notDone(j)) = nan;
end


%%
figure
scatter3(meshVar.dt_ctrl(:),meshVar.tc(:),1e-3*PP_array)
xlabel('control update period (s)')
ylabel('control horizon (s)')
zlabel('mean power absorbed (kW)')

figure
scatter3(meshVar.dt_ctrl(:),meshVar.tc(:),1/60*dur_array)
xlabel('control update period (s)')
ylabel('control horizon (s)')
zlabel('optimization duration (min)')

return

%% remove data point
% tc_rm = [30 30];
% dt_rm = [0.61 0.5];

% tc_rm = [23 ];
% dt_rm = [0.61 ];

for i = 1:length(tc_rm)
    [~, Itc] = min(abs(tc-tc_rm(i)));
    [~, Idt] = min(abs(dt_ctrl-dt_rm(i)));
    iVar = find(meshVar.tc==tc(Itc) & meshVar.dt_ctrl==dt_ctrl(Idt));
    PP_array(iVar) = NaN;
    dur_array(iVar) = NaN;
end

%%
for j = 1:length(dt_ctrl)
    for k = 1:length(tc)
        i = length(tc)*(j-1) + k;
        PP_2D(k,j) = PP_array(i);
        dur_2D(k,j) = dur_array(i);
    end
end
bottomEdge = 1;
leftEdge = 3;
width = 3.5625; % one column: 3+9/16, two column: 7.5
height = 2.5;
fontSize = 8;
lineWidth = 1;

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

n_plots = 1;
ax1 = subplot(n_plots,1,1);
ax1.FontName = 'times';
ax1.FontSize = fontSize-1;

hold on

[C,c1] = contour(ax1,dt_ctrl,tc,1e-3*PP_2D,'-','ShowText','on')
c1.LineWidth = lineWidth;
c1.LineColor = 'k';
clabel(C,c1,'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
xlabel('control update period (s)',...
    'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('control horizon (s)',...
    'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
title('Mean Power Capture (kW)',...
    'Interpreter','latex','FontSize',fontSize,'fontname','Times')


fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

n_plots = 1;
ax1 = subplot(n_plots,1,1);
ax1.FontName = 'times';
ax1.FontSize = fontSize-1;

hold on

[C,c1] = contour(ax1,dt_ctrl,tc,dur_2D/60/60,'-','ShowText','on')
c1.LineWidth = lineWidth;
c1.LineColor = 'k';
clabel(C,c1,'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
xlabel('control update period (s)',...
    'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('control horizon (s)',...
    'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
title('Duration (hr)',...
    'Interpreter','latex','FontSize',fontSize,'fontname','Times')
