
files = ls;
nfiles = size(files,1);
for j = 1:nfiles
    display(['file ',num2str(j),' of ',num2str(nfiles)])
    if strfind(files(j,:),"data_study_loadScheduleConstraints")

        load(files(j,:))
%         PP_array(iVar) = model_OLloadcontrol(tMPC(1),y0,Tpto,tMPC(end),par,1);
        PP_array(iVar) = PP;
        dur_array(iVar) = dur;
    end

end


%%
if 0
files = ls;
nfiles = size(files,1);
notDone = 1:100;
for j = 1:nfiles

    if strfind(files(j,:),"data_study_loadScheduleConstraints")
        load(files(j,:))
        [r,c,val] = find(notDone==iVar);
        notDone = [notDone(1:c-1), notDone(c+1:end)];
        
    end

end

try
jobArrayStr = num2str(notDone(1));
for j = 2:length(notDone)
    jobArrayStr = append(jobArrayStr,[',',num2str(notDone(j))]);
        PP_array(notDone(j)) = nan;
        dur_array(notDone(j)) = nan;
end
catch
end

end

%% Transform data to 3D variable mesh
I = length(deltat_Tmax);
J = length(T_max);
K = length(lbFrac);
for i = 1:I
    for j = 1:J
        for k = 1:K
            m = J*K*(i-1) + K*(j-1) + k;
            PP_3D(i,j,k) = PP_array(m);
        end
    end
end


return

%%
figure
scatter3(meshVar.deltat_Tmax(:),meshVar.T_max(:),1e-3*PP_array)
xlabel('max rate of change, from min to max (s)')
ylabel('Torque, max (Nm)')
zlabel('mean power absorbed (kW)')
title([num2str(lbFrac) ,' turndown'])

% figure
% scatter3(meshVar.deltat_Tmax(:),meshVar.T_max(:),1/60*dur_array)
% xlabel('max rate of change, from min to max (s)')
% ylabel('Torque, max (Nm)')
% zlabel('optimization duration (min)')

return

%% Plot power versus torque with each PTO turndown speed having same market
% and each lower bound fraction on its own figure; 
black = [0 0 0];
maroon = [122 0 25]/256;
gold = [255 204 51]/256;
blue = [0 75 135]/256;
orange = [226 100 55]/256;
green = [63 150 87]/256;
color = [maroon; gold; blue; orange; green]

bottomEdge = 1;
leftEdge = 3;
width = 3.5625; % one column: 3+9/16, two column: 7.5
height = 2.75;
fontSize = 9;
lineWidth = 1;

load("data_coulombPTO_dampingStudy_31-Aug-2022_2_slim.mat")
clearvars leg
for k = 1:K
    fig = figure;
    fig.Units = 'inches';
    fig.Position = [leftEdge bottomEdge width height ];

    n_plots = 1;
    ax1 = subplot(n_plots,1,1);
    ax1.FontName = 'times';
    ax1.FontSize = fontSize-1;

    hold on
    
    for i = 1:I
        scatter(1e-6*T_max,1e-3*PP_3D(i,:,k),50, ...
            'filled','x','LineWidth',2,'MarkerEdgeColor',color(i,:))
        
        legLabels(i) = convertCharsToStrings( ...
            ['min. load trans. time = ',num2str(deltat_Tmax(i)),'s']);
    end
    plot(1e-6*T_c_data(iSS,:),1e-3*PP_w_data(iSS,:),'-k','LineWidth',1)
    legLabels(i+1) = "Coulomb damping";
    xlabel('torque, max (MNm)', ...
    'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
    ylabel('power (kW)', ...
    'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
    title(['Mean Power Capture:',newline, ...
           'Minimum Load Fraction of ',num2str(lbFrac(k))],...
    'Interpreter','latex','FontSize',fontSize,'fontname','Times')

    leg = legend(legLabels)
    leg.FontSize = fontSize-1;
    leg.FontName = 'Times';
%     rect = [0.56, 0.2, 0.25, 0.15];
%     set(leg, 'Position', rect)
    set(leg, 'Location', 'best')
    ylim([0 90])
end

%% Determine index from 2D variable mesh
J = length(T_max);
K = length(deltat_Tmax);
j = 11; % T_max index
k = 4; % deltat_Tmax index
i = K*(j-1) + k

T_max(j)
meshVar.T_max(i)
deltat_Tmax(k)
meshVar.deltat_Tmax(i)

%% add Coulomb damping results
SS=7;
plot(1e-6*T_c_data(iSS,:),1e-3*PP_w_data(iSS,:),'k')
legend(['turndown rate = ',num2str(deltat_Tmax(1)),' seconds'],...
    ['turndown rate = ',num2str(deltat_Tmax(K-1),2),' seconds'],...
    ['No turndown'])

%% 
figure
plot(tMPLS,Tpto)
yLim = ylim;
ylim([0 yLim(end)])

%% Package power data (add i-th data set based based on lbFrac)
i = 2;
data(i).PP_2D = PP_3D;
data(i).deltat_Tmax = meshVar.deltat_Tmax;
data(i).T_max = meshVar.T_max;
data(i).lbFrac = lbFrac;
data(i).par = par;

%% Package power data (add data set for Coulomb damping)
i = 4;
data(i).PP_2D = PP_w_data(7,:);
data(i).deltat_Tmax = nan;
data(i).T_max = T_c_data(7,:);
data(i).lbFrac = 1;
data(i).par = nan;

%% Package load Schedule data
i = 3;
data(i).tMPLS = tMPLS;
data(i).Tpto = Tpto;
data(i).lbFrac = lbFrac;
data(i).T_max = T_max(j);
data(i).deltat_Tmax = deltat_Tmax(k);

%% Plot average power results for all cases
k = 1; %index for turndown rate
deltat_Tmax = data(1).deltat_Tmax(k,1);

legLabels = cell(1,length(data));
figure
for i = 1:length(data)-1
    scatter(1e-6*data(i).T_max(k,:),1e-3*data(i).PP_2D(k,:),'filled')
    legLabels{i} = [num2str(data(i).lbFrac),' turndown',];
    hold on
end
i = length(data);
plot(1e-6*data(i).T_max,1e-3*data(i).PP_2D,'k')
legLabels{i} = ['fixed torque'];

xlabel('torque, max (MNm)')
ylabel('power (kW)')
legend(legLabels)

title(['Mean Power Absorbed:',newline,'Comparison for Turndown Rate of ',num2str(deltat_Tmax),' Seconds'])

%% Plot average power results for one turndown fraction but multiple ramp rates
k = 3; %index for turndown ratio
deltat_Tmax = data(1).deltat_Tmax(k,1);

legLabels = cell(1,length(data(k).T_max(:,1)));
figure
for i = 1:length(data(k).T_max(:,1))
    scatter(1e-6*data(k).T_max(i,:),1e-3*data(k).PP_2D(i,:),'filled')
    legLabels{i} = [num2str(data(k).deltat_Tmax(i,1),3),' sec. ramp rate',];
    hold on
end
legLabels{i+1} = ['fixed torque'];
i = length(data);
plot(1e-6*data(i).T_max,1e-3*data(i).PP_2D,'k')


xlabel('torque, max (MNm)')
ylabel('power (kW)')
legend(legLabels)

title(['Mean Power Absorbed:',newline,'Comparison for Turndown Ratio of ',num2str(data(k).lbFrac)])

%% Plot load schedule results
k = 1; %index for turndown rate
deltat_Tmax = data(1).deltat_Tmax(k,1);

legLabels = cell(1,length(data)+1);
figure
for i = 1:length(data)
    plot(data(i).tMPLS,1e-6*data(i).Tpto,'linewidth',1.5)
    legLabels{i} = [num2str(data(i).lbFrac),' turndown',];
    hold on
end
i = length(data)+1;
plot(data(1).tMPLS([1 end]),1e-6*data(1).T_max*ones(2,1),'k','linewidth',1.5)
legLabels{i} = ['max torque'];

xlabel('time (s)')
ylabel('torque (MNm)')
legend(legLabels)

title(['Load Schedule:',newline,'Comparison for Turndown Rate of ',num2str(deltat_Tmax),' Seconds'])

xlim([1020 1180])
yLim = ylim;
ylim([0 yLim(end)])
