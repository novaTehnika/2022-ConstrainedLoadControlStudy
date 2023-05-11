%%
% Instructions: Add project directory to path, navigate to the to folder
% with data to be processed, and run this script.

%%
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
if 1
files = ls;
nfiles = size(files,1);
notDone = 1:length(meshVar.T_max(:));
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
test = 1;
for i = 1:I
    for j = 1:J
        for k = 1:K
            m = J*I*(k-1) + I*(j-1) + i;
            PP_3D(i,j,k) = PP_array(m);
            test = (meshVar.deltat_Tmax(m) == deltat_Tmax(i)) ...
                && (meshVar.T_max(m) == T_max(j)) ...
                && (meshVar.lbFrac(m) == lbFrac(k)) ...
                && test;
        end
    end
end

if ~test; error('indexing incorrect'); end
clearvars test

return


%% Calculate optimal power production for coulomb damping
% case 1: lower frac. load limit of 0
% case 2: lower frac. load limit of 0.25
% case 3: lower frac. load limit of 0.5
% case 4: lower frac. load limit of 0.75
% case 5: fixed load across all sea states

load("data_coulombPTO_dampingStudy_08-May-2023_slim.mat")

nTmax_Coulomb = 1000;
Tmax_Coulomb = linspace(min(T_c_data(iSS)),max(T_max(:)),nTmax_Coulomb);

% Map Coulomb damping power results to "Tmax_Coulomb" with extrapolation
% beyond the maximum the limit of "T_c_data" for each sea state
PP_w_extrap = interp1(T_c_data(iSS,:),PP_w_data(iSS,:),Tmax_Coulomb,'linear','extrap');
PP_w_extrap = (PP_w_extrap>0).*PP_w_extrap;

% initalize
PP_CDoptimal = zeros(length(lbFrac) + 1,nTmax_Coulomb);

% specify whether extrapolated data will count toward the weighted average
allowExtrap = 1;

% loop through lower fractional load limit
for ilbFrac = 1:length(lbFrac) + 1
    % loop through max torque
    for iTmax = 1:nTmax_Coulomb
        % test for lbFrac or fixed load
        if ilbFrac > length(lbFrac)
            % test for load out of bounds of "T_c_data"
            if T_c_data(iSS,end) < Tmax_Coulomb(iTmax)
                if allowExtrap
                    PPmax = PP_w_extrap(iTmax);
                else
                    PPmax = 0;
                end
            else
                PPmax = interp1(T_c_data(iSS,:),PP_w_data(iSS,:),Tmax_Coulomb(iTmax),'linear');
            end
        else
            % test for the T_c_data spanning the legal load range
            % cond1 - the upper bound for T_c_data exceeds the maximum load
            % cond2 - extrapolation is allowed
            if (T_c_data(iSS,end) >= Tmax_Coulomb(iTmax)) || allowExtrap
                % filter array to include 
                % loads less than or equal to max the load (cond1) and 
                % loads greater than or equal to the min load (cond2)
                cond1 = Tmax_Coulomb <= Tmax_Coulomb(iTmax);
                cond2 = Tmax_Coulomb >= lbFrac(ilbFrac)*Tmax_Coulomb(iTmax);
                PPfilt = PP_w_extrap(cond1 & cond2);
                PPmax = max(PPfilt);
            else
                PPmax = 0;
            end
        end
       PP_CDoptimal(ilbFrac,iTmax) = PPmax; 
    end
end

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

load("data_coulombPTO_dampingStudy_08-May-2023_slim.mat")
% load("data_coulombPTO_dampingStudy_31-Aug-2022_2_slim.mat")
clearvars leg

% Loop through lbfrac
for k = 1:K
    fig = figure;
    fig.Units = 'inches';
    fig.Position = [leftEdge bottomEdge width height ];

    n_plots = 1;
    ax1 = subplot(n_plots,1,1);
    ax1.FontName = 'times';
    ax1.FontSize = fontSize-1;

    hold on
    
    % loop through deltat_Tmax
    for i = 1:3
        scatter(1e-6*T_max,1e-3*PP_3D(i,:,k),50, ...
            'filled','x','LineWidth',2,'MarkerEdgeColor',color(i,:))
        
        legLabels(i) = convertCharsToStrings( ...
            ['min. load trans. time = ',num2str(deltat_Tmax(i)),'s']);
    end

    plot(1e-6*Tmax_Coulomb,1e-3*PP_CDoptimal(k,:),'-k','LineWidth',1)
    legLabels(i+1) = "optimal Coulomb damping";

    plot(1e-6*T_c_data(iSS,:),1e-3*PP_w_data(iSS,:),'--k','LineWidth',1)
    legLabels(i+2) = "fixed Coulomb damping";

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
