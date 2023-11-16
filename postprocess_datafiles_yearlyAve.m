
files = ls;
nfiles = size(files,1);
for j = 1:nfiles
display(['file ',num2str(j),' of ',num2str(nfiles)])
    if strfind(files(j,:),"data_study_yearlyAve")
        load(files(j,:))
%         PP_array(iVar) = model_OLloadcontrol(tMPC(1),y0,Tpto,tMPC(end),par,1);
        PP_array(iVar) = PP;
        dur_array(iVar) = dur;
    end

end
load('SSdata_HumboltBay_1D.mat', 'weight')
SS = 1:114;

%%
if 0

files = ls;
nfiles = size(files,1);
notDone = 1:4560;
for j = 1:nfiles
display(['file ',num2str(j),' of ',num2str(nfiles)])
    if strfind(files(j,:),"data_study_yearlyAve")
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

end
%% Transform data to 3D variable mesh
I = length(SS);
J = length(T_max);
K = length(lbFrac);
test = 1;
for i = 1:I
    for j = 1:J
        for k = 1:K
            m = J*I*(k-1) + I*(j-1) + i;
            PP_3D(i,j,k) = PP_array(m);
            SS_3D(i,j,k) = NDmeshVar.SS(m);
            T_max_3D(i,j,k) = NDmeshVar.T_max(m);
            lbFrac_3D(i,j,k) = NDmeshVar.lbFrac(m);
            test = (NDmeshVar.SS(m) == SS(i)) ...
                && (NDmeshVar.T_max(m) == T_max(j)) ...
                && (NDmeshVar.lbFrac(m) == lbFrac(k)) ...
                && test;
        end
    end
end

if ~test; error('indexing incorrect'); end
clearvars test

%% Compare data to prior results for SS7
% black = [0 0 0];
% maroon = [122 0 25]/256;
% gold = [255 204 51]/256;
% blue = [0 75 135]/256;
% orange = [226 100 55]/256;
% green = [63 150 87]/256;
% color = [maroon; gold; blue; orange; green];
% 
% bottomEdge = 1;
% leftEdge = 3;
% width = 3.5625; % one column: 3+9/16, two column: 7.5
% height = 2.75;
% fontSize = 8;
% lineWidth = 1;
% 
% PP_2D_ss7 = squeeze(PP_3D(7,:,:));
% figure;
% clearvars leg
% for k = 1:K
%     scatter(1e-6*T_max,1e-3*PP_2D_ss7(:,k),50, ...
%         'filled','x','LineWidth',2,'MarkerEdgeColor',color(k,:))
%     hold on
%     legLabels(k) = convertCharsToStrings( ...
%         ['turndown ratio = ',num2str(lbFrac(k))]);
% end
% xlabel('torque, max (MNm)', ...
% 'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
% ylabel('power (kW)', ...
% 'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
% title(['Mean Power Capture:',newline, ...
%        'SS 7'],...
% 'Interpreter','latex','FontSize',fontSize,'fontname','Times')
% 
% leg = legend(legLabels)
% 
% 
%%
% figure
% plot(tMPLS,Tpto)

%% Calculate weighted average power production 
for j = 1:J
    for k = 1:K
        PP_yearlyAve(j,k) = sum(weight/100.*PP_3D(:,j,k)');
        
    end
end

%% Calculate optimal weighted average power production for coulomb damping
% case 1: lower frac. load limit of 0
% case 2: lower frac. load limit of 0.25
% case 3: lower frac. load limit of 0.5
% case 4: lower frac. load limit of 0.75
% case 5: fixed load across all sea states

% load("data_coulombPTO_dampingStudy_20220927_slim.mat")
load("data_coulombPTO_dampingStudy_08-May-2023_slim.mat")

nTmax_Coulomb = 1000;
Tmax_Coulomb = linspace(min(T_c_data(:)),max(T_max(:)),nTmax_Coulomb);

% Map Coulomb damping power results to "Tmax_Coulomb" with extrapolation
% beyond the maximum the limit of "T_c_data" for each sea state
PP_w_extrap = zeros(114,nTmax_Coulomb);
for iSS = 1:114
    PP_w_extrap(iSS,:) = interp1(T_c_data(iSS,:),PP_w_data(iSS,:),Tmax_Coulomb,'linear','extrap');
    PP_w_extrap(iSS,:) = (PP_w_extrap(iSS,:)>0).*PP_w_extrap(iSS,:);
end

% initalize
PP_yearlyAveCoulomb = zeros(length(lbFrac) + 1,nTmax_Coulomb);
PPmax = 0;

% specify whether extrapolated data will count toward the weighted average
allowExtrap = 1;

% loop through lower fractional load limit
for ilbFrac = 1:length(lbFrac) + 1
    % loop through max torque
    for iTmax = 1:nTmax_Coulomb
        % loop through sea state
        for iSS = 1:114
            % test for lbFrac or fixed load
            if ilbFrac > length(lbFrac)
                % test for load out of bounds of "T_c_data"
                if T_c_data(iSS,end) < Tmax_Coulomb(iTmax)
                    if allowExtrap
                        PPmax = PP_w_extrap(iSS,iTmax);
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
                    PPfilt = PP_w_extrap(iSS,cond1 & cond2);
                    PPmax = max(PPfilt);
                else
                    PPmax = 0;
                end
            end
            PP_yearlyAveCoulomb(ilbFrac,iTmax) = PP_yearlyAveCoulomb(ilbFrac,iTmax) + weight(iSS)/100*PPmax;
        end
    end
end

%% Plot average power results for one turndown fraction but multiple turndown ratio
black = [0 0 0];
maroon = [122 0 25]/256;
gold = [255 204 51]/256;
blue = [0 75 135]/256;
orange = [226 100 55]/256;
green = [63 150 87]/256;
color = [maroon; gold; blue; orange; green];

bottomEdge = 1;
leftEdge = 3;
width = 3.5625; % one column: 3+9/16, two column: 7.5
height = 3.25;
fontSize = 9;
lineWidth = 1;

clearvars leg

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

n_plots = 1;
ax1 = subplot(n_plots,1,1);
ax1.FontName = 'times';
ax1.FontSize = fontSize-1;

hold on

% dummy plots for legend
iLeg = 0;

for k = 1:K
    scatter(-99*[1, 0.5],-99*[1, 0.5],50, ...
        'filled','s','LineWidth',2,...
        'MarkerEdgeColor',color(k,:),'MarkerFaceColor',color(k,:));
    iLeg = iLeg+1;
    legLabels(iLeg) = convertCharsToStrings( ...
        ['min. load frac. = ',num2str(lbFrac(k))]);
end

scatter(-99*[1 1],-99*[1 1],50, ...
        'filled','x','LineWidth',2,'MarkerEdgeColor','k');
iLeg = iLeg+1;
legLabels(iLeg) = convertCharsToStrings('model predictive control');

plot(-99*[1 1],-99*[1 1],'-k','LineWidth',1)
iLeg = iLeg+1;
legLabels(iLeg) = convertCharsToStrings('optimal Coulomb damping');

plot(-99*[1 1],-99*[1 1],'--k','LineWidth',1)
iLeg = iLeg+1;
legLabels(iLeg) = "fixed Coulomb damping";

% plot real data

for k = 1:K
    s(k) = scatter(1e-6*T_max,1e-3*PP_yearlyAve(:,k),50, ...
        'filled','x','LineWidth',2,...
        'MarkerEdgeColor',color(k,:),'MarkerFaceColor',color(k,:));
    s(k).HandleVisibility='off';
    
    p(k) = plot(1e-6*Tmax_Coulomb,1e-3*PP_yearlyAveCoulomb(k,:),'-','Color',color(k,:),'LineWidth',1);
    p(k).HandleVisibility='off';
end
p(K+1) = plot(1e-6*Tmax_Coulomb,1e-3*PP_yearlyAveCoulomb(end,:),'--k','LineWidth',1);
p(K+1).HandleVisibility='off';

xlabel('torque, max (MNm)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('power (kW)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
title(['Mean Power Capture, Yearly Average'],...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

leg = legend(legLabels)
leg.FontSize = fontSize-1;
leg.FontName = 'Times';
rect = [0.5, -0.2, 0.25, 0.15];
% set(leg, 'Position', rect)
set(leg, 'Location', 'best')
% set(leg, 'Location', 'southoutside')
xlim([0 1e-6*max(T_max)])
