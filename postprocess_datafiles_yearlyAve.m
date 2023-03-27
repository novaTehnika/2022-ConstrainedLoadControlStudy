
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


%% Transform data to 3D variable mesh
I = length(SS);
J = length(T_max);
K = length(lbFrac);
for i = 1:I
    for j = 1:J
        for k = 1:K
            m = J*I*(k-1) + I*(j-1) + i;
            PP_3D(i,j,k) = PP_array(m);
            SS_3D(i,j,k) = NDmeshVar.SS(m);
            T_max_3D(i,j,k) = NDmeshVar.T_max(m);
            lbFrac_3D(i,j,k) = NDmeshVar.lbFrac(m);
            III(i,j,k) = m;
        end
    end
end

%% Compare data to prior results for SS7
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
height = 2.75;
fontSize = 8;
lineWidth = 1;

PP_2D_ss7 = squeeze(PP_3D(7,:,:));
figure;
clearvars leg
for k = 1:K
    scatter(1e-6*T_max,1e-3*PP_2D_ss7(:,k),50, ...
        'filled','x','LineWidth',2,'MarkerEdgeColor',color(k,:))
    hold on
    legLabels(k) = convertCharsToStrings( ...
        ['turndown ratio = ',num2str(lbFrac(k))]);
end
xlabel('torque, max (MNm)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('power (kW)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
title(['Mean Power Capture:',newline, ...
       'SS 7'],...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

leg = legend(legLabels)


%%
figure
plot(tMPLS,Tpto)

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

load("data_coulombPTO_dampingStudy_20220927_slim.mat")

nTmax_Coulomb = 100;
Tmax_Coulomb = linspace(min(T_c_data(:)),max(T_max(:)),nTmax_Coulomb);

PP_yearlyAveCoulomb = zeros(length(lbFrac) + 1,nTmax_Coulomb);
PPmax = 0;

% loop through lower fractional load limit
for ilbFrac = 1:length(lbFrac) + 1
    % loop through max torque
    for iTmax = 1:nTmax_Coulomb
        % loop through sea state
        for iSS = 1:114
            % test for lbFrac or fixed load and out of bounds load
            if ilbFrac > length(lbFrac)
                if T_c_data(iSS,end) < Tmax_Coulomb(iTmax)
                    PPmax = 0;
                else
                    PPmax = interp1(T_c_data(iSS,:),PP_w_data(iSS,:),Tmax_Coulomb(iTmax),'linear');
                end
            else
                if T_c_data(iSS,end) < lbFrac(ilbFrac)*Tmax_Coulomb(iTmax)
                    PPmax = 0;
                else
                    PPmax = max(PP_w_data(iSS,T_c_data(iSS,:)<=Tmax_Coulomb(iTmax) & T_c_data(iSS,:)>=lbFrac(ilbFrac)*Tmax_Coulomb(iTmax)));
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
height = 2.75;
fontSize = 8;
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

for k = 1:K
    scatter(1e-6*T_max,1e-3*PP_yearlyAve(:,k),50, ...
        'filled','x','LineWidth',2,'MarkerEdgeColor',color(k,:))
    
    legLabels(k) = convertCharsToStrings( ...
        ['turndown ratio = ',num2str(lbFrac(k))]);
end
plot(1e-6*Tmax_Coulomb,1e-3*PP_yearlyAveCoulomb(:),'-k','LineWidth',1)
legLabels(k+1) = "Coulomb damping";
xlabel('torque, max (MNm)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('power (kW)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
title(['Mean Power Capture, Yearly'],...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

leg = legend(legLabels)
leg.FontSize = fontSize-1;
leg.FontName = 'Times';
rect = [0.56, 0.2, 0.25, 0.15];
set(leg, 'Position', rect)
% ylim([0 90])
