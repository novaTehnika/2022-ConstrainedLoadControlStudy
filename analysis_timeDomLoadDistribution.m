% Start with data loaded from load scheduling run

% perform 1D kernal density estimation on the load schedule
 % Create high resolution load schedule from control updates
t = 0:1e-2*par.dt_ctrl:tMPLS(end);
loadSchedule = size(t);

for it = 1:length(t)
    loadSchedule(it) = Tpto_ramp(t(it),Tpto,par.dt_ctrl);
end

 % check resulting time series
figure
scatter(tMPLS,Tpto,'rx')
hold on
plot(t,loadSchedule,'k-')
    
 % set up inputs to kde()
data = loadSchedule;
n = 1e3;
min = par.T_min;
max = par.T_max;

[bandwidth,density,xmesh,cdf]=kde(data,n,min,max);

% plot density against torque level as continuous curve
% 

figure
yyaxis left
plot(1e-6*xmesh,density)

yyaxis right
plot(1e-6*xmesh,cdf)


[bandwidth,density,xmesh,cdf]=kde(Tpto,n,min,max);