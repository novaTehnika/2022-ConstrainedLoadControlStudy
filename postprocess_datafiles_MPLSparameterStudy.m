
files = ls;
nfiles = size(files,1);
for j = 1:nfiles

    if strfind(files(j,:),"data_test_MPLSparameters")
        display(['file ',num2str(j),' of ',num2str(nfiles)])
        load(files(j,:))
%         PP_array(iVar) = model_OLloadcontrol(tMPC(1),y0,Tpto,tMPC(end),par,1);
        PP_array(iVar) = PP;
        dur_array(iVar) = dur;
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


%%
for j = 1:length(dt_ctrl)
    for k = 1:length(tc)
        i = length(tc)*(j-1) + k;
        PP_2D(k,j) = PP_array(i);
        dur_2D(k,j) = dur_array(i);
    end
end

figure
contourf(dt_ctrl,tc,1e-3*PP_2D,'ShowText','on')
xlabel('control update period (s)')
ylabel('control horizon (s)')
title('mean power absorbed (kW)')

figure
contourf(dt_ctrl,tc,dur_2D/60,'ShowText','on')
xlabel('control update period (s)')
ylabel('control horizon (s)')
title('optimization duration (min)')