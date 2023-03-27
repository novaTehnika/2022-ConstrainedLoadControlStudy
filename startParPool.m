function [] = startParPool(n_cpu)
    flag = 1;
    while(flag)
        try
            parpool('local',n_cpu);
            flag = 0;
        catch
            flag = 1;
        end
    
    end
end