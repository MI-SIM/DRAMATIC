%Calculation of TN and removal along branch point curve

%O2-RAMX
tnremove(xbp1,xbp2,5,4,'S_{O2}','r_{AMX,max}');

%O2-FWAS
tnremove(xbp3,xbp4,5,0,'S_{O2}','f_{WAS}');

%RAMX-FWAS
tnremove(xbp5,xbp6,0,4,'r_{AMX,max}','f_{WAS}');

function [Neff] = tnremove(x,y,flag,flag2,xparam,yparam)


for k = [1,2]
    if k == 1
        z = x;
    else
        z = y;
    end
    
    NH4 = ones(1,length(z(1,:)))*2;
    NO2 = z(1,:);
    
    if flag == 0
        SO2 = 0.3;
    else
        SO2 = z(flag,:);
    end
    
    if flag2 == 0
        p3 = 86/12.5;
    else
        p3 = z(flag2,:)/12.5;
    end
    
    NO3a = 12.5*0.34*(NO2./(0.5+NO2))*(SO2/(0.4+SO2));
    NO3b = 0.8772*p3.*(NH4./(0.03+NH4)).*(NO2./(0.005+NO2));
    
    NO3 = NO3a+NO3b;
    N2 = (2/0.17)*p3.*(NH4./(0.03+NH4)).*(NO2./(0.005+NO2));
    
    NH4_rem = 18;
    TN = NO2+NO3;
    
    Neff = 100*(NH4_rem-TN)/NH4_rem;
    
    %Plot
    par1 = z(4,:);
    par2 = z(5,:);
    colormap jet
    surface([par1;par1],[par2;par2],[Neff;Neff],[Neff;Neff],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',3);
    
    hold on
end
grid on
hcb=colorbar;
title(hcb,'N-removal (%)','fontsize',9)
xlabel(xparam,'fontsize',14)
ylabel(yparam,'fontsize',14)
