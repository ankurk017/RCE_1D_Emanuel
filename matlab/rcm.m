%{
file=input('Enter subdirectory name (default = output)','s');
if isempty(file)
   file='output';
end
%}
file = 'output';
filesep='/';
load(strcat(file,filesep,'time.out'))
load(strcat(file,filesep,'profile.out'))
load(strcat(file,filesep,'cldhov.out'))
load(strcat(file,filesep,'thov.out'))
load(strcat(file,filesep,'qhov.out'))
load(strcat(file,filesep,'rhhov.out'))
load(strcat(file,filesep,'mhov.out'))
load(strcat(file,filesep,'mdhov.out'))
load(strcat(file,filesep,'mphov.out'))
load(strcat(file,filesep,'buoy.out'))
load(strcat(file,filesep,'omhov.out'))
p=-profile(:,1);
time1=time(:,1);
timem=time(:,12);
tout=mean(time(:,13));
tk=profile(:,2)+273.15;
qm=0.001.*profile(:,3);
tv=tk.*(1+0.608.*qm);
aa=size(p);
a=aa(1,1);
z(1)=0.0;
for i=2:a,
    z(i)=z(i-1)+(287./9.81).*0.5.*(tv(i)+tv(i-1)).*log(p(i-1)./p(i));
end
zkm=0.001.*z;
i=1;
while i == 1
    clc
    it=menu('Select from the following:','Precipitation and Evaporation (evolution)',...
        'Surface air and sea surface temperatures (evolution)',...
        '500 hPa Temperature and specific humidity (evolution)',...
        'Top of atmos. shortwave and longwave fluxes (evolution)',...
        'Temperature (profile)',...
        'Specific humidity (profile)',...
        'Buoyancy of parcel lifted from lowest layer (profile)',...
        'Relative humidity (profile)',...
        'Heating rates (profile)',...
        'Convective mass fluxes (profile)',...
        'Convective entrainment and detrainment (profile)',...
        'Moist static energy (profile)',...
        'Cloud fraction (profile)',...
        'Cloud water concentration (profile)',...
        'Buoyancy of parcel lifted from lowest layer (Time-height)',...
        'Cloud fraction (time-height)',...
        'Temperature difference from initial profile (time-height)',...
        'Specific humidity difference from initial profile (time-height)',...
        'Relative humidity difference from initial profile (time-height)',...
        'Updraft convective mass flux (time-height)',...
        'Penetrative convective downdraft mass flux (time-height)',...
        'Unsaturated convective downdraft mass flux (time-height)',...
        'Net convective upward mass flux (time-height)',...
        'Omega (time-height)',...
        'Exit');
    if isempty(it)
        
        pause(0.05); % This is a total hack to suppress some weird behavior where it interprets blank input as 1.
    elseif ~isempty(it)
        
        figure(1)
        if it == 1
            precip=time(:,2);
            evap=time(:,3);
            h=plot(timem,precip,'-',timem,evap,'--');
            set(h,'linewidth',2)
            set(gca,'fontweight','bold')
            xlabel('Time (days)','fontweight','bold')
            ylabel('Precipitation and Evaporation (mm/day)','fontweight','bold')
            title('x variable=timem, yvariables=precip, evap','fontweight','bold')
            legend('Precipitation','Evaporation')
        elseif it == 2
            sst=time(:,5);
            ta=time(:,4);
            h=plot(time1,sst,'-',time1,ta,'--');
            set(gca,'fontweight','bold')
            set(h,'linewidth',2)
            xlabel('Time (days)','fontweight','bold')
            ylabel('SST and Ta (C)','fontweight','bold')
            title('x variable=time1, yvariables=sst, ta','fontweight','bold')
            legend('SST','Air temperature')
        elseif it == 3
            t500=time(:,6);
            q500=time(:,7);
            h=plot(time1,t500,'-',time1,q500,'--');
            set(h,'linewidth',2)
            set(gca,'fontweight','bold')
            xlabel('Time (days)','fontweight','bold')
            ylabel('500 mb T (C) and q (g/Kg)','fontweight','bold')
            title('x variable=time1, yvariables=t500, q500','fontweight','bold')
            legend('T','q')
        elseif it == 4
            sw=time(:,8);
            lw=time(:,9);
            h=plot(time1,sw,'-',time1,lw,'--');
            set(h,'linewidth',2)
            set(gca,'fontweight','bold')
            xlabel('Time (days)','fontweight','bold')
            ylabel('TOA Shortwave and Longwave (W/m^2)','fontweight','bold')
            title('x variable=time1, yvariables=sw, lw','fontweight','bold')
            legend('SW','LW')
        elseif it == 5
            tpro=profile(:,2);
            h=plot(tpro,p);
            set(h,'linewidth',2)
            set(gca,'fontweight','bold')
            u=get(gca,'YTick');
            u1=transpose(-u);
            v=cellstr(num2str(u1));
            set(gca, 'YTickLabel',v)
            xlabel('Temperature (C)','fontweight','bold')
            ylabel('Pressure (mb)','fontweight','bold')
            title('x variable=tpro, yvariable=p','fontweight','bold')
        elseif it == 6
            qpro=profile(:,3);
            h=plot(qpro,p);
            set(h,'linewidth',2)
            set(gca,'fontweight','bold')
            u=get(gca,'YTick');
            u1=transpose(-u);
            v=cellstr(num2str(u1));
            set(gca, 'YTickLabel',v)
            xlabel('Specific humidity (gm/Kg)','fontweight','bold')
            ylabel('Pressure (mb)','fontweight','bold')
            title('x variable=qpro, yvariable=p','fontweight','bold')
        elseif it == 7
            buoy1=profile(:,4);
            h=plot(buoy1,p);
            set(h,'linewidth',2)
            set(gca,'fontweight','bold')
            u=get(gca,'YTick');
            u1=transpose(-u);
            v=cellstr(num2str(u1));
            set(gca, 'YTickLabel',v)
            xlabel('Buoyancy (C)','fontweight','bold')
            ylabel('Pressure (mb)','fontweight','bold')
            title('x variable=buoy1, yvariable=p','fontweight','bold')
        elseif it ==8
            rh=100.*profile(:,5);
            h=plot(rh,p);
            set(h,'linewidth',2)
            set(gca,'fontweight','bold')
            u=get(gca,'YTick');
            u1=transpose(-u);
            v=cellstr(num2str(u1));
            set(gca, 'YTickLabel',v)
            xlabel('Relative humidity (%)','fontweight','bold')
            ylabel('Pressure (mb)','fontweight','bold')
            title('x variable=rh, yvariable=p','fontweight','bold')
        elseif it ==9
            cq=profile(:,6);
            rq=profile(:,7);
            tq=profile(:,8);
            aq=profile(:,9);
            h=plot(cq,p,rq,p,tq,p,aq,p);
            set(h,'linewidth',2)
            set(gca,'fontweight','bold')
            u=get(gca,'YTick');
            u1=transpose(-u);
            v=cellstr(num2str(u1));
            set(gca, 'YTickLabel',v)
            xlabel('Heating (C/day)','fontweight','bold')
            ylabel('Pressure (mb)','fontweight','bold')
            title('x variables=cq, rq, tq, aq; yvariable=p','fontweight','bold')
            legend('Convective','Radiative','Turbulent', 'Large-scale adiabatic')
        elseif it == 10
            mu=profile(:,10);
            mpen=profile(:,11);
            md=profile(:,12);
            h=plot(mu,p,mpen,p,md,p);
            set(h,'linewidth',2)
            set(gca,'fontweight','bold')
            u=get(gca,'YTick');
            u1=transpose(-u);
            v=cellstr(num2str(u1));
            set(gca, 'YTickLabel',v)
            xlabel('Mass flux (10^{-3} Kg m^{-2} s^{-1})','fontweight','bold')
            ylabel('Pressure (mb)','fontweight','bold')
            title('x variables=mu, mpen, md; yvariable=p','fontweight','bold')
            legend('Net upward','Penetrative','Unsaturated')
        elseif it == 11
            ent=profile(:,13);
            det=profile(:,14);
            h=plot(ent,p,det,p);
            set(h,'linewidth',2)
            set(gca,'fontweight','bold')
            u=get(gca,'YTick');
            u1=transpose(-u);
            v=cellstr(num2str(u1));
            set(gca, 'YTickLabel',v)
            xlabel('Entrainment and detrainment (10^{-3} Kg m^{-2} s^{-1})','fontweight','bold')
            ylabel('Pressure (mb)','fontweight','bold')
            title('x variables=ent, det; yvariable=p','fontweight','bold')
            legend('Entrainment','Detrainment')
        elseif it == 12
            mse=profile(:,15);
            mse=min(mse,mse(1)+20);
            h=plot(mse,p);
            set(h,'linewidth',2)
            set(gca,'fontweight','bold')
            u=get(gca,'YTick');
            u1=transpose(-u);
            v=cellstr(num2str(u1));
            set(gca, 'YTickLabel',v)
            xlabel('Moist static energy (K)','fontweight','bold')
            ylabel('Pressure (mb)','fontweight','bold')
            title('x variable=mse, yvariable=p','fontweight','bold')
        elseif it == 13
            clf=profile(:,16);
            h=plot(clf,p);
            set(h,'linewidth',2)
            set(gca,'fontweight','bold')
            u=get(gca,'YTick');
            u1=transpose(-u);
            v=cellstr(num2str(u1));
            set(gca, 'YTickLabel',v)
            xlabel('Cloud fraction','fontweight','bold')
            ylabel('Pressure (mb)','fontweight','bold')
            title('x variable=clf, yvariable=p','fontweight','bold')
        elseif it == 14
            clq=profile(:,17);
            h=plot(clq,p);
            set(h,'linewidth',2)
            set(gca,'fontweight','bold')
            u=get(gca,'YTick');
            u1=transpose(-u);
            v=cellstr(num2str(u1));
            set(gca, 'YTickLabel',v)
            xlabel('Cloud condensed water (g/Kg)','fontweight','bold')
            ylabel('Pressure (mb)','fontweight','bold')
            title('x variable=clq, yvariable=p','fontweight','bold')
        elseif it == 15
            buoy=max(buoy,-0.1);
            [ntemp,mtemp]=size(buoy);
            if ntemp > mtemp
                buoy=transpose(buoy);
            end
            %cfrac=transpose(buoy);
            cfrac=buoy;
            contourf(timem,zkm,cfrac)
            %shading interp
            xlabel('Time (days)','fontweight','bold')
            %ylabel('Pressure (mb)','fontweight','bold')
            ylabel('Altitude (km)','fontweight','bold')
            title('Lifted parcel buoyancy (K)','fontweight','bold')
            set(gca,'fontweight','bold')
            u=get(gca,'YTick');
            u1=transpose(-u);
            v=cellstr(num2str(u1));
            %set(gca, 'YTickLabel',v)
            set(gca,'YLim',[0 16])
            colorbar
        elseif it == 16
            cfrac=transpose(cldhov);
            contourf(timem,zkm,cfrac)
            xlabel('Time (days)','fontweight','bold')
            %ylabel('Pressure (mb)','fontweight','bold')
            ylabel('Altitude (km)','fontweight','bold')
            title('Cloud Fraction','fontweight','bold')
            set(gca,'fontweight','bold')
            u=get(gca,'YTick');
            u1=transpose(-u);
            v=cellstr(num2str(u1));
            %set(gca, 'YTickLabel',v)
            set(gca,'YLim',[0 16])
            colorbar
        elseif it == 17
            tht=transpose(thov);
            [m,n]=size(tht);
            ttemp=zeros(m,n);
            for j=1:n,
                ttemp(:,j)=tht(:,j)-tht(:,1);
            end
            tht=ttemp;
            contourf(timem,p,tht)
            xlabel('Time (days)','fontweight','bold')
            ylabel('Pressure (mb)','fontweight','bold')
            title('Temperature perturbation from initial state (K)','fontweight','bold')
            set(gca,'fontweight','bold')
            u=get(gca,'YTick');
            u1=transpose(-u);
            v=cellstr(num2str(u1));
            set(gca, 'YTickLabel',v)
            colorbar
        elseif it == 18
            tht=transpose(qhov);
            [m,n]=size(tht);
            ttemp=zeros(m,n);
            for j=1:n,
                ttemp(:,j)=tht(:,j)-tht(:,n);
            end
            tht=ttemp;
            contourf(timem,p,tht)
            xlabel('Time (days)','fontweight','bold')
            ylabel('Pressure (mb)','fontweight','bold')
            title(' Specific humidity perturbation from initial state (g/Kg)','fontweight','bold')
            set(gca,'fontweight','bold')
            u=get(gca,'YTick');
            u1=transpose(-u);
            v=cellstr(num2str(u1));
            set(gca, 'YTickLabel',v)
            colorbar
        elseif it == 19
            tht=transpose(rhhov);
            [m,n]=size(tht);
            ttemp=zeros(m,n);
            for j=1:n,
                ttemp(:,j)=tht(:,j)-tht(:,n);
            end
            tht=ttemp;
            contourf(timem,p,tht)
            xlabel('Time (days)','fontweight','bold')
            ylabel('Pressure (mb)','fontweight','bold')
            title('Relative humidity perturbation from initial state (%)','fontweight','bold')
            set(gca,'fontweight','bold')
            u=get(gca,'YTick');
            u1=transpose(-u);
            v=cellstr(num2str(u1));
            set(gca, 'YTickLabel',v)
            colorbar
        elseif it == 20
            tht=transpose(mhov);
            contourf(timem,p,tht)
            xlabel('Time (days)','fontweight','bold')
            ylabel('Pressure (mb)','fontweight','bold')
            title('Upward convective upward mass flux (10^{-3} Kg m^{-2} s^{-1})','fontweight','bold')
            set(gca,'fontweight','bold')
            u=get(gca,'YTick');
            u1=transpose(-u);
            v=cellstr(num2str(u1));
            set(gca, 'YTickLabel',v)
            colorbar
        elseif it == 21
            tht=transpose(mdhov);
            contourf(timem,p,tht)
            xlabel('Time (days)','fontweight','bold')
            ylabel('Pressure (mb)','fontweight','bold')
            title('Penetrative downdraft mass flux (10^{-3} Kg m^{-2} s^{-1})','fontweight','bold')
            set(gca,'fontweight','bold')
            u=get(gca,'YTick');
            u1=transpose(-u);
            v=cellstr(num2str(u1));
            set(gca, 'YTickLabel',v)
            colorbar
        elseif it == 22
            tht=transpose(mphov);
            contourf(timem,p,tht)
            xlabel('Time (days)','fontweight','bold')
            ylabel('Pressure (mb)','fontweight','bold')
            title('Unsaturated downdraft mass flux (10^{-3} Kg m^{-2} s^{-1})','fontweight','bold')
            set(gca,'fontweight','bold')
            u=get(gca,'YTick');
            u1=transpose(-u);
            v=cellstr(num2str(u1));
            set(gca, 'YTickLabel',v)
            colorbar
        elseif it == 23
            tht=transpose((mhov+mdhov+mphov));
            contourf(timem,p,tht)
            xlabel('Time (days)','fontweight','bold')
            ylabel('Pressure (mb)','fontweight','bold')
            title('Net convective upward mass flux (10^{-3} Kg m^{-2} s^{-1})','fontweight','bold')
            set(gca,'fontweight','bold')
            u=get(gca,'YTick');
            u1=transpose(-u);
            v=cellstr(num2str(u1));
            set(gca, 'YTickLabel',v)
            colorbar
        elseif it == 24
            tht=transpose((omhov));
            contourf(timem,p,tht)
            xlabel('Time (days)','fontweight','bold')
            ylabel('Pressure (mb)','fontweight','bold')
            title('Omega (hPa/hr)','fontweight','bold')
            set(gca,'fontweight','bold')
            u=get(gca,'YTick');
            u1=transpose(-u);
            v=cellstr(num2str(u1));
            set(gca, 'YTickLabel',v)
            colorbar
        elseif it == 25
            i = 0;
        end
        figure(gcf)
        drawnow;
        pause(0.05);
    end
end

