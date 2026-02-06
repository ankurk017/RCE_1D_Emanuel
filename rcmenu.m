file=input('Enter subdirectory name (default = output)','s');
if isempty(file)
   file='output';
end

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
it=input('Menu choices: \n \n Evolution: \n   1) Precipitation and Evaporation \n   2) Surface air and sea surface temperatures \n   3) 500 mb T and q \n   4) TOA shortwave and longwave \n   5) Potential intensity   \n  \n  Vertical Profiles:  \n   6) Temperature \n   7) Specific humidity \n   8) Buoyancy \n   9) Relative humidity \n  10) Heating rates  \n  11) Mass fluxes \n  12) Entrainment and detrainment  \n  13) Moist static energy  \n  14) Cloud fraction  \n  15) Cloud condensed water  \n \n Time-height sections: \n  16) Buoyancy \n  17) Cloud fraction \n  18) T diff \n  19) Q diff \n  20) RH diff \n  21) Upward mass flux \n  22) Penetrative mass flux \n  23) Unsaturated downdraft mass flux \n  24) Net mass flux     \n \n  0) exit \n ');
if it == 1
   precip=time(:,2);
   evap=time(:,3);
   h=plot(timem,precip,'-',timem,evap,'--');
   set(h,'linewidth',2)
   set(gca,'fontweight','bold')
   xlabel('Time (days)','fontweight','bold')         
   ylabel('Precipitation and Evaporation (mm/day)','fontweight','bold')
   legend('Precipitation','Evaporation',0)
elseif it == 2
   sst=time(:,5);
   ta=time(:,4);
   h=plot(time1,sst,'-',time1,ta,'--');
   set(gca,'fontweight','bold')
   set(h,'linewidth',2)
   xlabel('Time (days)','fontweight','bold')         
   ylabel('SST and Ta (C)','fontweight','bold')
   legend('SST','Air temperature',0)
elseif it == 3
   t500=time(:,6);
   q500=time(:,7);
   h=plot(time1,t500,'-',time1,q500,'--');
   set(h,'linewidth',2)
   set(gca,'fontweight','bold')
   xlabel('Time (days)','fontweight','bold')         
   ylabel('500 mb T (C) and q (g/Kg)','fontweight','bold')
   legend('T','q',0)
elseif it == 4
   sw=time(:,8);
   lw=time(:,9);
   h=plot(time1,sw,'-',time1,lw,'--');
   set(h,'linewidth',2)
   set(gca,'fontweight','bold')
   xlabel('Time (days)','fontweight','bold')         
   ylabel('TOA Shortwave and Longwave (W/m^2)','fontweight','bold')
   legend('SW','LW',0)
elseif it == 5
   vmax=time(:,10);
   h=plot(time1,vmax);
   set(h,'linewidth',2)
   set(gca,'fontweight','bold')
   xlabel('Time (days)','fontweight','bold')         
   ylabel('V_{max} (m/s)','fontweight','bold')
elseif it == 51
   pmin=time(:,11);
   h=plot(time1,pmin);
   set(h,'linewidth',2)
   set(gca,'fontweight','bold')
   xlabel('Time (days)','fontweight','bold')         
   ylabel('p_{min} (mb)','fontweight','bold')
elseif it == 6
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
elseif it == 7
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
elseif it == 8
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
elseif it ==9 
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
elseif it ==10 
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
   legend('Convective','Radiative','Turbulent', 'Large-scale adiabatic',0)
elseif it == 11
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
   legend('Net upward','Penetrative','Unsaturated',0)
elseif it == 12
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
   legend('Entrainment','Detrainment',0)
elseif it == 13
   mse=profile(:,15);
   h=plot(mse,p);
   set(h,'linewidth',2)
   set(gca,'fontweight','bold')
   u=get(gca,'YTick');
   u1=transpose(-u);
   v=cellstr(num2str(u1));
   set(gca, 'YTickLabel',v)
   xlabel('Moist static energy (K)','fontweight','bold')         
   ylabel('Pressure (mb)','fontweight','bold')
elseif it == 14
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
elseif it == 15
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
elseif it == 16
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
   pause
   tht1=cfrac(3:35,:);
   thtm=mean(max(tht1,0));
   h=plot(timem,thtm);
   set(h,'linewidth',3)
   set(gca,'fontweight','bold')
   xlabel('Time (days)','fontweight','bold')         
   ylabel('Vertically averaged buoyancy','fontweight','bold')
elseif it == 17
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
   pause
   tht1=cfrac(3:35,:);
   thtm=mean(tht1);
   h=plot(timem,thtm);
   set(h,'linewidth',3)
   set(gca,'fontweight','bold')
   xlabel('Time (days)','fontweight','bold')         
   ylabel('Vertically averaged cloud fraction','fontweight','bold')
elseif it == 18
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
   pause
   tht1=tht(3:35,:);
   thtm=mean(tht1);
   h=plot(timem,thtm);
   set(h,'linewidth',3)
   set(gca,'fontweight','bold')
   xlabel('Time (days)','fontweight','bold')         
   ylabel('Vertically averaged temperature perturbation from initial state (K)','fontweight','bold')
elseif it == 19
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
   pause
   tht1=tht(3:35,:);
   thtm=mean(tht1);
   h=plot(timem,thtm);
   set(h,'linewidth',3)
   set(gca,'fontweight','bold')
   xlabel('Time (days)','fontweight','bold')         
   ylabel('Vertically averaged specific humidity perturbation from initial state (g/Kg)','fontweight','bold')
elseif it == 20
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
   pause
   tht1=tht(3:35,:);
   thtm=mean(tht1);
   h=plot(timem,thtm);
   set(h,'linewidth',3)
   set(gca,'fontweight','bold')
   xlabel('Time (days)','fontweight','bold')         
   ylabel('Vertically averaged relative humidity perturbation from initial state (%)','fontweight','bold')
elseif it == 21
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
   pause
   tht1=tht(3:35,:);
   thtm=mean(tht1);
   h=plot(timem,thtm);
   set(h,'linewidth',3)
   set(gca,'fontweight','bold')
   xlabel('Time (days)','fontweight','bold')         
   ylabel('Vertically averaged upward convective mass flux','fontweight','bold')   
elseif it == 22
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
   pause
   tht1=tht(3:35,:);
   thtm=mean(tht1);
   h=plot(timem,thtm);
   set(h,'linewidth',3)
   set(gca,'fontweight','bold')
   xlabel('Time (days)','fontweight','bold')         
   ylabel('Vertically averaged penetrative downdraft mass flux','fontweight','bold')   
%  % 
elseif it == 23
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
   pause
   tht1=tht(3:35,:);
   thtm=mean(tht1);
   h=plot(timem,thtm);
   set(h,'linewidth',3)
   set(gca,'fontweight','bold')
   xlabel('Time (days)','fontweight','bold')         
   ylabel('Vertically averaged unsaturated downdraft mass flux','fontweight','bold')   
elseif it == 24
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
   pause
   tht1=tht(3:35,:);
   thtm=mean(tht1);
   h=plot(timem,thtm);
   set(h,'linewidth',3)
   set(gca,'fontweight','bold')
   xlabel('Time (days)','fontweight','bold')         
   ylabel('Vertically averaged net convective upward mass flux','fontweight','bold')   
elseif it == 0
break
end
end

