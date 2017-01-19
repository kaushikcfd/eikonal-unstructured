% This program gives the FIGURE 5 (a) and (b) of the 
% JCP revised version.

clear all;
close all;
npoints=1000;
vx=0.0;vy=0.6;

xmin=-0.0;xmax=1.1;
ymin=-0.5;ymax=1.5;

display('----------------------------------------------------------');

T=load('T1000.dat');

 x=[xmin:(xmax-xmin)/npoints:xmax];
 y=[ymin:(ymax-ymin)/npoints:ymax];

 hold off;
 %------------------------Before Caustic Region----------------------------
 wf_at_time=0.35;%Figure 5(a)
 figure;
 CFMMd=contour(x,y,T,wf_at_time,'b');
 M=size(CFMMd,2);
 CFMMd=CFMMd(:,[M:-1:2]);
 j=0;
 for i=1:M-1
     if(CFMMd(2,i)>-0.2 && CFMMd(2,i)<1.3)
         j=j+1;
         CFMM(:,j)=CFMMd(:,i);
     end
 end
 M=j;
%-----------------------Generating ray theory data---------------
 %Ray_Acoustics(xmin,xmax,ymin,ymax,vx,vy,wf_at_time,1,40000,0);
%----------------------------------------------------------------
 WFd=load('Linear_Wavefront035.dat');
 N=size(WFd);
 j=0;
 for i=1:N
     if(WFd(i,2)>-0.2 && WFd(i,2)<1.3)
         j=j+1;
         WF(j,:)=WFd(i,:);
     end
 end
 N=j;
 for i = 1:M
     for j=1:N
        D(j)=sqrt((CFMM(1,i)-WF(j,1))^2 + (CFMM(2,i)-WF(j,2))^2);
     end
     EN(i)=min(D);
 end
 fprintf('At T=0.35, Max(EN) = %e and Min(EN) = %e\n',max(EN), min(EN));
 
 figure;plot(CFMM(2,:),EN,'k','Linewidth',1.0);
 set(gca,'XLim',[-0.2 1.3]);
 xlabel('$y$','Interpreter','LaTex','fontsize',16,'fontname','Times','fontangle','italic');
 ylabel('$$|{\mbox{\boldmath{$e$}}}_j|$$','Interpreter','LaTex','fontsize',16,'fontname','Times','fontangle','italic');
 text(-0.08,0.00083,'\bf (a)','Fontsize',14);
 %------------------------After Caustic Region----------------------------
 wf_at_time=0.5;%Figure 5(a)
 figure;
 CFMMd=contour(x,y,T,wf_at_time,'b');
 M=size(CFMMd,2);
 CFMMd=CFMMd(:,[M:-1:2]);
 j=0;
 for i=1:M-1
     if(CFMMd(2,i)>-0.2 && CFMMd(2,i)<1.3)
         j=j+1;
         CFMM(:,j)=CFMMd(:,i);
     end
 end
 M=j;
%-----------------------Generating ray theory data---------------
 %Ray_Acoustics(xmin,xmax,ymin,ymax,vx,vy,wf_at_time,1,40000,0);
%----------------------------------------------------------------
 WFd=load('Linear_Wavefront05.dat');
 N=size(WFd);
 j=0;
 for i=1:N
     if(WFd(i,2)>-0.2 && WFd(i,2)<1.3)
         j=j+1;
         WF(j,:)=WFd(i,:);
     end
 end
 N=j;
 for i = 1:M
     for j=1:N
        D(j)=sqrt((CFMM(1,i)-WF(j,1))^2 + (CFMM(2,i)-WF(j,2))^2);
     end
     EN(i)=min(D);
 end
fprintf('At T=0.5, Max(EN) = %e and Min(EN) = %e\n',max(EN), min(EN));
 figure;plot(CFMM(2,:),EN,'k','Linewidth',1.0);
 set(gca,'XLim',[-0.2 1.3]);
 xlabel('$y$','Interpreter','LaTex','fontsize',16,'fontname','Times','fontangle','italic');
 ylabel('$$|{\mbox{\boldmath{$e$}}}_j|$$','Interpreter','LaTex','fontsize',16,'fontname','Times','fontangle','italic');
 text(-0.08,0.00063,'\bf (b)','Fontsize',14);
 
 display('-------------------End of the program-----------------------');