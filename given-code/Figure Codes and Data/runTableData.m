% vx=-0.0;vy=0.6; concave wavefront (focusing)
clear all


xmin=-0.0;xmax=1.1;
ymin=-0.5;ymax=1.5;



%-------------------------Reading data------------------------------------- 
% TN=[10:10:100];
% for i=1:10
%     str=num2str(TN(i));
%     str=strcat('T',str);
%     str=strcat(str,'.dat');
%     Tdummy=load(str);
%     
% end
% T(1,:,:)=Tdummy;
% return
%-------------------------Before focusing----------------------------------
npoints=200;
x=[xmin:(xmax-xmin)/npoints:xmax];
y=[ymin:(ymax-ymin)/npoints:ymax];

str=num2str(npoints);
str=strcat('T',str);
str=strcat(str,'.dat');
T=load(str);
str
wf_at_time=0.35;

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
 plot(CFMM(2,:),EN,'k','Linewidth',1.0);
 L1norm=0;L2norm=0;
 for i=1:M
     L1norm=L1norm+EN(i);
     L2norm=L2norm+EN(i)^2;
 end
 Dely=(ymax-ymin)/npoints;
 L2norm=sqrt(L2norm*Dely);
 L1norm=Dely*L1norm;
 Linfnorm=max(EN);
 fprintf('%d\t%f\t%f\t%f\n',npoints,L1norm,L2norm,Linfnorm);
 
%-------------------------After focusing-----------------------------------------

