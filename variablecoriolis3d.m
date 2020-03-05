close all; clear all; clc;
T=50000;

dx = 2000;
dy = 2000;
D=10;
g=9.8;
n=0;
dt=100;
eps = 0.1;
z=0;
x=-100000:dx:(100000-dx);
y=-100000:dy:(100000-dy);
u=zeros(length(x),length(y));
v=zeros(length(x),length(y));
u=zeros(length(x),length(y));
v=zeros(length(x),length(y));
uplus=u;
utot=u;
uminus=u;
vplus=v;
vtot=v;
vminus=v;
et=zeros(length(x),length(y));
etanew=et;
etaa=et; %etaaa
h=D+et;

fid=figure;
 
eta=et;
for k=2:length(et)-1
       for l=2:length(et)-1
           etanew(l,k)=eta(l,k)+0.5*exp((-(k-50)^2/(2*10^2))-(l-99)^2/(2*10^2)); %gør en gausisk
       end
end
f=5e-4;

  et=etanew;
   eta=et;
     us=u;
    vs=v;
    
    beta2=2e-11;
    fr(1:length(y))=f+beta2.*y;
    alfa(1:length(fr))=dt.*fr;
beta(1:length(alfa))=1/4 * alfa.^2;
    
    
figure(fid)
    hold on
    writerObj = VideoWriter('simulation.avi');
writerObj.FrameRate=60;
open(writerObj);
while n<T
  
   for k=2:length(et)-1
       for l=2:length(et)-1
       u(k,l) = us(k,l)-g*dt/dx*(et(k,l+1)- et(k,l)); % finder hastigheder
u(:,end)=0; %sætter sidste hastighed =0
u(:,end-1)=0;
u(:,1)=0;

%uc(k,l)=us(k,l)+1/2*dt*f*(vs(l,k)+vs(l+1,k)); %coriolis u

utot(k,l)=(u(k,l)-beta(k)*us(k,l)+alfa(k)*vs(k,l))/(1+beta(k));


uplus(k,l)=(utot(k,l)+abs(utot(k,l))) * 1/2; %finder uplus
uminus(k,l)=(utot(k,l)-abs(utot(k,l)))* 1/2; % finder uminus
       end
       for l=2:length(et)-1
       v(l,k) = vs(l,k)-g*dt/dy*(et(l+1,k)- et(l,k)); % finder hastigheder
v(end,:)=0; %sætter sidste hastighed =0
v(end-1,:)=0;
v(1,:)=0;

vtot(l,k)=(v(l,k)-beta(l)*vs(l,k)-alfa(l)*us(l,k))/(1+beta(l));

vplus(l,k)=(vtot(l,k)+abs(vtot(l,k))) * 1/2; %finder vplus
vminus(l,k)=(vtot(l,k)-abs(vtot(l,k)))* 1/2; % finder vminus
       end
       
    
   end
  
   u=utot;
       v=vtot;
         us=u;
    vs=v;
    eta(:,:)=et(:,:);
    for k=2:length(et)-1
        for l=2:length(et)-1
            h(k,l)=D+eta(k,l); %udregner højden
        end
    end
    
   for k =2:length(et)-1
      for l=2:length(et)-1
     
    
   et(k,l)=eta(k,l)-dt/dx*(uplus(k,l)*h(k,l)+uminus(k,l)*h(k,l+1)-uplus(k,l-1)*h(k,l-1)-uminus(k,l-1)*h(k,l) ) ...
       -dt/dy*(vplus(k,l)*h(k,l)+vminus(k,l)*h(k+1,l)-vplus(k-1,l)*h(k-1,l)-vminus(k-1,l)*h(k,l) ); %finder den nye eta i alle søjler 

   
      end
   end
   
   for k=2:length(et)-1
       for l=2:length(et)-1
   etaa(k,l)=(1-eps)*et(k,l) + 1/4 * eps*(et(k,l-1)+et(k,l+1))+ 1/4 * eps*(et(k-1,l)+et(k+1,l)); %ligger filter på
   end
   end
   et=etaa(1:length(et),1:length(et));
   et(end,:)=et(end-1,:);
et(:,end)=et(:,end-1);
%et(1,:)=et(2,:);
%et(:,1)=et(:,2);
%contourf(x(1:length(eta)),y(1:length(eta)),eta)
surf(x(1:length(eta)),y(1:length(eta)),eta)
  axis([min(x) max(x) min(y) max(y) -0.5 0.5]);
  colorbar
   frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
   shading interp
   pause(0.0001)
    hold off
    
    %z=z+1;
   %etasave(z,:)=et(50,:);
    n=n+dt; %finder ny tid
end
writeVideo(writerObj, frame);
close(writerObj); % Saves the movie.