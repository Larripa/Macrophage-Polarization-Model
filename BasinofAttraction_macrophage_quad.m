%% To plot Basin of Attraction for a given set of equations using ode45 solver
%% Quadstable Case
clc; clear all
warning('off') 

% The roots of the given governing equations (from quad-stability case)
r1 = [0.3066 ;1.8930] ;
r2 = [0.6259 ;0.3969] ;
r3 = [1.3820 ;1.5095] ;
r4 = [3.0609 ;0.1296];

% Initial conditions between 0 and 2 for two variables 
y1 = linspace(0,3.5,100) ;
y2 = linspace(0,2.5,100) ;
% Initialize the required matrices
Yr1 = [] ; Yr2 = [] ; Yr3 = [] ; Yr4 = [] ;
% Time Span
tmax=200;
tspan=[0,tmax];

tic 
for i = 1:length(y1)
     for j = 1:length(y2)
          y0 = [y1(i);y2(j)] ; %initial condition between 0 and 2 both for x1 and x2
          % Solve the system of Equations using ode45 solver
          
          [T,Y]=ode45(@macrophageODE, tspan,y0);
          FP=transpose(Y(end,:));

          % Locating the initial conditions according to error
          if norm(FP-r1)<1e-5
             Yr1 = [y0 Yr1]  ;
          elseif norm(FP-r2)<1e-5
             Yr2 = [y0 Yr2] ;
          elseif norm(FP-r3)<1e-5
             Yr3 = [y0 Yr3] ;
          elseif norm(FP-r4)<1e-5 % if not close to any of the roots
             Yr4 = [y0 Yr4] ;
          else
              x=[norm(FP-r1),norm(FP-r2),norm(FP-r3),norm(FP-r4)];
              I=find(x==min(x));
              if I==1
                  Yr1 = [y0 Yr1];
              elseif I==2
                  Yr2 = [y0 Yr2];
              elseif I==3
                  Yr3 = [y0 Yr3];
              elseif I==4
                  Yr4 = [y0 Yr4];
              end
          end
          
     end 
end
toc
warning('on') % Remove the warning off constraint
% Initialize figure
figure
%set(gcf,'color','w') 
hold on

plot(Yr1(1,:),Yr1(2,:),'.','color','g') ;
plot(Yr2(1,:),Yr2(2,:),'.','color','r') ;
plot(Yr3(1,:),Yr3(2,:),'.','color','m') ;
plot(Yr4(1,:),Yr4(2,:),'.','color','b') ;

plot(r1(1), r1(2),'k.','MarkerSize',20) ;
plot(r2(1), r2(2),'k*','MarkerSize',10) ;
plot(r3(1), r3(2),'ko','MarkerSize',10) ;
plot(r4(1), r4(2),'kd','MarkerSize',10) ;

legend('(0.3066, 1.8930)','(0.6259, 0.3969)','(1.3820, 1.5095)','(3.0609 ;0.1296)', 'fontsize', 12) ;
title('Basin of Attration for Quadstability','fontsize', 16);
xlabel('x_1', 'fontsize', 16) ;
ylabel('x_2', 'fontsize', 16) ;
xlim([0,3.5])
ylim([0,2.5])
function dy=macrophageODE(~,y)

dy=zeros(2,1);

    a1    = 15;
    a2    = 8;
    k1    = 1;
    k2    = 1;
    p1    = 0.5;
    p2    = 1; 
    n1    = 22;
    n2    = 6;
    b1    = 0.05;
    b2    = 0.05;
    q1    = 5.8; 
    q2    = 5.8;
    
    s1    = 5;
    s2    = 5;
    
    u=y(1);
    v=y(2);
    dy(1)=(a1*u^(n1)/(k1^(n1)+u^(n1))+s1)*(1/(1+v/p2))+b1-q1*u;
    dy(2)=a2*v^(n2)/(k2^(n2)+v^(n2))+s2*(1/(1+u/p1))+b2-q2*v;  
end

