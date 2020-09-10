%% To plot Basin of Attraction for a given set of equations using ode45 solver
%% Bistable Case

clc; clear all
warning('off') 

% The roots of the given governing equations 
 r1 = [0.5952  ;  0.3545] ; %for bi-stability case
 r2 = [1.3152  ;  0.2166] ;

% Initial conditions between 0 and 2 for two variables 
y1 = linspace(0,2,100) ;
y2 = linspace(0,2,100) ;

% Initialize the required matrices
Yr1 = [] ; Yr2 = [] ; 
% Time Span
tmax=80;
tspan=[0,tmax];

tic 
for i = 1:length(y1)
     for j = 1:length(y2)
          y0 = [y1(i);y2(j)] ; %initial condition between 0 and 2 both for x1 and x2
          % Solve the system of Equations using ode45 solver
          
          [T,Y]=ode45(@macrophageODE, tspan,y0);
          FP=transpose(Y(end,:));

          % Locating the initial conditions according to error
          if norm(FP-r1)<1e-8 
             Yr1 = [y0 Yr1]  ;
          elseif norm(FP-r2)<1e-8
             Yr2 = [y0 Yr2] ;
          else
              x=[norm(FP-r1),norm(FP-r2)];
              I=find(x==min(x));
              if I==1
                  Yr1 = [y0 Yr1];
              elseif I==2
                  Yr2 = [y0 Yr2];
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

plot(Yr1(1,:),Yr1(2,:),'.','color','r') ;
plot(Yr2(1,:),Yr2(2,:),'.','color','b') ;
plot(r1(1), r1(2),'k*','MarkerSize',20) ;
plot(r2(1), r2(2),'kd','MarkerSize',10) ;

legend('(0.5952, 0.3545)','(1.3152, 0.2166)', 'fontsize', 12) ;
title('Basin of Attration for Bistability','fontsize', 16);
xlabel('x_1', 'fontsize', 16) ;
ylabel('x_2', 'fontsize', 16) ;

function dy=macrophageODE(~,y)

dy=zeros(2,1);

    a1    = 5;
    a2    = 5;
    k1    = 1;
    k2    = 1;
    p1    = 0.5;
    p2    = 1; 
    n1    = 6;
    n2    = 6;
    b1    = 0.05;
    b2    = 0.05;
    q1    = 5; 
    q2    = 5;
    
    s1    =3.75;
    s2    =3.75;

    u=y(1);
    v=y(2);
    dy(1)=(a1*u^(n1)/(k1^(n1)+u^(n1))+s1)*(1/(1+v/p2))+b1-q1*u;
    dy(2)=a2*v^(n2)/(k2^(n2)+v^(n2))+s2*(1/(1+u/p1))+b2-q2*v;  
   
end

