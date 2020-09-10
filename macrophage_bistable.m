
function macrophage_bistable()
clear all
    s1    = 3.75;
    s2    = 3.75;
    %other parameters, 
    a1    = 5;
    a2    = 5;
    k1    = 1;
    k2    = 1;
    p1    = 0.5;
    p2    = 1; 
    n     = 6;
    m     = 6;
    b1    = 0.05;
    b2    = 0.05;
    q1    = 5; %mu
    q2    = 5;
    

tmax=40;
tspan=[0,tmax];
y0=[1.2,2];
[T,Y]=ode45(@macrophageODE, tspan,y0);
FP=Y(end,:);
figure(1)
plot(T,Y,'linewidth', 1.5)
xlim([0,tmax])
ylim([0,2])
legend('x_1','x_2','fontsize', 16)
title('Bistability: Switch Behavior','fontsize', 16)


xlabel('t', 'fontsize', 16) ;
% Computation of the nullclines R(I) and I(R) and,
%   additionally, the directional field of the ODE
y1 = linspace(0.1,2,20) ;
y2 = linspace(0.1,2,20) ;
[Y1,Y2] = meshgrid(y1,y2) ;

u = zeros(size(Y1)) ;
v = zeros(size(Y1)) ;

for j = 1:numel(Y1)
	z = macrophageODE(0,[Y1(j) ; Y2(j)]) ;
	u(j) = z(1) ;
	v(j) = z(2) ;
end

figure(2) ;
quiver( Y1, Y2, u, v, 2, 'r', 'linewidth', 1.5 ) ;
hold('on') ;

nullclineU = p2*(a1*y1.^n./(k1^n+y1.^n)+s1)./(q1*y1-b1)-p2;
nullclineV = p1*s2./(q2*y2-b2-a2*y2.^m./(y2.^m+k2^m))-p1;

% Graphical output starts here

plot(y1, nullclineU, 'b-', 'linewidth', 2)
plot(nullclineV, y2, 'g-', 'linewidth', 2)
plot(Y(:,1),Y(:,2),'k','linewidth',1.5)
xlim([0,2])
ylim([0,2])
hold('off');
xlabel('x_1', 'fontsize', 16) ;
ylabel('x_2', 'fontsize', 16) ;

grid;
legend('Direction field', 'Nullcline for x_1', 'Nullcline for x_2', ...
       'location', 'northwest','fontsize', 12) ;
title('Bistability: Solution Trajectory','fontsize', 16)
   
end


function dy=macrophageODE(~,y)

dy=zeros(2,1);

    s1    = 3.75;
    s2    = 3.75;
    %other parameters, 
    a1    = 5;
    a2    = 5;
    k1    = 1;
    k2    = 1;
    p1    = 0.5;
    p2    = 1; 
    n     = 6;
    m     = 6;
    b1    = 0.05;
    b2    = 0.05;
    q1    = 5; %mu
    q2    = 5;

    u=y(1);
    v=y(2);
    dy(1)=(a1*u^n/(k1^n+u^n)+s1)*(1/(1+v/p2))+b1-q1*u;
    dy(2)=a2*v^m/(k2^m+v^m)+s2*(1/(1+u/p1))+b2-q2*v;  

end

