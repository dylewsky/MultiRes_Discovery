clc; clear; close all;

%% parameters

a0=0;   b0=0;   r0=sqrt(2^2+0.8^2); theta0=atan2(0.8,2);

epsilon=0.01;
delta = 4;

T=24;



%% RK4 integration of the mixed system

h=epsilon/100;

TimeSpan=[0:h:T];   TimeSteps=length(TimeSpan)-1;

x=zeros(4,TimeSteps+1); x(:,1)=[a0; b0; r0; theta0];

for i=1:TimeSteps

    xx=x(:,i);

    a=xx(1); b=xx(2); r=xx(3); theta=xx(4);

    k1=[(3*r*b^2*cos(theta)+r*sin(theta)/epsilon + delta*(r*sin(theta))^3)/(3*b^2+1); ...

        (r*cos(theta)-r*sin(theta)/epsilon - delta*(r*sin(theta))^3)/(1+3*b^2); ...

        sin(theta)*(-a+b^3-(a+b)^3*r^2*cos(theta)*sin(theta)); ...

        ((-a+b^3)*cos(theta)+(a+b)^3*r^2*sin(theta)^3)/r];

    

    xx=x(:,i)+h/2*k1;

    a=xx(1); b=xx(2); r=xx(3); theta=xx(4);

    k2=[(3*r*b^2*cos(theta)+r*sin(theta)/epsilon + delta*(r*sin(theta))^3)/(3*b^2+1); ...

        (r*cos(theta)-r*sin(theta)/epsilon - delta*(r*sin(theta))^3)/(1+3*b^2); ...

        sin(theta)*(-a+b^3-(a+b)^3*r^2*cos(theta)*sin(theta)); ...

        ((-a+b^3)*cos(theta)+(a+b)^3*r^2*sin(theta)^3)/r];

    

    xx=x(:,i)+h/2*k2;

    a=xx(1); b=xx(2); r=xx(3); theta=xx(4);

    k3=[(3*r*b^2*cos(theta)+r*sin(theta)/epsilon + delta*(r*sin(theta))^3)/(3*b^2+1); ...

        (r*cos(theta)-r*sin(theta)/epsilon - delta*(r*sin(theta))^3)/(1+3*b^2); ...

        sin(theta)*(-a+b^3-(a+b)^3*r^2*cos(theta)*sin(theta)); ...

        ((-a+b^3)*cos(theta)+(a+b)^3*r^2*sin(theta)^3)/r];

    

    xx=x(:,i)+h*k3;

    a=xx(1); b=xx(2); r=xx(3); theta=xx(4);

    k4=[(3*r*b^2*cos(theta)+r*sin(theta)/epsilon + delta*(r*sin(theta))^3)/(3*b^2+1); ...

        (r*cos(theta)-r*sin(theta)/epsilon - delta*(r*sin(theta))^3)/(1+3*b^2); ...

        sin(theta)*(-a+b^3-(a+b)^3*r^2*cos(theta)*sin(theta)); ...

        ((-a+b^3)*cos(theta)+(a+b)^3*r^2*sin(theta)^3)/r];

    

    x(:,i+1)=x(:,i)+(k1+2*k2+2*k3+k4)/6*h;

end



%%

plot(TimeSpan, x);

title('Raw data');

xlabel('Time');



%%

figure

x1_=x(1,:)+x(2,:);

x2_=x(3,:).*cos(x(4,:));

y1_=x(3,:).*sin(x(4,:));

y2_=-x(1,:)+x(2,:).^3;

x_=x;



subplot(2,2,1);

plot(TimeSpan,x1_);

axes1=axis;

title('hidden slow x_1 (exact)');

subplot(2,2,3);

plot(TimeSpan,x2_);

axes3=axis;

title('hidden slow x_2 (exact)');



x1IC=a0+b0; x2IC=r0*cos(theta0);    y1IC=r0*sin(theta0);    y2IC=-a0+b0^3;

avgCoeff=(y1IC^2+epsilon*y2IC^2)/2;

X=zeros(2,TimeSteps+1);    X(:,1)=[x1IC; x2IC];

for i=1:TimeSteps

    XX=X(:,i);

    k1=[XX(2); -avgCoeff*XX(1)^3];

    

    XX=X(:,i)+h/2*k1;

    k2=[XX(2); -avgCoeff*XX(1)^3];



    XX=X(:,i)+h/2*k2;

    k3=[XX(2); -avgCoeff*XX(1)^3];



    XX=X(:,i)+h*k3;

    k4=[XX(2); -avgCoeff*XX(1)^3];

    

    X(:,i+1)=X(:,i)+(k1+2*k2+2*k3+k4)/6*h;

end

subplot(2,2,2);

plot(TimeSpan,X(1,:));

axis(axes1);

title('hidden slow X_1 (approximate, closed)');

subplot(2,2,4);

plot(TimeSpan,X(2,:));

axis(axes3);

title('hidden slow X_2 (approximate, closed)');

%% 
save('raw_data_2_hiRes.mat','x','TimeSpan')

%%
figure
plot(TimeSpan,x,'b', 'LineWidth', 1.5)
hold on
load('raw_data_2_hiRes.mat')
plot(TimeSpan,x,'k', 'LineWidth', 1.5)
xlim([TimeSpan(1) TimeSpan(end)])
title('Linear (Black) vs. Nonlinear (Blue)')
