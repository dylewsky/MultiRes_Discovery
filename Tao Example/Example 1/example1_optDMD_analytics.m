load 'raw_data.mat';
load 'hf_sindy_data.mat';
nSteps = 2^14;

x_trueVars = ab_to_xy(x);
figure
plot(TimeSpan(1:nSteps),x_trueVars(:,1:nSteps),'LineWidth',2)
title('Projection into x,y variables')
legend('x1','x2','y1','y2')
ylim([-10 10])

disp(ab_to_xy(w_sorted_comb)/norm(ab_to_xy(w_sorted_comb)))

disp(ab_to_xy(w_sorted))

% GT:
% y1 = r*sin(th)
% y2 = -a + b^3
% 
% y1_t = r_t*sin(th) + r*th_t*cos(th)
% y2_t = -a_t + 3*b^2*b_t

