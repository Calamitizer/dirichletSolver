% dirichletSolverScript.m
% J. Alex Ruble
% MA 428 Project 1

% Hardcode problem values
% region, domainDefFunc, nonhomFunc, res
region = {[0,5],[0,4]};
h = @(x) heaviside(x)+(x == 0)/2; % heaviside with h(0)==1
yLower = @(x) 0;
yUpper = @(x) h(-x+2) + h(-x+3) + 2*h(-x+5);
domainDefFunc = @(x,y) all([region{1}(1)<=x, x<=region{1}(2), yLower(x)<=y, y<=yUpper(x)]);
nonhomFunc = @(x,y) x.*exp(-x.^2-y.^2);
% nonhomFunc = @(x,y) x*(x-y)^3;
funcName = char(nonhomFunc);
funcName = funcName(7:end);
res = 2;

[U,Boundary,xP,yP,A] = dirichletSolver(region,domainDefFunc,nonhomFunc,res);

figure(1);
hold on;
set(0,'DefaultTextInterpreter', 'none');
surf(xP,yP,U,'FaceColor','interp','FaceLighting','gouraud','FaceAlpha',.8,'EdgeColor',[.25 .25 .25],'EdgeAlpha',.2)
mesh(xP,yP,Boundary,'FaceColor','flat','FaceColor',[1 0 0],'FaceAlpha',.25,'EdgeColor',[.8 0 0],'EdgeAlpha',.75)
fs = 14;
title(['-L(u) = ',funcName,', h = 2^-',int2str(res)],'FontSize',fs);
xlabel('x','FontSize',fs);
ylabel('y','FontSize',fs);
zlabel('u(x,y)','FontSize',fs);
legend('u(x,y)','Boundary');

figure(2);
hold on;
spy(A>0,'r.',1);
spy(A<0,'b.',1);
title(['spy(A), h=2^-',int2str(res)]);