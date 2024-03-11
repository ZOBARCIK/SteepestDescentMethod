clc;
clear;
clear all;

format long
%Ahmet Arda Zobar
%Steepest Descent Algorithm

%---------------------STEP 1-------------------------
%variables from the midterm assignment document.  
% w1 
% w2
%taken as x(1) and x(2)


fnc = @(x) ((x(1)^2 + x(2)^2 + 1.636*x(1)*x(2))./2)-0.8182*x(1)-0.354*x(2)+0.1250; %function itself with an handle

gx= @(x)  x(1)+(1.6365/2)*x(2) -0.8182; %x gradient (x partial derivative)
gy= @(x) x(2)+(1.6365/2)*x(1) -0.354; %y gradient (y partial derivative) 

x=[0 0];         

N=150; %(max iteration)

ep1=10^-10;
ep2=10^-10;
ep3=10^-10;

k=1; %iteration index

%---------------------------------STEP 2-------------------------------

while k<=N
    x_e = x;
    %first x is x_e
    
    %%%%%%%%%%%
    n=1
   grad= [gx(x); gy(x)]; %gradient vector {\nabla f(xk)}
   dir= grad; %direction vectoor
pk=dir;  %direction vector pk= {\nabla f(xk)}
npk=n*pk;
%%%%%%%

    x= x - npk;
    
    xdump(k) = x(1);
    ydump(k) = x(2);
    zdump(k) = fnc(x);
    gdump(k) = grad(1);
    
    
    k=k+1;
    df=fnc(x)-fnc(x_e); %delta f
    dx= x - x_e;
    if abs(df)<=ep1
        break
    elseif abs(norm(dx))<= ep2
        break
    elseif abs(norm(pk)) <= ep3
       
        break
    end

end
fprintf('delta_f \n')
disp(df)

fprintf('delta_x \n')
disp(dx)

disp('number of iterations have been done')
disp(k)

disp('the values of w1 and w2')
disp(x(1))
disp(x(2))

disp('final f(x)')
disp(fnc(x))

%----------------------------PLOTTING THE FUNCTION -------
x=linspace(-50,50,1000); 
y=linspace(-50,50,1000);
xlabel('x axis')
ylabel ('y axis')
zlabel('z ekseni')
fnc = @(x,y) ((x.^2 + y.^2 + 1.636.*x.*y)./2)-0.8182.*x-0.354.*y+0.1250; %function itself with an handle
gx= @(x,y)  x+(1.6365/2).*y -0.8182; %x gradient (x partial derivative)
gy= @(x,y) y+(1.6365/2).*x -0.354; %y gradient (y partial derivative) 
%grad=gx + gy;

figure(1)
[X,Y] = meshgrid(x,y);
f=fnc(X,Y);
surf(x,y,f);
shading flat;


figure(2)
plot(xdump,ydump,'LineWidth',2,'color','0,0,0')
legend('x[k] vs y[k]')
xlabel('x[k]')
ylabel('y[k]')

figure(3)
hold on
plot(zdump,'LineWidth',2)
yline(0,'r','LineWidth',3)
plot(gdump,'-.k','LineWidth',2)
legend('z[k]','zero','\nabla grad[k]')
axis([1 max([length(gdump),length(zdump)]) min([min(zdump),min(gdump)])-0.1 max([max(zdump),max(gdump)])+0.1])

figure(4)
plot(xdump,'k','LineWidth',1.2)
hold on
plot(ydump,'-.','LineWidth',1.2,'color','0,0,0')
legend('x[k]','y[k]')
xlabel('k')
axis([1 max([length(xdump),length(ydump)]) min([min(xdump),min(ydump)])-0.1 max([max(xdump),max(ydump)])+0.1])
%-----------------------------------

