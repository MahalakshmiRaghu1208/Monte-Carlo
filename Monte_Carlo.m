close all;
clear all;
%variables 
t=0;
dt=0.01;

m=1;%mass

g=9.80665;
%no of samples
n=1000;
V=[];
Vx=[];
Vy=[];
X=[];
Y=[];
T=[];
Tt=[];

mean=2.5;
sd=0.3;
mean1=45;
sd1=5;
Th=[];
for i=1:n
 % velocity and angle components
v=normrnd(mean,sd);
theta=normrnd(mean1,sd1);
%velocity is a gaussian distribution
vx=v*cosd(theta);
vy=v*sind(theta);
Th=[Th theta];
%position components
x=0;
y=0;
%finding tmax and tf
tmax=vy/g;
tf=2*tmax;
T=t:dt:tf;% creating timestep matrix
for k=1:length(T)
    V(k,i)= v;
    Vx(k,i)=vx;
    Vy(k,i)=vy;
    X(k,i)=x;
    Y(k,i)=y;
    vy=vy-g*dt;
    % simple 2d propagation with euler
    x=x+vx*dt;
    y=y+vy*dt;
    v=((vx.^2)+(vy.^2)).^0.5;
end
Tt=[Tt length(T)];
end

%Downrange achieved for each sample
D=[];
figure(1)
for l=1:n
plot(X(1:Tt(l),l),Y(1:Tt(l),l)); 
D=[D X(Tt(l),l)];
hold on;
end
D=D';

%Assume success criterion is mean and deviation
smean=0.6;
ssd=0.10;
%defining the success matrix
S=[];
S1=[];
F1=[];
S2=[];
F2=[];

for k=1:n
    if D(k,1)<=(smean+ssd) && D(k,1)>=(smean-ssd) 
    S=[S D(k,1)];
    end
end
S=S';
S=sort(S);

%the velocity contributing to success/failure
for k=1:n
    if D(k,1)<=(smean+ssd) && D(k,1)>=(smean-ssd) 
    S1=[S1 V(1,k)];
    else
    F1=[F1 V(1,k)];
    end
end
S1=S1';
F1=F1';

%the theta contributing to success/failure
for k=1:n
    if D(k,1)<=(smean+ssd) && D(k,1)>=(smean-ssd) 
    S2=[S2 Th(1,k)];
    else
    F2=[F2 Th(1,k)];
    end
end
S2=S2';
F2=F2';
points=1:1:n;
figure(2)
for k=1:n
    if D(k,1)<=(smean+ssd) && D(k,1)>=(smean-ssd) 
    plot(V(1,k),Th(k),'b*')
    hold on;
    else
    plot(V(1,k),Th(k),'r*')
    end
end
xlabel("Velocity")
ylabel("Theta")

figure(7)
for k=1:n
    if D(k,1)<=(smean+ssd) && D(k,1)>=(smean-ssd) 
    plot(V(1,k),points(k),'b*')
    hold on;
    else
    plot(V(1,k),points(k),'r*')
    end
end
xlabel("Velocity")
ylabel("Simulation run#")

figure(8)
for k=1:n
    if D(k,1)<=(smean+ssd) && D(k,1)>=(smean-ssd) 
    plot(Th(k),points(k),'b*')
    hold on;
    else
    plot(Th(k),points(k),'r*')
    end
end
xlabel("Theta")
ylabel("Simulation run#")



%Probability fit for the Downrange
pd=fitdist(D(:,1),'Normal');
mn=pd.mu;
sig=pd.sigma;
disp("the mean of the output data:"+ mn);
disp("the standard deviation of the output data:"+ sig);
p=pdf('Normal',D(:,1),mn,sig);
r=max(p);

figure(3);
%plotting the function 
plot(D(:,1),p,'b.',[mn mn],[0 r],'r-');%basic plot
hold on;
M=[];
for k=1:length(S)
    ap=interp1(D(:,1),p,S(k),'nearest');
    M=[M ap];
end
M=M';
N=[];
O=[];
for l=1:length(S)
    N=[N S(l)];
    O=[O M(l)];
end
A=[N(1) N N(length(N))];
B=[0 O 0];
fill(A,B,'r','FaceAlpha',0.3);
xlabel('Downrange achieved(m)');
ylabel('Probability density function');



% pd1=fitdist(S1(:,1),'Normal');
% mn1=pd1.mu;
% sig1=pd1.sigma;
% h1=1.06*sig1*(n^(-1/5));
% 
% pd2=fitdist(F1(:,1),'Normal');
% mn2=pd2.mu;
% sig2=pd2.sigma;
% h2=1.06*sig2*(n^(-1/5));

d=0.01;
gp=mean-4*sd:d:mean+4*sd;
[f1,xi1]=ksdensity(S1,gp);
[f2,xi2]=ksdensity(F1,gp);

gp1=mean1-4*sd1:d:mean1+4*sd1;
[f3,xi3]=ksdensity(S2,gp1);
[f4,xi4]=ksdensity(F2,gp1);


figure(4)
plot(xi1,f1,'LineWidth',2)
hold on;
plot(xi2,f2,'LineWidth',2)
xlabel('velocity variation')
ylabel('PDF')

figure(5)
plot(xi3,f3,'LineWidth',2)
hold on;
plot(xi4,f4,'LineWidth',2)
xlabel('Theta variation')
ylabel('PDF')

% f1=f1/length(S1);
% f2=f2/length(F1);
% 
% f3=f3/length(S2);
% f4=f4/length(F2);


fun1=[];
for i=1:length(gp)
    fun1=[fun1 abs(f1(i)-f2(i))];
end 

fun2=[];
for i=1:length(gp1)
    fun2=[fun2 abs(f3(i)-f4(i))];
end 


J1=sum(fun1,"all")*d;
disp("the cost function of the output data due to velocity variation:"+ J1);

J2=sum(fun2,"all")*d;
disp("the cost function of the output data due to theta variation:"+ J2);

%finding the success confidence percentage 
l1=smean+ssd;
l2=smean-ssd;

Z1=(l1-mn)/sig;
Z2=(l2-mn)/sig;

Confidence=normcdf(Z1)-normcdf(Z2);
disp("the confidence level of the success"+ Confidence);
