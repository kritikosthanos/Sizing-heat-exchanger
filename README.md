# Sizing-heat-exchanger

clear all 

Fcph=2.17;%kJ/s K 

Q=[54.06 447 661.5 71.35 1304 ;%kJ/s[=] 

0.15 0.29 0.24 0.21 0.29;%kW/m2 K [=]U 

298.15 298.15 773.15 298.15 298.15;%K [=]Tina 

724.05 701.95 848.15 796.95 338.15];%K[=] Toutb 

 

 

for j=1:120 

w=perms(1:5); 

% j=52; 

n=1; 

T(n,j)=1894.15; 

for i=1:5 

    

%     cp=heatcapacity(T(n)) 

    n=n+1; 

    Tin=Q(3,w(j,i)); 

    Tout=Q(4,w(j,i)); 

    U=Q(2,w(j,i)); 

     

     

    T(n,j)=T(n-1,j)-(Q(1,w(j,i))/(Fcph)); 

    DT(i)=((T(n,j)-Tout)-(T(n-1,j)-Tin))/log((T(n,j)-Tout)/(T(n-1,j)-Tin)); 

    A(n-1,j)=Q(1,w(j,i))/(U*DT(i)); 

 

end 

 

Aol(j)=sum(A(:,j)); 

Amin=min(Aol) 

 

end 

 

 

 

 

10.5 Διαστασιολόγηση Εναλλάκτη Θερμότητας Νερού 

 

Κύριος κώδικας: 

clc 

clear all 

%sustaseis  

%     1        2      3      4    5      6      7 

%     CH4      H2O    CO     H2   CO2    O2     N2 

y=[0 1 0 0 0 0 0]; 

yout=[0 0.2264 0 0 0.045281 0.079724 0.648591 ]; 

N=0.466;% 

Nout=5.33e-2;%kgmole/s 

dt=[0.407 0.532 0.657 0.732 0.982 1.232 1.732]*0.0254;%inches->m 

  

%temperature-pressure 

T=25+273.15;%K 

Tout=1325.5;%K 

P=100000; 

Pout=100000; 

%lenght 

zi=0;zf=10; 

z=linspace(zi,zf,100000); 

hout=enthalpy(yout,Tout);%Kj/kmol 

h=enthalpy(y,T); 

  

x0=zeros(5,1); 

x0(1)=0;%A_m2 

x0(4)=N*h;%Kj/s 

x0(5)=Nout*hout;%Kj/s 

x0(2)=T;%K 

x0(3)=Tout;%K 

x0(6)=P; 

x0(7)=Pout; 

  

%mass matrix 

MASS=eye(7); MASS(2,2)=0; MASS(3,3)=0; 

Teks(length(z),length(dt))=1; 

options =odeset('reltol',1e-4,'abstol',1e-10,'mass',MASS); 

for i=1:7 

[z,x]=ode15s(@heat_exchanger,z,x0,options,dt(i)); 

Teks(:,i)=x(:,2); 

AA(:,i)=x(:,1); 

end 

i=1 

 

Για όλα τα συστατικά: 

clc 

clear all 

%sustaseis  

%     1        2      3      4    5      6      7 

%     CH4      H2O    CO     H2   CO2    O2     N2 

y=[0 1 0 0 0 0 0]; 

yout=[0 0.2264 0 0 0.045281 0.079724 0.648591 ]; 

N=0.466;% 

Nout=5.33e-2;%kgmole/s 

dt=[0.407 0.532 0.657 0.732 0.982 1.232 1.732]*0.0254;%inches->m 

  

%temperature-pressure 

T=25+273.15;%K 

Tout=1325.5;%K 

P=100000; 

Pout=100000; 

%lenght 

zi=0;zf=10; 

z=linspace(zi,zf,100000); 

hout=enthalpy(yout,Tout);%Kj/kmol 

h=enthalpy(y,T); 

  

x0=zeros(5,1); 

x0(1)=0;%A_m2 

x0(4)=N*h;%Kj/s 

x0(5)=Nout*hout;%Kj/s 

x0(2)=T;%K 

x0(3)=Tout;%K 

x0(6)=P; 

x0(7)=Pout; 

  

%mass matrix 

MASS=eye(7); MASS(2,2)=0; MASS(3,3)=0; 

Teks(length(z),length(dt))=1; 

options =odeset('reltol',1e-4,'abstol',1e-10,'mass',MASS); 

for i=1:7 

[z,x]=ode15s(@heat_exchanger,z,x0,options,dt(i)); 

Teks(:,i)=x(:,2); 

AA(:,i)=x(:,1); 

end 

i=1 

 

Υπολογισμός της ενθαλπίας: 

function h=enthalpy(y,T) 

clc 

% A user defined function to calculate 

% the specific enthalpy of the gases of 

% the SR/SOFC/UB systems 

% 

% Input arguments : y = vector of mole fractions of the gas mixture [-] 

%                   T = temperature [K] 

% Output argument : h = specific enthalpy of gas mixture [kJ/kmol] 

  

% Components 

%     1        2      3      4    5      6      7 

%     CH4      H2O    CO     H2   CO2    O2     N2 

  

h=0; 

href=[-74870   -241830 -110530 0 -393520 0 0]; % enthalpy of formation of the ideal das at 298.15 K in [kJ/kmol] 

Tref=298.15;                                   % temperature at which enthapy equals enthalpy of formation of ig 

% cp = C + B * T + A * T^2 

% constants A,B,C are fitted values to raw data from NIST Chemistry webbook 

% 

A=[0           0         2.5359e-6     1.5807e-6   -2.2051e-05  -7.1534e-06    0      ]; 

B=[0.05346     0.01133   0.00300006    -0.00050303  0.05228      0.016494    0.0053629]; 

C=[19.676      29.704    27.774        29.058       23.853       24.991         27.237]; 

for j=1:7 

    h= h + y(j) * ((A(j)/3) * (T^3-Tref^3) + (B(j)/2) * (T^2-Tref^2) + C(j) * (T-Tref) + href(j)); 

end 

end 

 
