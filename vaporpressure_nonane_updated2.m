clc
clear
R=8.314;   %universal gas constant J/molK
Tc=594.6;  %critical temperature of C9H20 in Kelvin
Pc=22.8*1.01325e5; %critical pressure of nonane in Pascal
w=0.444;   %accentric factor of nonane

%Using the peng-robinson equation of state
f=0.37464+1.54226*w-0.26992*w^2;   %calculate f as a function of accentric factor
b=0.07780*R*Tc/Pc;                 %calculate the constant b in the PR equation of state
T=linspace(220,535,250);           %specify temperature range from which vapor pressure is to be found
Tr=linspace(0,0,250);              %initialize reduced temperature
for i=1:250
Tr(i)=T(i)/Tc;                          %calculate the reduced temperature
end

Psat=linspace(0,0,250);            %initialize vaporization pressure
Psat(1)=0.46;



for i=2:250
    P=Psat(i-1);           %guess a random value for the vaporization pressure
    fl=1.2;
    fg=1;
    while abs((fl/fg)-1)>0.0001
      Liquidprops=PengRobinson(T(i),P,Tc,Pc,w,1);
      Gasprops=PengRobinson(T(i),P,Tc,Pc,w,0);
      Zl=Liquidprops(1);
      phil=Liquidprops(2);
      Zg=Gasprops(1);
      phig=Gasprops(2);
      fl=phil*P;
      fg=phig*P;
      P=(sqrt(fl*fg)/Psat(i-1))*P;
   end
    Psat(i)=P;

end
Psatmpa=Psat/1e6;
% PengRobinson.m : calculates the compressibility factor and fugacity coefficient  
% of a pure compound with the Peng Robinson equation of state (PR EOS)
% Authors: Abraham Adu-Mills
%
% function result = PengRobinson(T,P,Tc,Pc,w,Liquido)
% Parameters: T,P,w,Tc,Pc,w,MW,Liquido
% T: Temperature [=] K                                          
% P: Presure [=] Pa                                             
% Tc: critical temperature [=] K                               
% Pc: critical presure [=] Pa                                   
% w: accentic factor
% Liquido:  if Liquido == 1, then calculates liquid fugacity;  
%           if Liquido == 0 then calculates vapor fugacity
% Example:
% [Z fhi] = PengRobinson(273,2*1.013*1e5,304.21,7.382*1e6,0.044,1)
% Plot solution
figure(1)
hold on
plot(T,Psatmpa, 'LineWidth',1,'Color', 'Black');
xlabel('Temperature(K)')
ylabel('Pressure(Mpa)')
title('Vaporization Curve of Nonane (C9H20) from 220K to 535K')
hold off

function output = PengRobinson(T,P,Tc,Pc,w,Liquido)
R = 8.314; % gas constant [=] J/(mol K)
% Reduced variables
Tr = T/Tc ;
ZR=0;
% Parameters of the EOS for a pure component
m = 0.37464 + 1.54226*w - 0.26992*w^2;
alfa = (1 + m*(1 - sqrt(Tr)))^2;
a = 0.45724*(R*Tc)^2/Pc*alfa;
b = 0.0778*R*Tc/Pc;
A = a*P/(R*T)^2;
B = b*P/(R*T);
% Compressibility factor
Z = roots([1 -(1-B) (A-3*B^2-2*B) -(A*B-B^2-B^3)]);

for i = 1:3
   if isreal(Z(i))
   	ZR = [Z(i)];  
   else
       disp('complex Z value detected')
   end
end
if Liquido == 1
    Z = min(ZR);   
else
    Z = max(ZR);
end
% Fugacity coefficient
fhi = exp(Z - 1 - log(Z-B) - A/(2*B*sqrt(2))*log((Z+(1+sqrt(2))*B)/(Z+(1-sqrt(2))*B)));
if isreal(fhi)
 
    output = [Z fhi];
else
    disp('No real solution for "fhi" is available in this phase')
    output =['N/A' 'N/A' 'N/A'];
end
end 

    
