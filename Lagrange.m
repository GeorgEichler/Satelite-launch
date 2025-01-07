function [m0, m, mEmpty, mFuel, mEnd] = Lagrange(vEnd, Isp, epsilon)
%Input variables:
%vEnd - Delta V
%Isp - specific Impulse for each stage (1 x s) vector
%epsilon - structural ratios for each stage (1 x s) vector
%
%All output masses are with an unity payload mass (simply scalable for
%arbitrary payload mass by multyplying it with payload mass)
%Output variables:
%m0 - total mass of rocket
%m - mass for each stage (1 x s) vector
%mEmpty - empty mass for each stage (1 x s) vector
%mFuel - fuel mass for each stage (1 x s) vector



g0 = 0.00981; %gravitational acceleration on sea-level
c = g0*Isp; %exhaust velocity in km/s
%Set payload to unity and safe it in m0
mPayload = 50000;
mEnd = mPayload;
m0 = mPayload;

%numerically solve the constraint equation
mu1 = fzero(@(x) constraint_equation(x, vEnd, Isp, epsilon), [1e-3,1e3]);

%Calculate the mass ratio for each stage 
%Defined by n(i) = m0_i/(m_0i - mFuel_i)
%where m0_i is total mass of rocket after discarding stages 1,2,...,i-1
%and mFuel_i is the mass of the fuel of stage i
n = (mu1*c - 1) ./ (mu1*epsilon.*c);

%Get the number of stage
s = length(Isp);
%preinitialise the mass for each stage, i.e.
%m(i) = mass of stage i 
m = zeros(1,s);

for i = s:-1:1
    m(i) = (1 - n(i)) / (epsilon(i)*n(i) - 1)*mPayload;
    
    mPayload = mPayload + m(i);
end

mEmpty = epsilon.*m;
mFuel = m - mEmpty;

m0 = m0 + sum(m);

return;



%Equation ensuring final velocity equal our target velocity
function f = constraint_equation(x, v, Isp, epsilon)

g0 = 9.81;
c = g0*Isp/1000;

if x <= 0 || any((x*c - 1) <= 0)
    f = -1;
else
    f = sum( c.* log( (x*c - 1) ./(x*epsilon.*c) ) ) - v;
end

return;