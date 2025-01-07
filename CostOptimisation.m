g0 = 9.81;

%Setting common specific impulse and structural ratio = mEmpty/(mEmpty+mFuel)
Isp = 300;
epsilon = 0.1;
numStages = 10;
stages = 1:numStages;

%Using LEO orbital velocity and add recommended velocity loss
vfinal = (7.8 + 1.5)*10^3; %final velocity in m/s
r = zeros(1,numStages);

for n = 1:numStages
    %Calculate payload ratio for each stage by the rearranged equation 
    %from Orbital Mechanics for Engineering Students, p. 676
    r (n) = ( (exp(-vfinal/(g0*Isp*n)) -epsilon ) / (1 - epsilon))^n;
end

%Get only the reasonable indices of r
positive_indices = r >= 0;

stages_restricted = stages(positive_indices);
r_restricted = r(positive_indices);

p1 = scatter(stages_restricted,r_restricted,'DisplayName', ...
    ['$I_{sp}$ = ', num2str(Isp),'s', ', $\epsilon$ = ', num2str(epsilon), ...
    ', $v_{final}$ = ', num2str(vfinal / 1000), ' km/s']);
hold on;
plot(stages_restricted,r_restricted);
xlabel('Number of stages');
ylabel('Payload/total mass ratio');
xlim([1 numStages]);
legend(p1,'Location','southeast','Interpreter','latex');

%Set a payload, 140000 is the about the one of the Saturn V
mPayload = 140000;
%payload for each stage, last entry is payload of whole rocket
%and saving the payload for 1, 2,..., numStages rocket
%E.g. for a 3 stage rocket the payload for stage 2 is given by 
%the actual mPayload + mass of stage 3
PayloadMasses = zeros(numStages,numStages+1);
%Mass of each stage and for all n stage rockets
m = zeros(numStages,numStages);
mEmpty = zeros(numStages,numStages);
mPropellant = zeros(numStages,numStages);

for k = 1:numStages
    
    m0 = mPayload/r(k);
    PayloadMasses(k,1) = m0;
    PayloadMasses(k,2) = (mPayload * m0^(numStages-1))^(1/numStages);
    
    if k ~= 1
        for n = 2:k
            PayloadMasses(k,n+1) = (PayloadMasses(k,2))^n/(m0)^(n-1);
        end
    end
    
    for n = 1:k
        m(k,n) = PayloadMasses(k,n) - PayloadMasses(k,n+1);
        mEmpty(k,n) = epsilon*m(k,n);
        mPropellant(k,n) = (1-epsilon)*m(k,n);
    end
    
    %total mass of each rocket
    totalMass = sum(m,2);
end

%cost per kg of propellant
pricePropellant = 15.99;
%sum is taken for each row
costPropellant = pricePropellant*sum(mPropellant,2);
%turn it into a row vector
costPropellant = costPropellant';

%cost per kg of empty stage
priceStage = 40;
%sum is taken for each row
costStage = priceStage*sum(mEmpty,2);
%turn it into a row vector
costStage = costStage';

%Another approach using fixed costs assuming 5 engines
%Uses data from the engine word document on the oneDrive folder
%this seems to be the more realistic option
priceEngine = 15*10^6;
costStage = (1:numStages)*5*priceEngine;

costRatio = (costPropellant + costStage) ./ r;

costRatio_restricted = costRatio(positive_indices);

figure;
p1 = scatter(stages_restricted,costRatio_restricted,'DisplayName', ...
    ['$I_{sp}$ = ', num2str(Isp),'s', ', $\epsilon$ = ', num2str(epsilon), ...
    ', $v_{final}$ = ', num2str(vfinal / 1000), ' km/s']);
hold on;
plot(stages_restricted,costRatio_restricted);
xlabel('Number of stages');
ylabel('Cost ratio');
xlim([1 numStages]);
legend(p1,'Location','northwest','Interpreter','latex');