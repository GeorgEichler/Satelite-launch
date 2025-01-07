g0 = 9.81;

%Setting common specific impulse and structural ratio = mEmpty/(mEmpty+mFuel)
Isp = 300;
epsilon = 0.15;

%r is payload ratio 
r = 0.05;
numStages = 10;
stages = 1:numStages;
Deltav = zeros(1,numStages);

for n = 1:numStages
    %Calculating Deltav, equation taken from Orbital mech. for Engineering, p. 676
    Deltav(n) = Isp*g0*n*log(1 / ( r^(1/n)*(1-epsilon) +epsilon ));
end

%theoretical limiting velocity, taken from book from above, p. 677
vMax = Isp*g0*(1-epsilon)*log(1/r);
disp(vMax);

%plotting Deltav in relation to maximal possible velocity
p1 = scatter(stages,Deltav/vMax, 'DisplayName', ...
    ['$\epsilon = $',num2str(epsilon),',r = ',num2str(r)]);
hold on;
plot(stages,Deltav/vMax);
xlabel('Number of stages','DisplayName','Smooth version');
ylabel('Final velocity for stage similar rocket');
legend(p1,'Location','northwest','Interpreter','latex');

figure;

vfinal = 10*10^3; %final velocity in m/s
r = zeros(1,numStages);


for n = 1:numStages
    %Rearranging equation for Deltav from above to r
    r (n) = ( (exp(-vfinal/(g0*Isp*n)) -epsilon ) / (1 - epsilon))^n;
end

%Get only the reasonable indices of r
positive_indices = r >= 0;

stages_restricted = stages(positive_indices);
r_restricted = r(positive_indices);

%plot payload ratio for each stage
p1 = scatter(stages_restricted,r_restricted,'DisplayName', ...
    ['$I_{sp}$ = ', num2str(Isp),'s', ', $\epsilon$ = ', num2str(epsilon), ...
    ', $v_{final}$ = ', num2str(vfinal / 1000), ' km/s']);
hold on;
plot(stages_restricted,r_restricted);
xlabel('Number of stages');
ylabel('Payload/total mass ratio');
xlim([1 numStages]);
legend(p1,'Location','southeast','Interpreter','latex');