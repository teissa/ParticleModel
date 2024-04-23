% This code from Eisssa and Kilpatrick 2023 builds a heterogeneous particle
% landscape and runs simulations of the sotchastic differential equation of
% a particle within the describe potential landscape. 


%%%%%% building the landscape
n=4; %number of wells
A=1; %amplitude
ut=@(x)-A/n*cos(n*x);
% now build a discretized version of the prior and check sampling from it
sampres = 1e3;  % discretization of the prior function
xsamp = linspace(-pi,pi,sampres); %angle values
figure;
hold on, plot(xsamp,ut(xsamp),'k','LineWidth',2);  yticks([])
xlim([-pi pi])% and compare to true distribution
set(gca,'fontsize',30);set(gca,'ticklabelinterpreter','latex')
xlabel('Target($\theta$)','fontsize',30,'interpreter','latex');
ylabel('$U_4(\theta)$','fontsize',30,'interpreter','latex');

%%%%%% simulate the trajectories of particles within the landscape
D = 0.05;   % noise 
T = 5;      % delay time (seconds)
Nsim = 20;   % number of sims
dt = 0.01;  % timestep of stochastic sims
nt = round(T/dt)+1; % number of timesteps per sim
tim=0:dt:T;

% use the random function to sample from the above distribution and plot
% sampling
randind=randi(256,1,Nsim);
randinp = xsamp(randind); % this identifies the angle theta for each simulation
xfs = linspace(-pi,pi,256); %angles discretized to find distortion
dx = xfs(2)-xfs(1); %angle step size

figure; hold on;
ylabel('Time(s)','Interpreter','Latex')
xlabel('Target','Interpreter','Latex')
set(gca,'FontSize',30); set(gca,'TickLabelInterpreter','Latex')
axis([-pi pi 0 T ])
for l=1:Nsim
        inp = randinp(l); %angle theta from distribution creted above
        x=zeros(1,nt); %starting point of particle
        x(1)=inp;
        for j=2:nt
            x(j)=mod(x(j-1)-dt*A*sin(n*x(j-1))+sqrt(dt*2*D)*randn+pi,2*pi)-pi; %movement of the particle, normalized by circular distance
        end
        
        plot(x,tim,'Color','k')
end
xticks([-2,0,2]); xticklabels({'-100','0','100'})

%note, the basic plotting function above will look distorted for particles
%that are near the edges (-pi and pi), but the values will be correct in
%the matrix. 