clear
close all

%create the prior and parameterize it to build the sampling 

D = 0.05;   % noise corresponding to corrosion of parametric WM rep, diffusion
T = 5;      % delay time (seconds)
dt = 0.01;  % timestep 
nt = round(T/dt)+1; % number of timesteps per sim

% now parameterize the prior
m = 4;  % number of peaks on the periodic prior
Am = 1; % amplitude of the prior
n=4;
A=1;
offset=0; %for testing landscapes that are already heterogeneous but offset
pr = @(x) exp(Am*cos(m*x)); %prior function
ut=@(x)-A*cos(n*x-offset); %static landscape
nf = integral(pr,-pi,pi); %integral of the prior (for normalizing)
cpr = @(x)integral(pr,-pi,x)/nf; %cumulative probability distribution function


% build sampling
Ntrial = 100;   % number of trials 

% now build a discretized version of the prior and check sampling from it
sampres = 5e3;  % discretization of the prior function
xsamp = linspace(-pi,pi,sampres); %angle values
Fsamp = zeros(1,sampres);
for j=1:sampres, Fsamp(j) = cpr(xsamp(j)); end %builds the cumulative probability 

% now we make an object that we can sample from 
pdist = makedist('PiecewiseLinear','x',xsamp,'Fx',Fsamp);

% use the random function to sample from the above distribution and plot
randinp = random(pdist,1,Ntrial); % this identifies the angle theta for each simulation

figure; hold on, plot(xsamp,pr(xsamp),'k','LineWidth',2);  yticks([])
plot(randinp,ones(1,Ntrial),'.')
legend('Prior','Samples')
xlim([-pi pi])% and compare to true distribution
set(gca,'fontsize',30);set(gca,'ticklabelinterpreter','latex')
xlabel('$\theta$','fontsize',30,'interpreter','latex');
ylabel('$P(\theta)$','fontsize',30,'interpreter','latex');

figure; histogram(randinp,-pi:pi/16:pi)
xlim([-pi pi])% and compare to true distribution
set(gca,'fontsize',30);set(gca,'ticklabelinterpreter','latex')
xlabel('$\theta$','fontsize',30,'interpreter','latex');
ylabel('Sample Counts','fontsize',30,'interpreter','latex');

%describes the wells in the updating landscape
kap =8;  % amplitude of von mises to be used
sh = 1/4;   % shift
sc = 5;     % scaling
xs = linspace(-pi,pi,1000); dx = xs(2)-xs(1);

%initialize parameters
pot=zeros(Ntrial+1,length(xs));
pot(1,:) = 1/2/pi+0*xs; %initial potential landscape (flat)
 figure, hold on,  xlim([-pi pi]); yticks([])
set(gca,'fontsize',24); set(gca,'ticklabelinterpreter','Latex')
xlabel('$\theta$','fontsize',30,'interpreter','latex');
ylabel('$U(\theta)$','fontsize',30,'interpreter','latex');
plot(xs,pot(1,:),'linewidth',4,'Color','k');


%to track the responses and compute error
r_out = zeros(Ntrial,1); %responses at end of trial
r_dist = r_out; %distortion
for l=1:Ntrial
    inp = randinp(l);
    x=inp; 

    pgrad = gradient(pot(l,:),dx); %find the derivative of the current potential landscape
    
    % simulate the delayed response trial
    for j=1:nt-1 
        x=mod(x-dt*interp1(xs,pgrad,x)+sqrt(dt*2*D)*randn+pi,2*pi)-pi; 
    end
    r_out(l) = x;%response at end with drift and diffusion
    
    % add the distortion value for this trial
    r_dist(l) = (1-cos(x-inp)); %distortion for this trial


    %update the potential landscape using the previous trial input
    pot(l+1,:) = pot(l,:)+sc*(sh-exp(kap*(cos(xs-inp)-1)))/l;
    pot(l+1,:) = pot(l+1,:)/dx/sum(pot(l+1,:));
    
end

%plot heatmap
figure;  pcolor(xs,1:Ntrial+1,pot);  shading flat; colorbar
 xlabel('$\theta$','fontsize',30,'interpreter','latex');
ylabel('Trial','fontsize',30,'interpreter','latex'); xticks([-2,0,2]); xticklabels({'-100','0','100'})
set(gca,'fontsize',30); set(gca,'ticklabelinterpreter','Latex'); set(gca,'Ydir','reverse')


%plot final landscape
 figure, hold on,  xlim([-pi pi]); yticks([])
set(gca,'fontsize',24); set(gca,'ticklabelinterpreter','Latex')
xlabel('$\theta$','fontsize',30,'interpreter','latex');
ylabel('$U(\theta)$','fontsize',30,'interpreter','latex');
plot(xs,pot(Ntrial+1,:),'linewidth',4,'Color','k');
set(gca,'FontSize',30); set(gca,'TickLabelInterpreter','Latex')
