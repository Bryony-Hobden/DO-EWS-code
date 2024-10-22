%% The only script you need to run/ make changes to
clear all
%%
%set seed
rng(4)

%%
%y(1) = tg, y(2) = tnaw, y(3) = phi, y(4)= Ia, y(5) = lambda, y(6) = C
%initial conditions
y0 = [-26.09; 1.13; 2.9; -32.86; 0; 0];

%load pdfs for speed sampling
load("pdfs.mat", 'scaledbenthic', 'filteredbenthic','pdfs', 'benthic_sp_2')

p.scaledbenthic=scaledbenthic;
p.filteredbenthic=filteredbenthic;
p.pdfs=pdfs;
p.benthic_sp_2 =benthic_sp_2;

% DO_solve_IE(y0, Ig_0, dt, epsilon, sigma, years(KY))
p.Ig_0=-34;
p.dt=0.0001;
p.epsilon=0.01;
p.sigma=0;
p.years=60;
p.beta= 1;
p.tau=10;

% p.TNAW_const = 0 TNAW varies across last glacial period, p.TNAW_const = 1
% TNAW is fixed across last glacial period at a value specified by p.TNAW =

p.TNAW_const = 1;
p.TNAW = 20;

%% DO simulation 
[t,y, starfunc, Ig, Cdot] = DO_solve_IE(y0, p);
%%
save("simulation.mat");
%% Comparison plot of NGRIP and simulation

NGRIP20y = readtable('NGRIP.csv','NumHeaderLines',5, 'VariableNamingRule','preserve');  % skips the first 5 rows of data
NGRIP20y = table2array(NGRIP20y);
is_start_times= [11700, 14680, 23340, 27780, 28900, 32500, 33740, 35480, 38220, 40160, 41460, 43340, 46860, 54220, 55800, 58040, 59080, 64100, 69620, 72340, 76440, 84760, 104040, 108280, 115380];
is_end_times = [12900, 23100, 27540, 28600, 32040, 33360, 34740, 36580, 39900, 40800, 42240, 44280, 48340, 55400, 56500, 58560, 63840, 69400, 70380, 74100, 77760, 87600, 105440, 110640, 119140];


is_durations = circshift(is_start_times,-1) - is_end_times;
is_durations = is_durations(1:end-1);

%% Identify where 'DO events' are in simulation
diff=0;

DO_events = [];

for i = 1:(length(Ig)-50)
    diff = Ig(i+50)-Ig(i);
    if diff>1
        DO_events = [DO_events 120-(i/10000)];
    end
end

DO_events = round(DO_events,1);
DO_events = unique(DO_events);

%% Plot of NGRIP vs simulation isotope ratio
figure
x0=10;
y0=10;
width=1000;
height=400;
set(gcf,'position',[x0,y0,width,height])
set(gcf,'color','#E7ECEF');
subplot(2, 1, 1)
hold on
set( gca, 'xdir', 'reverse' )
set(gca,'FontSize',10, 'FontName', 'Outfit')
plot(NGRIP20y(:,1), NGRIP20y(:,2))
xline(is_start_times/1000, Color=[.7 .7 .7], LineWidth=1)
xlim([0 120])
ylim([-49, -33])
ylabel('NGRIP \delta^{18}O')
hold off


subplot(2, 1, 2)
hold on
set( gca, 'xdir', 'reverse' )
set(gca,'FontSize',10, 'FontName', 'Outfit')
plot(-t, Ig,"LineWidth",1)
xline(DO_events, Color=[.7 .7 .7], LineWidth=1)
xlim([0 120])
xlabel('Time [kyr B2k]')
ylabel('Simulated I_G')
ylim([-49, -33])
hold off

%% Creating stable/untsbale manfifolds for phase portrait
stable1= linspace(4.5,2, 1000);
stable2= linspace(-3,-1, 1000);
stable3 = linspace(-2, 1, 1000);
unstable1 = linspace(1, 3, 1000);
unstable2 = linspace(-1,2, 1000);

for i =1:length(stable1)
    
    stable1c(i) = -(8/135)*(stable1(i))^3 +(4/45)*(stable1(i))^2 +(16/45)*stable1(i)+ (11/27);
    stable2c(i) = -(8/135)*(stable2(i))^3 +(4/45)*(stable2(i))^2 +(16/45)*stable2(i)+ (11/27);

    stable3c(i) = -2*(stable3(i)-1/2)*(stable3(i)-3/2);
    unstable1c(i) = -2*(unstable1(i)-1/2)*(unstable1(i)-3/2);
    unstable2c(i) = -(8/135)*(unstable2(i))^3 +(4/45)*(unstable2(i))^2 +(16/45)*unstable2(i)+ (11/27);
end
%% Phase portrait plots
a = (2-log(2))/2 + 0.01*(20- 10);
cnull = @(x) exp(2*(x-a));
xspan = linspace(-5,5,1000);      
figure
set(gcf,'color','#E7ECEF');
set(gca,'FontSize',10, 'FontName', 'Outfit')
hold on
plot(y(6,:), y(5,:), Color=[0.8500 0.3250 0.0980], LineWidth=1)
plot(stable1c, stable1, LineStyle='-', Color='black', LineWidth=1)
plot(stable2c, stable2, LineStyle='-', Color='black',LineWidth=1)
plot(stable3c, stable3, LineStyle='-', Color='black',LineWidth=1)
plot(unstable1c, unstable1, LineStyle='--',Color='black',LineWidth=1)
plot(unstable2c, unstable2, LineStyle='--',Color='black',LineWidth=1)
plot(y2(6,1:200000), y2(5,1:200000), Color = 'black', LineWidth=1.5)
yline(0)
xline(0)
plot(xspan, cnull(xspan), Color=[0 0.4470 0.7410], LineWidth=1)
%scatter(1,2,'filled', 'd', 'MarkerEdgeColor','none', 'MarkerFaceColor',[0 0 0])
xlabel('C')
ylabel('\lambda')
xlim([-0.1 1.2])
scatter(0.905, 1.357, 'filled','o','MarkerFaceColor','black')
%scatter(0.985986, 1.75966, 'filled','o','MarkerFaceColor','black')
%scatter(0.971, 2.318, 'filled','o',Color=[0 0.4470 0.7410])
grid on
hold off