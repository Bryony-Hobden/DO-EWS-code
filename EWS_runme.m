clear all
close all

%%
load simulation.mat

%%
NGRIP20y = readtable('NGRIP.csv','NumHeaderLines',5, 'VariableNamingRule','preserve');  % skips the first 5 rows of data
NGRIP20y = table2array(NGRIP20y);
is_start_times= [11700, 14680, 23340, 27780, 28900, 32500, 33740, 35480, 38220, 40160, 41460, 43340, 46860, 54220, 55800, 58040, 59080, 64100, 69620, 72340, 76440, 84760, 104040, 108280, 115380];
is_end_times = [12900, 23100, 27540, 28600, 32040, 33360, 34740, 36580, 39900, 40800, 42240, 44280, 48340, 55400, 56500, 58560, 63840, 69400, 70380, 74100, 77760, 87600, 105440, 110640, 119140];


is_durations = circshift(is_start_times,-1) - is_end_times;
is_durations = is_durations(1:end-1);
%% Code to identify simulated DO events 
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
%% Comparison of NGRIP vs simulated Ig time series
figure
set(gcf,'color','#E7ECEF');
subplot(2, 1, 1)
hold on
set( gca, 'xdir', 'reverse' )
set(gca,'FontSize',10, 'FontName', 'Outfit')
set(gca(), ...
    'Layer','top')
plot(NGRIP20y(:,1), NGRIP20y(:,2))
xlim([0 120])
xline(is_start_times/1000, Color=[.7 .7 .7], LineWidth=1)
ylim([-49, -33])
ylabel('NGRIP \delta^{18}O')
hold off

subplot(2, 1, 2)
hold on
set( gca, 'xdir', 'reverse' )
set(gca,'FontSize',10, 'FontName', 'Outfit')
set(gca(), ...
    'Layer','top')
rectangle('position', [71.03 -49 77.1-71.03 49-33], EdgeColor='red', FaceColor='#FFC0CB')
rectangle('position', [16.79  -49 20.07-16.69 49-33], EdgeColor='blue', FaceColor='#add8e6')
rectangle('position', [44.42 -49 47.27-44.42 49-33], EdgeColor='red',FaceColor='#FFC0CB')
rectangle('position', [23.83 -49 29.15-23.83 49-33], EdgeColor='blue',FaceColor='#add8e6')
plot(120:-0.0001:0, Ig,"LineWidth",1)
xlim([0 120])
xlabel('Time [kyr B2k]')
ylabel('Simulated I_G')
xline(DO_events(1:19), Color=[.7 .7 .7], LineWidth=1)
xline(DO_events(21:end), Color=[.7 .7 .7], LineWidth=1)
ylim([-49, -33])
hold off
%% Pre function stuff
compiled_t_ig = [120:-0.0001:0; Ig]; 
a_int1 = find(round(compiled_t_ig(1,:),4) ==77.1297);
b_int1 = find(round(compiled_t_ig(1,:),4) ==71.0304);
INT1 = compiled_t_ig(:, a_int1:b_int1);

%% Function stuff
[x_5, W_5] = bry_func(INT1, 0.005, 80);
%%
plot(t(a_int1:b_int1), y(6, a_int1:b_int1))
%% interstadial 1
figure
subplot(3,1,1)
set( gca, 'xdir', 'reverse' )
set(gca,'FontSize',10, 'FontName', 'Outfit')
set(gca(), ...
    'Layer','top')
hold on
patch([x_5(2,1) x_5(2,80) x_5(2,80) x_5(2,1)], [-42 -42 -37 -37],  [.8 .8 .8], edgecolor = 'none')
plot(INT(1,:), INT(2,:))
%plot(x_5(2,:), x_5(8,:))
xlim([71.0304 77.1297])
ylim([-41 -37.5])
ylabel('Simulated I_g')


subplot(3,1,2)
set( gca, 'xdir', 'reverse' )
set(gca,'FontSize',10, 'FontName', 'Outfit')
set(gca(), ...
    'Layer','top')
hold on
patch([x_5(2,1) x_5(2,80) x_5(2,80) x_5(2,1)], [0 0 0.5 0.5],  [.8 .8 .8], edgecolor = 'none')
plot(INT(1,:), y(6, a_int:b_int))
xlim([71.0304 77.1297])
ylim([0 0.5])
ylabel('C')


subplot(3,1,3)
set( gca, 'xdir', 'reverse' )
set(gca,'FontSize',10, 'FontName', 'Outfit')
set(gca(), ...
    'Layer','top')
hold on
patch([x_5(2,1) x_5(2,80) x_5(2,80) x_5(2,1)], [0 0 0.21 0.21],  [.8 .8 .8], edgecolor = 'none')
patch([x_5(2,80:end) fliplr(x_5(2,80:end))], [(x_5(6,80:end))/1000 fliplr(x_5(7,80:end))/1000], [0.97 0.71 0.64], EdgeColor = 'none')
plot(x_5(2,:), x_5(5,:)/1000, color = 'r', LineWidth= 1.2)
plot(x_5(2,:), x_5(6,:)/1000, Color='r', LineWidth= 0.8)
plot(x_5(2,:), x_5(7,:)/1000, Color='r', LineWidth= 0.8)
xlim([71.0304 77.1297])
ylim([0 0.18])
ylabel('\alpha (1/yr) ')
xlabel('Time (kyr b2k)')
%% interstadial 2
a_int2 = find(round(compiled_t_ig(1,:),4) == 47.2651);
b_int2 = find(round(compiled_t_ig(1,:),4) == 44.4166);
INT2 = compiled_t_ig(:, a_int2:b_int2);

[x_5_int2, W_5_int2] = bry_func(INT2, 0.005, 80);

subplot(3,1,1)
set( gca, 'xdir', 'reverse' )
set(gca,'FontSize',10, 'FontName', 'Outfit')
set(gca(), ...
    'Layer','top')
hold on
patch([x_5_int2(2,1) x_5_int2(2,80) x_5_int2(2,80) x_5_int2(2,1)], [-42 -42 -37 -37],  [.8 .8 .8], edgecolor = 'none')
plot(INT2(1,:), INT2(2,:))
%plot(x_5(2,:), x_5(8,:))
xlim([44.4166 47.2651])
ylim([-42 -37.5])
ylabel('Simulated I_g')


subplot(3,1,2)
set( gca, 'xdir', 'reverse' )
set(gca,'FontSize',10, 'FontName', 'Outfit')
set(gca(), ...
    'Layer','top')
hold on
patch([x_5_int2(2,1) x_5_int2(2,80) x_5_int2(2,80) x_5_int2(2,1)], [0 0 0.5 0.5],  [.8 .8 .8], edgecolor = 'none')
plot(INT2(1,:), y(6, a_int2:b_int2))
xlim([44.4166 47.2651])
ylim([0 0.5])
ylabel('C')


subplot(3,1,3)
set( gca, 'xdir', 'reverse' )
set(gca,'FontSize',10, 'FontName', 'Outfit')
set(gca(), ...
    'Layer','top')
hold on
patch([x_5_int2(2,1) x_5_int2(2,80) x_5_int2(2,80) x_5(2,1)], [0 0 0.21 0.21],  [.8 .8 .8], edgecolor = 'none')
patch([x_5_int2(2,80:end) fliplr(x_5_int2(2,80:end))], [(x_5_int2(6,80:end))/1000 fliplr(x_5_int2(7,80:end))/1000], [0.97 0.71 0.64], EdgeColor = 'none')
plot(x_5_int2(2,:), x_5_int2(5,:)/1000, color = 'r', LineWidth= 1.2)
plot(x_5_int2(2,:), x_5_int2(6,:)/1000, Color='r', LineWidth= 0.8)
plot(x_5_int2(2,:), x_5_int2(7,:)/1000, Color='r', LineWidth= 0.8)
xlim([44.4166 47.2651])
ylim([0 0.13])
ylabel('\alpha (1/yr) ')
xlabel('Time (kyr b2k)')

%% Pre function stuff
% stadial 1
a_std = find(round(compiled_t_ig(1,:),4) == 20.0644);
b_std = find(round(compiled_t_ig(1,:),4) == 16.7931);
STD = compiled_t_ig(:, a_std:b_std);

%% Function stuff
% 5 year sampling
[x_5_s, W_5_s] = bry_func(STD, 0.005, 80);

%% proper plots stadial

figure
title('Stadial')
subplot(3,1,1)
set(gca(), ...
    'Layer','top')
set( gca, 'xdir', 'reverse' )
hold on
patch([x_5_s(2,1) x_5_s(2,80) x_5_s(2,80) x_5_s(2,1)], [-48 -48 -46.3 -46.3],  [.8 .8 .8], edgecolor = 'none')
plot(STD(1,:), STD(2,:))
xlim([16.8021 20.0629])
ylim([-48 -46.3])
ylabel('Simulated I_g')
%xlabel('Time (Kyrs)')

subplot(3,1,2)
set( gca, 'xdir', 'reverse' )
set(gca(), ...
    'Layer','top')
hold on
patch([x_5_s(2,1) x_5_s(2,80) x_5_s(2,80) x_5_s(2,1)], [-1 -1 1.5 1.5],  [.8 .8 .8], edgecolor = 'none')
plot(STD(1,:), y(6, a_std:b_std))
xlim([16.8021 20.0629])
ylim([0 1.2])
ylabel('C')
%xlabel('Time (Kyrs)')

subplot(3,1,3)
set( gca, 'xdir', 'reverse' )
set(gca(), ...
    'Layer','top')
hold on
patch([x_5_s(2,1) x_5_s(2,80) x_5_s(2,80) x_5_s(2,1)], [0 0 0.1 0.1],  [.8 .8 .8], edgecolor = 'none')
patch([x_5_s(2,80:end) fliplr(x_5_s(2,80:end))], [(x_5_s(6,80:end))/1000 fliplr(x_5_s(7,80:end))/1000], [0.97 0.71 0.64], EdgeColor = 'none')
plot(x_5_s(2,:), x_5_s(5,:)/1000, color = 'r', LineWidth= 1.2)
plot(x_5_s(2,:), x_5_s(6,:)/1000, Color='r', LineWidth= 0.8)
plot(x_5_s(2,:), x_5_s(7,:)/1000, Color='r', LineWidth= 0.8)
xlim([16.8021 20.0629])
ylabel('\alpha (1/yr) ')
xlabel('Time (kyr b2k)')

%% stadial 2
a_std2 = find(round(compiled_t_ig(1,:),4) == 29.15);
b_std2 = find(round(compiled_t_ig(1,:),4) == 23.7155);
std2 = compiled_t_ig(:, a_std2:b_std2);

[x_5_std2, W_5_std2] = bry_func(std2, 0.005, 80);

figure
subplot(3,1,1)
set( gca, 'xdir', 'reverse' )
set(gca,'FontSize',10, 'FontName', 'Outfit')
set(gca(), ...
    'Layer','top')
hold on
patch([x_5_std2(2,1) x_5_std2(2,80) x_5_std2(2,80) x_5_std2(2,1)], [-47.5 -47.5 -43 -43],  [.8 .8 .8], edgecolor = 'none')
plot(std2(1,:), std2(2,:))
%plot(x_5(2,:), x_5(8,:))
xlim([23.7155 29.15])
ylim([-47.5 -43])
ylabel('Simulated I_g')


subplot(3,1,2)
set( gca, 'xdir', 'reverse' )
set(gca,'FontSize',10, 'FontName', 'Outfit')
set(gca(), ...
    'Layer','top')
hold on
patch([x_5_std2(2,1) x_5_std2(2,80) x_5_std2(2,80) x_5_std2(2,1)], [0 0 1.2 1.2],  [.8 .8 .8], edgecolor = 'none')
plot(std2(1,:), y(6, a_std2:b_std2))
xlim([23.7155 29.15])
ylim([0.5 1.1])
ylabel('C')


subplot(3,1,3)
set( gca, 'xdir', 'reverse' )
set(gca,'FontSize',10, 'FontName', 'Outfit')
set(gca(), ...
    'Layer','top')
hold on
patch([x_5_std2(2,1) x_5_std2(2,80) x_5_std2(2,80) x_5_std2(2,1)], [0 0 0.21 0.21],  [.8 .8 .8], edgecolor = 'none')
patch([x_5_std2(2,80:end) fliplr(x_5_std2(2,80:end))], [(x_5_std2(6,80:end))/1000 fliplr(x_5_std2(7,80:end))/1000], [0.97 0.71 0.64], EdgeColor = 'none')
plot(x_5_std2(2,:), x_5_std2(5,:)/1000, color = 'r', LineWidth= 1.2)
plot(x_5_std2(2,:), x_5_std2(6,:)/1000, Color='r', LineWidth= 0.8)
plot(x_5_std2(2,:), x_5_std2(7,:)/1000, Color='r', LineWidth= 0.8)
xlim([23.7155 29.15])
ylim([0 0.13])
ylabel('\alpha (1/yr) ')
xlabel('Time (kyr b2k)')