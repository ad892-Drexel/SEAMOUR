%% Run this file to test the SEAMOUR Model
% Last updat 08/13/2024
clear
clc
close all
%%
run('init_file.m')
worker_model='SEAMOUR_Sim';
open(worker_model)
set_param(worker_model,'FastRestart','off')
set_param(worker_model, 'SimulationMode', 'Normal');
set_param(worker_model,'SimMechanicsOpenEditorOnUpdate','on')
%%

A0_color=[228,26,28]/255;
A1_color=[55,126,184]/255;
A2_color=[77,175,74]/255;
A3_color=[152,78,163]/255;

%%
for j=1:4

    if j==1
        clc 
        load('A0_Pitch_Char.mat')
        load('A0_Yaw_Char.mat')
        load('A0_Roll_Char.mat')
        disp('Executing Characteristic Stroke...')
    end

    if j==2
        clc 
        load('A1_Pitch_Learned.mat')
        load('A1_Yaw_Learned.mat')
        load('A1_Roll_Learned.mat')
        disp('Executing Learned Stroke A...')
    end

    if j==3
        clc    
        load('A2_Pitch_Learned.mat')
        load('A2_Yaw_Learned.mat')
        load('A2_Roll_Learned.mat')
        disp('Executing Learned Stroke B...')
    end

    if j==4

        clc 
        load('A3_Pitch_Learned.mat')
        load('A3_Yaw_Learned.mat')
        load('A3_Roll_Learned.mat')
        disp('Executing Learned Stroke C...')

    end


    %%
    cycles=4;
    for i=1:cycles-1
        Pitch_Traj=[Pitch_Traj, Pitch_Traj];
        Yaw_Traj=[Yaw_Traj, Yaw_Traj];
        Roll_Traj=[Roll_Traj, Roll_Traj];
    end



    %% Time Vector
    Time=simulation_ts:simulation_ts:length(Roll_Traj)*simulation_ts;

    %% Send Trajectories to Motors in Simscape

    fore_pitch_left=linspace(0,0,10);
    assignin('base','fore_pitch_left',[Time; -Pitch_Traj]');
    fore_yaw_left=linspace(0,0,10);
    assignin('base','fore_yaw_left',[Time;   Yaw_Traj]');
    fore_roll_left=linspace(0,0,10);
    assignin('base','fore_roll_left',[Time;  -Roll_Traj]');

    fore_pitch_right=linspace(0,0,10);
    assignin('base','fore_pitch_right',[Time; Pitch_Traj]');
    fore_yaw_right=linspace(0,0,10);
    assignin('base','fore_yaw_right',[Time;   -Yaw_Traj]');
    fore_roll_right=linspace(0,0,10);
    assignin('base','fore_roll_right',[Time;  Roll_Traj]');

    %% Simulate
    warning('off', 'Simulink:Engine:BlockDiagramChangedOnDisk'); % Removes annoying meaningless error signals
    out=sim(worker_model);
    warning('off', 'Simulink:Engine:BlockDiagramChangedOnDisk'); % Removes annoying meaningless error signals
    %%
    if j==1
        TrialData_A0=out.Out.signals.values;
    end

    if j==2
        TrialData_A1=out.Out.signals.values;
    end

    if j==3
        TrialData_A2=out.Out.signals.values;
    end

    if j==4
        TrialData_A3=out.Out.signals.values;
    end

end

%%
% figure
% subplot(3,1,1)
% plot(TrialData(:,1)) %x
% subplot(3,1,2)
% plot(TrialData(:,2)) %y
% subplot(3,1,3)
% plot(TrialData(:,3)) %z
% 
ts=0.001;
time=ts:ts:length(TrialData_A0(:,8))*ts;

figure('Position',[100 100 1300 300])
hold on
plot(TrialData_A0(:,1),TrialData_A0(:,3),'-','Color',A0_color,'LineWidth',3)
plot(TrialData_A1(:,1),TrialData_A1(:,3),'--','Color',A1_color,'LineWidth',3)
plot(TrialData_A2(:,1),TrialData_A2(:,3),':','Color',A2_color,'LineWidth',3)
plot(TrialData_A3(:,1),TrialData_A3(:,3),'-.','Color',A3_color,'LineWidth',3)
hold off

xlim([0,2.55])
ylim([-1.5,1.5])
yticks(-1.5:.75:1.5)
fontsize(18,"points");
fontname('Times New Roman')


%%
figure('Position',[100 100 1300 300])
hold on
plot(time,-rad2deg(TrialData_A0(:,8)),'-','Color',A0_color,'LineWidth',3)
plot(time,-rad2deg(TrialData_A1(:,8)),'--','Color',A1_color,'LineWidth',3)
plot(time,-rad2deg(TrialData_A2(:,8)),':','Color',A2_color,'LineWidth',3)
plot(time,-rad2deg(TrialData_A3(:,8)),'-.','Color',A3_color,'LineWidth',3)

xline(3,'k--','LineWidth',1)
xline(6,'k--','LineWidth',1)
xline(9,'k--','LineWidth',1)

ylim([-30 70])
xlim([0 12])
% xlabel('Time (s)','FontAngle','italic')
% ylabel('Pitch (^o)','FontAngle','italic')
% title('Pitch Change')
hold off

ylim([-30 70])
xlim([0 12])
xticks(0:3:12)
fontsize(18,"points");
fontname('Times New Roman')
%legend({'Char','A','B','C'},Location="northwest",NumColumns=2)

ind_stroke_A0=[];
ind_stroke_A1=[];
ind_stroke_A2=[];
ind_stroke_A3=[];

for i=1:4
    strokes=-TrialData_A0(((i-1)*3000)+1:i*3000,8);
    ind_stroke_A0=[ind_stroke_A0 strokes];

    strokes=-TrialData_A1(((i-1)*3000)+1:i*3000,8);
    ind_stroke_A1=[ind_stroke_A1 strokes];

    strokes=-TrialData_A2(((i-1)*3000)+1:i*3000,8);
    ind_stroke_A2=[ind_stroke_A2 strokes];

    strokes=-TrialData_A3(((i-1)*3000)+1:i*3000,8);
    ind_stroke_A3=[ind_stroke_A3 strokes];
   
end

deltas_A0=[ind_stroke_A0(end,1)-ind_stroke_A0(1,1),...
           ind_stroke_A0(end,2)-ind_stroke_A0(1,2),...
           ind_stroke_A0(end,3)-ind_stroke_A0(1,3),...
           ind_stroke_A0(end,4)-ind_stroke_A0(1,4)];

deltas_A1=[ind_stroke_A1(end,1)-ind_stroke_A1(1,1),...
           ind_stroke_A1(end,2)-ind_stroke_A1(1,2),...
           ind_stroke_A1(end,3)-ind_stroke_A1(1,3),...
           ind_stroke_A1(end,4)-ind_stroke_A1(1,4)];

deltas_A2=[ind_stroke_A2(end,1)-ind_stroke_A2(1,1),...
           ind_stroke_A2(end,2)-ind_stroke_A2(1,2),...
           ind_stroke_A2(end,3)-ind_stroke_A2(1,3),...
           ind_stroke_A2(end,4)-ind_stroke_A2(1,4)];

deltas_A3=[ind_stroke_A3(end,1)-ind_stroke_A3(1,1),...
           ind_stroke_A3(end,2)-ind_stroke_A3(1,2),...
           ind_stroke_A3(end,3)-ind_stroke_A3(1,3),...
           ind_stroke_A3(end,4)-ind_stroke_A3(1,4)];

disp(deltas_A0)
disp(deltas_A1)
disp(deltas_A2)
disp(deltas_A3)

%%

figure('Position',[100 100 1300 300])
hold on
plot(time,-rad2deg(TrialData_A0(:,11)),'-','Color',A0_color,'LineWidth',3)
plot(time,-rad2deg(TrialData_A1(:,11)),'--','Color',A1_color,'LineWidth',3)
plot(time,-rad2deg(TrialData_A2(:,11)),':','Color',A2_color,'LineWidth',3)
plot(time,-rad2deg(TrialData_A3(:,11)),'-.','Color',A3_color,'LineWidth',3)

xline(3,'k--','LineWidth',1)
xline(6,'k--','LineWidth',1)
xline(9,'k--','LineWidth',1)

ylim([-5 15])
xlim([0 12])
% xlabel('Time (s)','FontAngle','italic')
% ylabel('Pitch (^o)','FontAngle','italic')
% title('Pitch Change')
hold off

ylim([-5 15])
xlim([0 12])
xticks(0:3:12)
fontsize(18,"points");
fontname('Times New Roman')
%%
legend({'Char','A','B','C'},Location="northwest",NumColumns=2)

ind_stroke_A0_vel=[];
ind_stroke_A1_vel=[];
ind_stroke_A2_vel=[];
ind_stroke_A3_vel=[];

for i=1:4
    strokes=-TrialData_A0(((i-1)*3000)+1:i*3000,11);
    ind_stroke_A0_vel=[ind_stroke_A0_vel strokes];

    strokes=-TrialData_A1(((i-1)*3000)+1:i*3000,11);
    ind_stroke_A1_vel=[ind_stroke_A1_vel strokes];

    strokes=-TrialData_A2(((i-1)*3000)+1:i*3000,11);
    ind_stroke_A2_vel=[ind_stroke_A2_vel strokes];

    strokes=-TrialData_A3(((i-1)*3000)+1:i*3000,11);
    ind_stroke_A3_vel=[ind_stroke_A3_vel strokes];
   
end


%%

%%
legend({'Char','A','B','C'},Location="northwest",NumColumns=2)

% ind_stroke_A0_vel=[];
% ind_stroke_A1_vel=[];
% ind_stroke_A2_vel=[];
% ind_stroke_A3_vel=[];
% 
% for i=1:4
%     strokes=-TrialData_A0(((i-1)*3000)+1:i*3000,4);
%     ind_stroke_A0_vel=[ind_stroke_A0_vel strokes];
% 
%     strokes=-TrialData_A1(((i-1)*3000)+1:i*3000,4);
%     ind_stroke_A1_vel=[ind_stroke_A1_vel strokes];
% 
%     strokes=-TrialData_A2(((i-1)*3000)+1:i*3000,4);
%     ind_stroke_A2_vel=[ind_stroke_A2_vel strokes];
% 
%     strokes=-TrialData_A3(((i-1)*3000)+1:i*3000,4);
%     ind_stroke_A3_vel=[ind_stroke_A3_vel strokes];
% 
% end
% 
% disp(mean(ind_stroke_A0_vel))
% disp(mean(ind_stroke_A1_vel))
% disp(mean(ind_stroke_A2_vel))
% disp(mean(ind_stroke_A3_vel))
% 
% reshaped_matrix = reshape((TrialData_A0(1:end-1,4)), [], 4);
% mean_cycle=mean(sum(reshaped_matrix')/4);
% disp('A0 Angular Pitch Change Per Cycle')
% disp(mean_cycle)
% 
% reshaped_matrix = reshape((TrialData_A1(1:end-1,4)), [], 4);
% mean_cycle=mean(sum(reshaped_matrix')/4);
% disp('A1 Angular Pitch Change Per Cycle')
% disp(mean_cycle)
% 
% reshaped_matrix = reshape((TrialData_A2(1:end-1,4)), [], 4);
% mean_cycle=mean(sum(reshaped_matrix')/4);
% disp('A2 Angular Pitch Change Per Cycle')
% disp(mean_cycle)
% 
% reshaped_matrix = reshape((TrialData_A3(1:end-1,4)), [], 4);
% mean_cycle=mean(sum(reshaped_matrix')/4);
% disp('A3 Angular Pitch Change Per Cycle')
% disp(mean_cycle)
% 
% disp(mean(ind_stroke_A0_vel))
% disp(mean(ind_stroke_A1_vel))
% disp(mean(ind_stroke_A2_vel))
% disp(mean(ind_stroke_A3_vel))
% 
% reshaped_matrix = reshape(-rad2deg(TrialData_A0(1:end-1,11)), [], 4);
% mean_cycle=mean(sum(reshaped_matrix')/4);
% disp('A0 Angular Pitch Change Per Cycle')
% disp(mean_cycle)
% 
% reshaped_matrix = reshape(-rad2deg(TrialData_A1(1:end-1,11)), [], 4);
% mean_cycle=mean(sum(reshaped_matrix')/4);
% disp('A1 Angular Pitch Change Per Cycle')
% disp(mean_cycle)
% 
% reshaped_matrix = reshape(-rad2deg(TrialData_A2(1:end-1,11)), [], 4);
% mean_cycle=mean(sum(reshaped_matrix')/4);
% disp('A2 Angular Pitch Change Per Cycle')
% disp(mean_cycle)
% 
% reshaped_matrix = reshape(-rad2deg(TrialData_A3(1:end-1,11)), [], 4);
% mean_cycle=mean(sum(reshaped_matrix')/4);
% disp('A3 Angular Pitch Change Per Cycle')
% disp(mean_cycle)

%%
% 
% 
% LS=2;
% 
% figure('Position',[200 0 450 300])
% hold on
% plot(time,(TrialData_A0(:,4)),'-','Color',A0_color,'LineWidth',LS)
% plot(time,(TrialData_A1(:,4)),'--','Color',A1_color,'LineWidth',LS)
% plot(time,(TrialData_A2(:,4)),':','Color',A2_color,'LineWidth',LS)
% plot(time,(TrialData_A3(:,4)),'-.','Color',A3_color,'LineWidth',LS)
% 
% xline(3,'k--','LineWidth',1)
% xline(6,'k--','LineWidth',1)
% xline(9,'k--','LineWidth',1)
% 
% hold off
% % xlabel('Time (s)','FontAngle','italic')
% % ylabel(x_dot_string,'Interpreter','latex','FontAngle','italic')
% % title('Surge Velocity')
% fontsize(12,"points");
% fontname('Times New Roman')
% 
% ylim([-.25 .75])
% yticks(-.25:.25:.75)
% xlim([0 12])
% xticks(0:3:12)
% %%
% 
% LS=2;
% 
% figure('Position',[200 0 450 300])
% hold on
% plot(time,(TrialData_A0(:,6)),'-','Color',A0_color,'LineWidth',LS)
% plot(time,(TrialData_A1(:,6)),'--','Color',A1_color,'LineWidth',LS)
% plot(time,(TrialData_A2(:,6)),':','Color',A2_color,'LineWidth',LS)
% plot(time,(TrialData_A3(:,6)),'-.','Color',A3_color,'LineWidth',LS)
% 
% xline(3,'k--','LineWidth',1)
% xline(6,'k--','LineWidth',1)
% xline(9,'k--','LineWidth',1)
% 
% hold off
% %legend({'Char','A1','A2','A3'},'Location','NorthWest')
% % xlabel('Time (s)','FontAngle','italic')
% % ylabel(z_dot_string,'Interpreter','latex','FontAngle','italic')
% % title('Heave Velocity')
% 
% fontsize(12,"points");
% fontname('Times New Roman')
% 
% ylim([-.25 .75])
% yticks(-.25:.25:.75)
% xlim([0 12])
% xticks(0:3:12)
% 
% %%
% LS=2;
% 
% figure('Position',[200 0 450 300])
% hold on
% plot(time,sqrt(TrialData_A0(:,6).^2+TrialData_A0(:,4).^2),'-','Color',A0_color,'LineWidth',LS)
% plot(time,sqrt(TrialData_A1(:,6).^2+TrialData_A0(:,4).^2),'--','Color',A1_color,'LineWidth',LS)
% plot(time,sqrt(TrialData_A2(:,6).^2+TrialData_A0(:,4).^2),':','Color',A2_color,'LineWidth',LS)
% plot(time,sqrt(TrialData_A3(:,6).^2+TrialData_A0(:,4).^2),'-.','Color',A3_color,'LineWidth',LS)
% 
% xline(3,'k--','LineWidth',1)
% xline(6,'k--','LineWidth',1)
% xline(9,'k--','LineWidth',1)
% 
% hold off
% %legend({'Char','A1','A2','A3'},'Location','NorthWest')
% % xlabel('Time (s)','FontAngle','italic')
% % ylabel('Velocity (m/s)','FontAngle','italic')
% % title('Combined Velocity')
% 
% fontsize(12,"points");
% fontname('Times New Roman')
% 
% ylim([-.25 .75])
% yticks(-.25:.25:.75)
% xlim([0 12])
% xticks(0:3:12)
% 
% %% Print Results
% 
% % Assuming velx and velz are defined as follows
% velx = TrialData_A0(:,4); % Extracting the 4th column for velx
% velz = TrialData_A0(:,6); % Extracting the 6th column for velz
% 
% % Calculating the mean and maximum values for velx and velz
% mean_velx = mean(velx);
% max_velx = max(velx);
% mean_velz = mean(velz);
% max_velz = max(velz);
% 
% % Calculating the magnitude of each pair (velx, velz)
% magnitudes = sqrt(velx.^2 + velz.^2);
% 
% % Calculating mean and maximum of the magnitudes
% mean_magnitude = mean(magnitudes);
% max_magnitude = max(magnitudes);
% 
% % Displaying the results using fprintf
% disp('----------A0----------')
% fprintf('Mean of velx: %f\nMaximum of velx: %f\nMean of velz: %f\nMaximum of velz: %f\nMean magnitude: %f\nMaximum magnitude: %f\n', ...
%         mean_velx, max_velx, mean_velz, max_velz, mean_magnitude, max_magnitude);
% 
% % Assuming velx and velz are defined as follows
% velx = TrialData_A1(:,4); % Extracting the 4th column for velx
% velz = TrialData_A1(:,6); % Extracting the 6th column for velz
% 
% % Calculating the mean and maximum values for velx and velz
% mean_velx = mean(velx);
% max_velx = max(velx);
% mean_velz = mean(velz);
% max_velz = max(velz);
% 
% % Calculating the magnitude of each pair (velx, velz)
% magnitudes = sqrt(velx.^2 + velz.^2);
% 
% % Calculating mean and maximum of the magnitudes
% mean_magnitude = mean(magnitudes);
% max_magnitude = max(magnitudes);
% 
% % Displaying the results using fprintf
% disp('----------A1----------')
% fprintf('Mean of velx: %f\nMaximum of velx: %f\nMean of velz: %f\nMaximum of velz: %f\nMean magnitude: %f\nMaximum magnitude: %f\n', ...
%         mean_velx, max_velx, mean_velz, max_velz, mean_magnitude, max_magnitude);
% 
% 
% % Assuming velx and velz are defined as follows
% velx = TrialData_A2(:,4); % Extracting the 4th column for velx
% velz = TrialData_A2(:,6); % Extracting the 6th column for velz
% 
% % Calculating the mean and maximum values for velx and velz
% mean_velx = mean(velx);
% max_velx = max(velx);
% mean_velz = mean(velz);
% max_velz = max(velz);
% 
% % Calculating the magnitude of each pair (velx, velz)
% magnitudes = sqrt(velx.^2 + velz.^2);
% 
% % Calculating mean and maximum of the magnitudes
% mean_magnitude = mean(magnitudes);
% max_magnitude = max(magnitudes);
% 
% % Displaying the results using fprintf
% disp('----------A2----------')
% fprintf('Mean of velx: %f\nMaximum of velx: %f\nMean of velz: %f\nMaximum of velz: %f\nMean magnitude: %f\nMaximum magnitude: %f\n', ...
%         mean_velx, max_velx, mean_velz, max_velz, mean_magnitude, max_magnitude);
% 
% % Assuming velx and velz are defined as follows
% velx = TrialData_A3(:,4); % Extracting the 4th column for velx
% velz = TrialData_A3(:,6); % Extracting the 6th column for velz
% 
% % Calculating the mean and maximum values for velx and velz
% mean_velx = mean(velx);
% max_velx = max(velx);
% mean_velz = mean(velz);
% max_velz = max(abs(velz));
% 
% % Calculating the magnitude of each pair (velx, velz)
% magnitudes = sqrt(velx.^2 + velz.^2);
% 
% % Calculating mean and maximum of the magnitudes
% mean_magnitude = mean(magnitudes);
% max_magnitude = max(magnitudes);
% 
% % Displaying the results using fprintf
% disp('----------A3----------')
% fprintf('Mean of velx: %f\nMaximum of velx: %f\nMean of velz: %f\nMaximum of velz: %f\nMean magnitude: %f\nMaximum magnitude: %f\n', ...
%         mean_velx, max_velx, mean_velz, max_velz, mean_magnitude, max_magnitude);
% 
% 
% 
% 
% % 
% % figure
% % subplot(3,1,1)
% % plot(rad2deg(TrialData(:,7)))
% % subplot(3,1,2)
% % plot(-rad2deg(TrialData(:,8))) %Pitch
% % subplot(3,1,3)
% % plot(rad2deg(TrialData(:,9)))
% 
