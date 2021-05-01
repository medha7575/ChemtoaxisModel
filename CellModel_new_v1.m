clc

clear all; % clearing the memory
close all; % closing all figures

% Initiate the x, y coordinates of the cell center
cellCenter_x = 100;
cellCenter_y = 100;
r = 30; % radius of the cell
k= 0.5; % spring constant of the boundary
lembda= 0.0003; % area modulus of the complete cell
% points = zeros(ceil(2*3.14*30)+10,2); %initiating discretisation of the boundary

Area=zeros(1100,2); % Storing the time history of cell area
A_0 = 0; % variable to store the rest area of the cell when the springs are of rest length
totalPoints = 1; % total number of boundary points
bond_length=0; % current bond length of a particular boundary edge
%
angle= 359; % to have a complete shape round the circle
 force_time= 105;  % time for which internal forces are applied
force=0.2;   % internal forces acting at certain part of cell
%  force_n=2;      % for creating new folder to store data
 foldername=strcat( 'force_B','force_n02','_for_',num2str(force_time),'steps');

run_time=25000;
mkdir('NEWresults',foldername);

absoluteFolderPath = strcat('NEWresults\', foldername);

%  Kon= 0.2; % recruitment rate of FBP17 molecules
% Koff= 1.5; % detahcmnet constant for FBP17 molecules
% fpro_0= 0.2; %protusion constant for actin protusive forces
% foldername=strcat('tensionffect','lembda00009','k0001','kon02koff1point5','fbpPOINTS9');
% mkdir('Simulation_result',foldername);
% absoluteFolderPath = strcat('Simulation_result\', foldername);

for theta = 0:1:angle'
    rad = (pi)*(theta/180);
    x = cellCenter_x + r*cos(rad);
    y = cellCenter_y + r*sin(rad);
    
    points(totalPoints,1) = x;
    points(totalPoints,2) = y;
    
    totalPoints = totalPoints+1;
end

totalPoints = totalPoints-1;

storehouse= zeros(run_time,7);
Variedforces=zeros(totalPoints,2);



% CALCULATING REST TANGENT LENGTH OF CELL %
for i= 1:totalPoints
    next_i = mod(i,totalPoints) + 1;
    bond_length = bond_length + sqrt((points(i,1)-points(next_i,1))^2+(points(i,2)-points(next_i,2))^2);
end
tn_rest= bond_length/totalPoints;

% CALCULATING REST AREA OF CELL %
pointa=points';
A_0 = 1*polyarea(pointa(1,:),pointa(2,:));
clear pointa;
%------------------%


%SAVING INITIAL SHAPE OF THE CELL%
% filename= strcat('plot_','initial');
% fig = figure;
% plot(points(1:totalPoints,1),points(1:totalPoints,2),'blue');
% %print(fig,filename,'-dpng')
% close all

% STORING LATEST POINTS IN OLD-ARRAY%
oldpoints = points;


%   TotalFBP=0;
%     oldfbp= zeros(totalPoints, 1);
%      for i= 1: 20
%         oldfbp(i)= 4;
%     end
%     for i= 60:70
%         oldfbp(i)= 2;
%     end
%     
%     for i= 80:90
%         oldfbp(i)= 2;
%     end
%     for i= 320:330
%         oldfbp(i)= 2;
%     end
%     
%     for i= 340:350
%         oldfbp(i)= 2;
%     end
%     for i= 110:120
%         oldfbp(i)= 1;
%     end
%     for i= 130:140
%         oldfbp(i)= 1;
%     end
%     for i= 190:200
%         oldfbp(i)= 1;
%     end 
%     for i= 210:220
%         oldfbp(i)= 1;
%     end  

for step = 1:run_time % Simulation step
    %Step1: Calculation of forces at all position
    %Step2: Change of pixel positions
    
    
%     step  %print the step
    
    % CALCULATING CURRENT AREA OF CELL %
    pointa=points';
    A_c = polyarea(pointa(1,:),pointa(2,:));
    clear pointa;
    
    % SAVING THE CURRENT AREA IN A ARRAY%
    Area(step,1)= step;
    Area(step,2)= A_c;
    
    %--------------------%
    
    % initial distribution of fbp
%     TotalFBP=0;
%     oldfbp= zeros(totalPoints, 1); % defining intial fbp molecules at each bead
    %     for TotalFBP= 1: 200
    %         i= round((totalPoints-1)*rand(1,1))+ 1;
    %         oldfbp(i)= oldfbp(i)+1;
    %     end
   
    
    for i  = 1: totalPoints
        
        % CALCULATING TANGENTIAL FORCE ACTING ON THE CELL%
        
        previous_i= mod((i+ totalPoints-2),totalPoints)+1;
        next_i = mod(i,totalPoints) + 1;
        
        % length of the tangent towards the right
        tn_rx= oldpoints(next_i,1) - oldpoints(i,1);
        tn_ry= oldpoints(next_i,2) - oldpoints(i,2);
        tn_r = sqrt((tn_rx)^2 + (tn_ry)^2);
        
        %       Forces on the right of ith point due to tension in spring
        
        Ftn_rx= k*(tn_r- tn_rest)*(tn_rx/tn_r);
        Ftn_ry= k*(tn_r- tn_rest)*(tn_ry/tn_r);
        %
        
        % length of the tangent towards the left
        tn_lx= oldpoints(i,1) - oldpoints(previous_i,1);
        tn_ly= oldpoints(i,2) - oldpoints(previous_i,2);
        tn_l = sqrt((tn_lx)^2 + (tn_ly)^2);
        
        % Forces on the left of ith point due to tension in spring
        Ftn_lx= k*(tn_l- tn_rest)*(tn_lx/tn_l);
        Ftn_ly= k*(tn_l- tn_rest)*(tn_ly/tn_l);
        %
        %Net force acting on the ith bead due to tension
        ftn_x = (Ftn_rx - Ftn_lx);
        ftn_y = (Ftn_ry - Ftn_ly);
        
        
        %  [Ftn_rx Ftn_ry Ftn_lx Ftn_ly ftn_x ftn_y]
        
        % CALCULATING INTERNAL PRESSURE FORCE ACTING ON THE CELL%
        Fip = lembda*(A_0-A_c);
        %---------------------%
        
        %         %check for differential application of force
        
        
        if (i >=170 && i <=190 && step <= force_time)
            Fip= Fip+force;
        end
        %
        %CALCULATING NORMAL USING TANGENTS AT THE BEAD
        
        %    Area_FIP=         [A_0 A_c Fip]
        
        % unit length of the tangent
        tn_rx_unit= (tn_rx/tn_r);
        tn_ry_unit= (tn_ry/tn_r);
        
        tn_lx_unit= (tn_lx/tn_l);
        tn_ly_unit= (tn_ly/tn_l);
        
        %        Normal at the bead
        Nx= -(tn_rx_unit-tn_lx_unit);
        Ny= -(tn_ry_unit-tn_ly_unit);
        
        %CALCULATING NORMALIZED UNIT NORMAL VECTOR FOR APPLICATION OF FORCE
        
        norm = sqrt(Nx^2+Ny^2);
        
        Nx_unit = Nx/norm;
        Ny_unit = Ny/norm;
        
        %Force acting on the bead due to internal pressure
        fip_x = Fip*Nx_unit;
        fip_y = Fip*Ny_unit;
        
        
        %                   FIP=  [fip_x fip_y]
        
        % CALCULATING ACTIN PROTUSIVE FORCE ACTING ON THE CELL DUE TO FBP17
        
        %
%         ftn = sqrt(ftn_x^2+ftn_y^2);  % absolute tension at the ith bead
%         
%         
%         %         if (i >=0 && i <=45 && step <= force_time)
%         %             ftn= ftn-force;
%         %         end
%         %
%         
%         fbp(i)= oldfbp(i) + Kon+ oldfbp(i)/(oldfbp(i)+0.25) -Koff*(ftn/(1+ftn));
%         %         fbp(i)= oldfbp(i) + Kon -Koff*(ftn/(1+ftn));
%         
%         if fbp(i) < 0
%             fbp(i)=0;
%         end
%         
%         
%         fpro= fpro_0*(fbp(i)/(1+fbp(i)));
%         
%         fpro_x=fpro*Nx_unit;
%         fpro_y=fpro*Ny_unit;
%         NewFBPtotal=    TotalFBP+fbp(i)-oldfbp(i);
        %
        %
        
%         Variedforces(i,1)=ftn;
%         Variedforces(i,2)=fpro;
        
        % CREATING A RANDOM NOISE FOR INITIAL PERTURBATION OF THE CELL
        
        %         if(step ==1)
        %             noise = rand(1,1)*1;
        %         else noise = 0;
        %         end
        
        % NET FORCE ACTING ON THE BEAD DUE TO COMBINATION OF ALL THE FORCES
        fx= fip_x+ ftn_x;%+fpro_x;%+noise;
        fy= fip_y+ ftn_y;%+fpro_y;%+noise;
        
        % % TO ADD NOISE AT EVERY POINT AT ALL STEPS
        %         fx = fx+ 0.05*(-1+rand(1,1)*2)*fx;
        %         fy = fy+ 0.05*(-1+rand(1,1)*2)*fy;
        % %
        
        
        % NEW POSITION OF THE BEAD
        xnew = points(i,1)+fx;
        ynew = points(i,2)+fy;
        %
        % UPDATING THE NEW POINTS IN THE ARRAY
        points(i,1) = xnew;
        points(i,2) = ynew;
        %
        %FOR PLOTTING EVERY POINT CHANGE WITHIN A TIME STEP
        %         plot(points(1:totalPoints,1),points(1:totalPoints,2),'blue');
        %         pause(0.2)
        %         clf
       
        
    end
    
    % SAVING THE NEWLY GENERATED POINTS IN OLD_POINTS FOR THE NEXT STEP
    for kount  = 1 : totalPoints
        oldpoints(kount,1) = points(kount,1);
        oldpoints(kount,2) = points(kount,2);
    end
    
    %CALCULATING CENTROID OF THE CELL
    x_sum=0;
    y_sum=0;
    for p= 1: totalPoints
        x_sum = x_sum + oldpoints(p,1);
        y_sum = y_sum + oldpoints(p,2);
        cent_x= x_sum/totalPoints;
        cent_y= y_sum/totalPoints;
    end
    
    
    
%     if (mod(step,1000)==0)
     
%     cd(absoluteFolderPath)
%     
%     filename5= strcat('Ftn_variation',num2str(step));
%            fig = figure;
%     plot(points(1:360,1), Variedforces(:,1),'blue*');
%     print(fig,filename5,'-dpng');
%     
%     
%      filename6= strcat('Fpro_variation',num2str(step));
%            fig = figure;
%     plot(points(1:360,1), Variedforces(:,2),'blue*');
%     print(fig,filename6,'-dpng');
%     
%     cd ..
%     cd ..
%     end
% %     
    
    displacement=sqrt((cent_x-cellCenter_x)^2 + (cent_y-cellCenter_y)^2);
    
    
    storehouse(step,1)=cent_x;
    storehouse(step,2)=cent_y;
    storehouse(step,3)=displacement;
    storehouse(step,4)=A_c;
    storehouse(step,5)=Fip;
    %     storehouse(step,6)=fpro;
    %     storehouse(step,7)= NewFBPtotal;
    %     storehouse(step,6)=force_time;
    % %     storehouse(step,7)=force;
    
    
    % PLOTTING THE CELL AFTER EACH STEP
    points(totalPoints+1,1) = points(1,1);
    points(totalPoints+1,2) = points(1,2);
    points;
    %     if(mod(step,100)==0)
%     plot(points(1:totalPoints,1),points(1:totalPoints,2),'blue*');
%     hold on
%     plot(cent_x, cent_y,'red*');
%     pause(0.00001)
%     clf
%     %     end
    %---------------------------%
    
    % FOR SAVING PLOTS IN A FILE%
    
    
    cd(absoluteFolderPath)
    if (mod(step,10)==0 && step<= 100)
        
        filename= strcat('plot_',num2str(step));
        fig = figure;
        plot(points(1:totalPoints,1),points(1:totalPoints,2),'blue*');
        hold on
        plot(cent_x, cent_y,'red*');
        print(fig,filename,'-dpng');
        close all
        
    end
%     if (mod(step,25)==0 && step<= 1000)
%         
%         filename= strcat('plot_',num2str(step));
%         fig = figure;
%         plot(points(1:totalPoints,1),points(1:totalPoints,2),'blue*');
%         hold on
%         plot(cent_x, cent_y,'red*');
%         print(fig,filename,'-dpng');
%         close all
%     end
    if (mod(step,50)==0 && step<= 100000)
        
        filename= strcat('plot_',num2str(step));
        fig = figure;
        plot(points(1:totalPoints,1),points(1:totalPoints,2),'blue*');
        hold on
        plot(cent_x, cent_y,'red*');
        print(fig,filename,'-dpng');
        
        close all
    end
    cd ..
    cd ..
    %     %
        %     %                         step
        %     %                         Area_FIP=         [A_0 A_c Fip]
        %     %                         [ ftn_x ftn_y fip_x fip_y ]%fpro_x fpro_y]
        %     %                         Center= [cent_x,cent_y]
        %     %
        %
        %
end
cd (absoluteFolderPath)

filename1= 'cell_movement';
fig = figure;
plot(storehouse(:,1),storehouse(:,2),'blue*');
print(fig,filename1,'-dpng');

filename2= 'final_displacement';
plot(storehouse(:,3),'blue*');
print(fig,filename2,'-dpng');

filename3= 'Area_variation';
plot(storehouse(:,4),'blue*');
print(fig,filename3,'-dpng');

filename4= 'Fip_variation';
plot(storehouse(:,5),'blue*');
print(fig,filename4,'-dpng');

% filename5= 'Fpro_variation';
% plot(storehouse(:,6),'blue*');
% print(fig,filename5,'-dpng');
%
% filename6= 'fbp';
% plot(storehouse(:,7),'blue*');
% print(fig,filename6,'-dpng');
%

cd ..


cd ..
close all

