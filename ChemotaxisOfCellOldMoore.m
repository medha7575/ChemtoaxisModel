clc

clear all; % clearing the memory
close all; % closing all figures



% Initiate the x, y coordinates of the cell center
cellCenter_x = 125;
cellCenter_y = 80;
r = 40; % radius of the cell
k= 0.6; % spring constant of the boundary
lembda= 0.0009; % area modulus of the complete cell
% points = zeros(ceil(2*3.14*30)+10,2); %initiating discretisation of the boundary

Area=zeros(1100,2); % Storing the time history of cell area
A_0 = 0; % variable to store the rest area of the cell when the springs are of rest length
totalPoints = 1; % total number of boundary points
bond_length=0; % current bond length of a particular boundary edge
%
angle= 359; % to have a complete shape round the circle
bound=200;
EGF = zeros(bound,bound);
run_time=10000;
foldername=strcat('Chemotaxis','for_mu_P5','egf5','Moore','step100');
mkdir(foldername);
absoluteFolderPath = foldername;

% 
% for row= 80:84
%    for col = 76:80
%     
%     
%     EGF(row,col) = 500;
%    end
% end


for theta = 0:1:angle
    rad = (pi)*(theta/180);
    x = cellCenter_x + r*cos(rad);
    y = cellCenter_y + r*sin(rad);
    
    points(totalPoints,1) = x;
    points(totalPoints,2) = y;
    
    totalPoints = totalPoints+1;
end

totalPoints = totalPoints-1;

storehouse= zeros(run_time,7);
%Variedforces=zeros(totalPoints,2);



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

    
    
    col_l=75;
    col_h=77;
    row_l=82;
    row_h=84;
    

for step = 1:run_time % Simulation step
    %Step1: Calculation of forces at all position
    %Step2: Change of pixel positions
    
    
    if step==1
for row= row_l:row_h
   for col = col_l:col_h
    EGF(row,col) = 5;
   end
end
    end


    if step<2000
    if mod(step,100)==0
        col_l=col_l-1;
        col_h=col_h-1;
        row_l=row_l-1;
        row_h=row_h-1;
for row= row_l:row_h
   for col = col_l:col_h      
    EGF(row,col) = 5;
   end
end
    end
    end
    
      
    if step==2000
        row_l=110;
        row_h=112;
        col_l=80;
        col_h=82;
for row= row_l:row_h
   for col = col_l:col_h      
    EGF(row,col) = 5;
   end
end
    end
if step>2000
    if mod(step,100)==0
        col_l=col_l-1;
        col_h=col_h-1;
        row_l=row_l+1;
        row_h=row_h+1;
for row= row_l:row_h
   for col = col_l:col_h      
    EGF(row,col) = 5;
   end
end
    end
end
    
% if step==8000
%         row_l=50;
%         row_h=52;
%         col_l=125;
%         col_h=127;
% for row= row_l:row_h
%    for col = col_l:col_h      
%     EGF(row,col) = 10;
%    end
% end
%     end
% if step>8000
%     if mod(step,150)==0
%         col_l=col_l+1;
%         col_h=col_h+1;
%         row_l=row_l-1;
%         row_h=row_h-1;
% for row= row_l:row_h
%    for col = col_l:col_h      
%     EGF(row,col) = 5;
%    end
% end
%     end
%     end
    
    % assigning EGF gradient
    
    EGF_old = EGF;
    for x = 2:bound-1
        for y = 2:bound-1 
            EGF(x,y) = (1/8)*(EGF_old(x-1,y)+EGF_old(x+1,y)+EGF_old(x,y-1)+EGF_old(x,y+1)+EGF_old(x-1,y-1)+EGF_old(x-1,y+1)+EGF_old(x+1,y-1)+EGF_old(x+1,y+1));
        end
    end

    %step  %print the step
    
    % CALCULATING CURRENT AREA OF CELL %
    pointa=points';
    A_c = polyarea(pointa(1,:),pointa(2,:));
    clear pointa;
    
    % SAVING THE CURRENT AREA IN A ARRAY%
    Area(step,1)= step;
    Area(step,2)= A_c;
    
    %--------------------%
    
        
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
        
        
      

        %% APPLYING CHEMOTAXIS
        
        fxcom=0;
        fycom=0;
        
        colNo = round(oldpoints(i,1));
        rowNo = round(oldpoints(i,2));
        
        for row= -1:1
            for col=-1:1
                
                if (not(row==0 & col==0))
                   
                chemoDiff=EGF(rowNo,colNo) - EGF(rowNo+row,colNo+col);
                
                x_dir= -col/(sqrt(row^2+col^2));
                y_dir= -row/(sqrt(row^2+col^2));

%                 [colNo rowNo col row x_dir y_dir chemoDiff]
                fxcom=fxcom+chemoDiff*x_dir;
                fycom=fycom+chemoDiff*y_dir;
                end
            end 
        end
        
        %pause(100);
        
        mu=0.5;
        
       % [step sum(sum(EGF))]
    
        
        fx_chem= mu*fxcom;
        fy_chem= mu*fycom;
        
        f_chem_mag = sqrt(fx_chem^2+fy_chem^2);
        
           [step f_chem_mag]
%         
%         f_chem_mag_sat = f_chem_mag/(f_chem_mag+0.5);
%         
%         prefactor = f_chem_mag_sat/f_chem_mag;
%       
        if(f_chem_mag > 0)
        fx_chem = fx_chem/(10*f_chem_mag);
        fy_chem = fy_chem/(10*f_chem_mag);
        end
        
        % CREATING A RANDOM NOISE FOR INITIAL PERTURBATION OF THE CELL
        
        %         if(step ==1)
        %             noise = rand(1,1)*1;
        %         else noise = 0;
        %         end
        
      
        
        % NET FORCE ACTING ON THE BEAD DUE TO COMBINATION OF ALL THE FORCES
        
        if(mod(step,100)< 10)
        fx= fip_x+ ftn_x+fx_chem;%+noise;
        fy= fip_y+ ftn_y+fy_chem;%+noise;
        
        else
        fx= fip_x+ ftn_x;
        fy= fip_y+ ftn_y;
        end    
          [colNo rowNo fip_x fip_y ftn_x ftn_y fx_chem fy_chem f_chem_mag]
        
        %fx= fx_chem;%+noise;
        %fy= fy_chem;%+noise;
        
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
    
    step;
    
    
    
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
    
      
    displacement=sqrt((cent_x-cellCenter_x)^2 + (cent_y-cellCenter_y)^2);
    
    
    storehouse(step,1)=cent_x;
    storehouse(step,2)=cent_y;
    storehouse(step,3)=displacement;
    storehouse(step,4)=A_c;
    storehouse(step,5)=Fip;
     storehouse(step,6)=col_h;
    storehouse(step,7)=row_h;
    
    %   
    
    
    % PLOTTING THE CELL AFTER EACH STEP
    points(totalPoints+1,1) = points(1,1);
    points(totalPoints+1,2) = points(1,2);
    points;
    
%     image(EGF*10);
%     hold on
%     plot(points(1:totalPoints,1),points(1:totalPoints,2),'r*');
% %     plot(cent_x, cent_y,'red*');
%     pause(.000001)
%     
    
    %---------------------------%
    
    % FOR SAVING PLOTS IN A FILE%
    
%     
        cd(absoluteFolderPath)
%        
if (mod(step,100)==0)
    
    filename= strcat('plot_',num2str(step));
    fig = figure;
    image(EGF*10);
    hold on
    plot(points(1:totalPoints,1),points(1:totalPoints,2),'r*');
    print(fig,filename,'-dpng');
    
    close all
end
cd ..
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

filename5= 'chemoGradient';
plot(storehouse(:,5),storehouse(:,6),'blue*');
print(fig,filename5,'-dpng');

% %
% 
cd ..
% 
% 
% cd ..
% close all
% 