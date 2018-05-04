% Authors ~ 
    % Suyash Sardar 

% Function Calculates the following ~
    % 1.Pressure Distribution across thetha and radius
    % 2.Load Carrying Capacity of the bearing

% Inputs ~ 
    %[n ~ Attitude Ratio]
    %[x_nodes ~ Number of Nodes in X direction]
    %[z_nodes ~ Number of Nodes in Z direction]
    %[L_B ~ Length to Width Ratio]

% Outputs ~  
    %[ h_bar ~ Height at various nodes]
    %[ p_bar ~ Pressure at various nodes]
    %[ Load_capacity ~ Load carrying capacity of the bearing]
    
% Trial run for function
% [h_bar,p_bar,Load_capacity] = two_de_polar(2,20,20,1);
    
function [h_bar,p_bar,Load_capacity] = two_de_polar(n,thetha_nodes,r_nodes,L_B)

flag =0;
iter =0;
thetha_t = 60 * (pi/180);

B_L = 1/ L_B;
dthetha  = 1/ (thetha_nodes-1);
dr  = 1/ (r_nodes-1);

% Creating Mesh
thetha=0:dthetha:1;
r=0:dr:1;
[Thetha,R] = meshgrid(thetha,r);


p_bar = zeros(thetha_nodes,r_nodes);
%h_bar = n + n * (1 - Thetha);
h_bar = n - (n-1) * Thetha;

while flag ~=1
    
    p_bar_prev=p_bar;
    iter = iter + 1;
    
    for i=2:r_nodes-1
        for j=2:thetha_nodes-1
            
            % Updating Pressure Matrix
            
            A=((-1.5 * (1-n) / h_bar(i,j) * R(i,j)^2 * thetha_t^2) * (p_bar(i+1,j)-p_bar(i-1,j))) / dthetha;
            B=(1/ (R(i,j)^2 * thetha_t^2) * (p_bar(i+1,j) + p_bar(i-1,j)) / (dthetha^2));
            C=(p_bar(i,j+1) - p_bar(i,j-1)) / (2 * dr * R(i,j));
            D=(p_bar(i,j+1) + p_bar(i,j-1)) / (dr^2);
            E=((n-1)/ (h_bar(i,j)^3));
            F=2*((1/(dthetha^2 * R(i,j)^2 * thetha_t^2))+((1/dr^2)));
            p_bar(i,j)=(A+B+D+E)/F;
            p_bar(i,j)=p_bar_prev(i,j)+(p_bar(i,j)-p_bar_prev(i,j))*0.9;
            
        end
    end
    
    % Checking For Convergence
    convergence= (sum(sum(p_bar - p_bar_prev))/sum(sum(p_bar)));
    sprintf("iter: %d conv: %f",iter, convergence)
    if convergence<=1e-4
             flag=1;
    end

    % Plotting  pressure distribution 
  drawnow
  surf(Thetha,R,p_bar); 
  title([ 'PRESSURE DISTRIBUTION' '    ' 'for' '    ''Attitude Ratio:' '    ' num2str(n)])
  xlabel('Non-dimentional Thetha');
  ylabel('Non-dimentional Radius');
  zlabel('Non-dimentional Pressure');
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Calculating Load Carrying Capacity 
  
  % Trapezoidal 2D Rule
  
  % Four Corner Points of the Meshgrid 
  Load_capacity = (p_bar(1,1) + p_bar(thetha_nodes,1) + p_bar(1,r_nodes) + p_bar(thetha_nodes,r_nodes)) ...
      * (dthetha * dr) / 4 ; 
  
  % Four Sides Except Corner Points of the Meshgrid
  Load_capacity = Load_capacity + (sum(p_bar(2:thetha_nodes-1,1)) + sum(p_bar(2:thetha_nodes-1,r_nodes))...
      + sum(p_bar(1,2:r_nodes-1)) + sum(p_bar(thetha_nodes,2:r_nodes-1))) * (dthetha * dr / 2) ;
  
  % Central Points (i.e : All points except the 4 sides of the Meshgrid) 
  Load_capacity = Load_capacity + (sum(sum(p_bar(2:thetha_nodes-1, 2:r_nodes-1)))) * (dthetha * dr);
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 % Displaying Results
end
disp(' ')
t_time=clock;
disp(['================================ ',date,' ================================'])
disp(['============= Steady State Analysis of Hydrodynamic Slider Bearings ============'])
disp(['================================= Time ',num2str(t_time(4)),':',num2str(t_time(5)),' ================================='])
disp('*****************************************************************************')
sprintf("Load Carrying Capacity (Non-Dimensionalized Value) : %f", Load_capacity)
disp('*****************************************************************************')