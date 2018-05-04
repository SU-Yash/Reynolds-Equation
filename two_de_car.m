% Authors ~ 
    % Suyash Sardar 

% Function Calculates the following ~
    % 1.Pressure Distribution across width and length
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
% [h_bar,p_bar,Load_capacity] = two_de_car(2,20,20,1);
    
function [h_bar,p_bar,Load_capacity] = two_de_car(n,x_nodes,z_nodes,L_B)

flag =0;
iter =0;

B_L = 1/ L_B;
dx  = 1/ (x_nodes-1);
dz  = 1/ (z_nodes-1);

% Creating Mesh
x=0:dx:1;
z=0:dz:1;
[X,Z] = meshgrid(x,z);


p_bar = zeros(x_nodes,z_nodes);
h_bar = n - (n-1) * X;

while flag ~=1
    
    p_bar_prev=p_bar;
    iter = iter + 1;
    
    for i=2:z_nodes-1
        for j=2:x_nodes-1
            
            % Updating Pressure Matrix
            A=(p_bar(i+1,j)+p_bar(i-1,j))/(dx^2);
            B=(B_L^2)*(p_bar(i,j+1)+p_bar(i,j-1))/(dz^2);
            C=(n-1)/(h_bar(i,j)^3);
            D=(1.5/h_bar(i,j))*(1-n)*(p_bar(i+1,j)-p_bar(i-1,j))/dx;
            E=2*((1/dx^2)+((1/dz^2)*(B_L^2)));
            p_bar(i,j)=(A+B+C+D)/E;
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
  surf(X,Z,p_bar); 
  title([ 'PRESSURE DISTRIBUTION' '    ' 'for' '    ''Attitude Ratio:' '    ' num2str(n)])
  xlabel('Non-dimentional Length');
  ylabel('Non-dimentional Width');
  zlabel('Non-dimentional Pressure');
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Calculating Load Carrying Capacity 
  
  % Trapezoidal 2D Rule
  
  % Four Corner Points of the Meshgrid 
  Load_capacity = (p_bar(1,1) + p_bar(x_nodes,1) + p_bar(1,z_nodes) + p_bar(x_nodes,z_nodes)) ...
      * (dx * dz) / 4 ; 
  
  % Four Sides Except Corner Points of the Meshgrid
  Load_capacity = Load_capacity + (sum(p_bar(2:x_nodes-1,1)) + sum(p_bar(2:x_nodes-1,z_nodes))...
      + sum(p_bar(1,2:z_nodes-1)) + sum(p_bar(x_nodes,2:z_nodes-1))) * (dx * dz / 2) ;
  
  % Central Points (i.e : All points except the 4 sides of the Meshgrid) 
  Load_capacity = Load_capacity + (sum(sum(p_bar(2:x_nodes-1, 2:z_nodes-1)))) * (dx * dz);
  
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



