% Authors ~ 
    % Suyash Sardar

% Function Calculates the following ~
    % 1.Pressure Distribution across width
    % 2.Load Carrying Capacity of the bearing
    % 3.Shear Stress distribution in the bearing
    % 4.Friction Force 
    % 5.Friction Coefficient

% Inputs ~ 
    %[n ~ Attitude Ratio]
    %[nodes ~ Number of Nodes]

% Outputs ~  
    %[ h_bar ~ Height at various nodes]
    %[ p_bar ~ Pressure at various nodes]
    %[ tau_bar ~ Shear stress at various nodes]
    %[ Load_capacity ~ Load carrying capacity of the bearing]
    %[ Friction_force ~ Friction force generated in the bearing]
    %[ myu ~ Friction coefficient corresponding to the given load and friction force]
    
% Trial run for function
% [h_bar,dx,p_bar,tau_bar,Load_capacity,Friction_force,myu] = one_de(2,20);
    
function [h_bar,dx,p_bar,tau_bar,Load_capacity,Friction_force,myu] = one_de(n,nodes)

x_bar = linspace(0,1,nodes); % Discritizing x-direction in nx nodes  
h_bar = n - (n-1)*x_bar;  
flag = 0; 
iter=0;
p_bar = zeros(1,nodes); % vector containing pressures at nx nodes
tau_bar = zeros(1,nodes);% vector containing shear stresses at nx nodes

% Following values are calculated to ease further computation.
dx = x_bar(2)-x_bar(1); 
h_sqr = power(h_bar,2);
h_cub = power(h_bar,3);
dx_sqr = power(dx,2);
%dx_mat = ones(nx,1)*dx; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculating the Pressure Distribution at nx nodes
while(flag~=1)
    p_bar_prev = p_bar;
    iter = iter+ 1;
    
    % updating the pressure vector
    for i = 2 : (nodes-1)
        p_bar(i) = (dx_sqr/ (2 * h_cub(i))) * (((h_cub(i)* (p_bar(i+1) + p_bar(i-1))) / dx_sqr) +...
        (((p_bar(i+1) - p_bar(i-1)) * 3 * h_sqr(i) * (1-n)) /  (2 * dx)) - (1-n));
        
    end
    
     % checking for convergence
    convergence = (abs(sum(p_bar))*dx - abs(sum(p_bar_prev))*dx) / abs((sum(p_bar_prev))*dx);
    sprintf("iter: %d conv: %f",iter, convergence)
   
    if convergence < power(10,-4)
        flag = 1;
    end
    
    
  drawnow
  plot(p_bar);
  title([ 'PRESSURE DISTRIBUTION' '    ' 'for' '    ''Attitude Ratio:' '    ' num2str(n)])
  ylabel('Non-dimentional pressure');
  xlabel('Non-dimentional width');
       
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculating Load Carrying Capacity
    
    %W = dot(p_bar,dx_mat); 
    
    %Using Trapeziodal Rule for Approximating Integration
    
    %temp = (p_bar(1) + p_bar(nxnodes)) * (dx/2);  % First and last terms : f(a) * dx/2 + f(b) * dx/2 
    %Load_capacity = temp + sum(p_bar(2:nxnodes-1))*dx; % center terms
    Load_capacity = sum(p_bar)* dx;
    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculating the Shear Stress Distribution at nx nodes

for i = 1 : nodes
    
    % Forward Difference Method to Approximate the Derivative at Node 1
    if (i == 1)
        tau_bar(i) = (((3 * h_bar(i)) * (p_bar(i+1) - p_bar(i)) / dx) + (1 / h_bar(i)));
    end
    
    % Central Difference Method to Approximate the central Derivatives 
    if (i > 1 && i < nodes)
        tau_bar(i) = (((3 * h_bar(i)) * (p_bar(i+1) - p_bar(i-1)) / (2*dx)) + (1 / h_bar(i)));
    end
    
    % Backward Difference Method to Approximate the Derivative at node nx
    if (i == nodes)
        tau_bar(i) = (((3 * h_bar(i)) * (p_bar(i) - p_bar(i-1)) / dx) + (1 / h_bar(i)));
    end
       
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculating the Friction Force 

% Using Trapezoidal Rule to Approximate the Integration

temp = (tau_bar(1) + tau_bar(nodes)) * (dx/2);  % First and last terms : f(a) * dx/2 + f(b) * dx/2 
Friction_force = temp + sum(tau_bar(2:nodes-1))*dx; % center terms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculating the Friction Coefficient 

myu = Friction_force / (6 * Load_capacity); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t_time=clock;
disp(['================================ ',date,' ================================'])
disp(['============= Steady State Analysis of Hydrodynamic Slider Bearings ============'])
disp(['================================= Time ',num2str(t_time(4)),':',num2str(t_time(5)),' ================================='])    
sprintf("Load Carrying Capacity (Non-Dimensionalized Value) : %f", Load_capacity)
disp('*****************************************************************************')
sprintf("Friction Force Acting  (Non-Dimensionalized Value) : %f", Friction_force)
disp('*****************************************************************************')
sprintf("Coefficient of Friction (Non-Dimensionalized Value) : %f", myu)
disp('*****************************************************************************')



end
