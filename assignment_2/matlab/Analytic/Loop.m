% This script automatically loops a number of times over the analytical
% interface interpolation for a number of mesh refinement steps
% It collects the errors for the NN and RBF interpolation schemes and plots
% the results as a function of the number of structure elements in log-log
% plots

maxref = 8;
vec_N          = zeros(maxref,2);
err_vec_Df_NN  = zeros(maxref,1);
err_vec_Df_RBF = zeros(maxref,1); 
err_vec_Ps_NN  = zeros(maxref,2);
err_vec_Ps_RBF = zeros(maxref,2); 
err_vec_dW_NN  = zeros(maxref,2);
err_vec_dW_RBF = zeros(maxref,2);
for iref=1:8
    
    % Determine the number of structure ans fluid mesh points (multiply by
    % 2 for every new refinement level)
    Ns = 5*2^(iref-1)+1;
    Nf = 7*2^(iref-1)+1;
    
    Analyt;
    
    err_vec_Df_NN(iref)    = err_Df_NN;
    err_vec_Df_RBF(iref)   = err_Df_RBF;
    err_vec_Ps_NN(iref,:)  = [err_Ps_NN  err_Ps_NN_cv ];
    err_vec_Ps_RBF(iref,:) = [err_Ps_RBF err_Ps_RBF_cv];
    err_vec_dW_NN(iref,:)  = [err_dW_NN  err_dW_NN_cv ];
    err_vec_dW_RBF(iref,:) = [err_dW_RBF err_dW_RBF_cv];
    vec_N(iref,:)          = [Ns-1 Nf-1];
end

%% POST PROCESSING

% show the results for the interpolation of displacement from structure to
% flow
figure(3)
hold off;
loglog(vec_N(:,1),err_vec_Df_NN ,'r-x');
hold on;
loglog(vec_N(:,1),err_vec_Df_RBF,'g-x');
legend('NN','RBF');
xlabel('Number of structure elements')
ylabel('Error in fluid point displacement')

% show the results for the interpolation of pressures from flow to
% structure using both the consistent and conservative approaches
figure(4)
hold off;
loglog(vec_N(:,1),err_vec_Ps_NN(:,1) ,'r-x');
hold on;
loglog(vec_N(:,1),err_vec_Ps_RBF(:,1),'g-x');
loglog(vec_N(:,1),err_vec_Ps_NN(:,2) ,'r--o');
loglog(vec_N(:,1),err_vec_Ps_RBF(:,2),'g--o');
legend('NN','RBF','NN - conservative','RBF - conservative')
xlabel('Number of structure elements')
ylabel('Error in structure point pressure')

% show the results for the error in conservation of work over the interface
% for both the consistent and conservative approaches
figure(5)
hold off;
loglog(vec_N(:,1),err_vec_dW_NN(:,1) ,'r-x');
hold on;
loglog(vec_N(:,1),err_vec_dW_RBF(:,1),'g-x');
loglog(vec_N(:,1),err_vec_dW_NN(:,2) ,'r--o');
loglog(vec_N(:,1),err_vec_dW_RBF(:,2),'g--o');
legend('NN','RBF','NN - conservative','RBF - conservative')
xlabel('Number of structure elements')
ylabel('Error in work at the interface')
