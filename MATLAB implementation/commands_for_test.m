classdef commands_for_test
    properties (Access = public)
        % constraints variables
        lyap_const_test = true;
        ip_bound_test = false;
        state_bound_test = true;
        obs_bound_test = true;
        velocity_profile_test = true;

        % tuning variables
        ip_hard_const = true; %limit u using min and max clamping to U_min and U_max
        Ks_multiplier = 1;
        K_Ys_multiplier = 1;
        gamma_value  = 1;
        slack_cost = 10;
        R_multiplier = 100;
        Q_multiplier = 10;%only multiplier for position vector
        
        %display variables
        disp_t = true;
        disp_u = false;
        disp_goal_free = false;
        
        choosing_targets_corners = true;
    end
end