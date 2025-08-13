function [parsed_argins,sim_params,log] = initialize_cpssm(cpssm_root_dir, user_inputs)
%INITIALIZE_CPSSM Summary of this function goes here
%   Detailed explanation goes here

%% instantiate simulation parameters with constant values
sim_params = get_const_sim_params();

%% handle input args
[parsed_argins, log] = parse_input_args(cpssm_root_dir, ...
    sim_params.const.all_constellations, user_inputs{:});

% set basic simulation parameters directly from the input arguments: IPP
% altitude, drift velocity and receiver parameters, etc.
sim_params = set_sim_params_from_argin(sim_params, parsed_argins);
end

