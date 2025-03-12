function rhof_veff_ratio_struct = get_all_rhof_veff_ratios(training_data_general_params, mask_angle_deg)
    ephs = extract_all_rinex_eph(training_data_general_params.date_time);

    mask_angle_rad = mask_angle_deg * pi/180;
    sats_in_los_struct = get_sats_in_los(ephs,training_data_general_params, mask_angle_rad);
    
    rhof_veff_ratio_struct = struct();
end