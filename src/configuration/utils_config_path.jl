#=
function config_path(parameter::run_parameter, t::Int)
  FRAC        = parameter.FRAC
  COMP        = parameter.COMP
  STEP        = parameter.STEP
  ZLIM        = parameter.ZLIM
  CHI         = parameter.CHI
  N_vpp       = parameter.N_vpp
  interpolate = parameter.interpolate
  sigma       = parameter.sigma
  #--------------------------------------------------------------------
  dir_IO   = "IO"     
  dir_meta = "Meta"   
  dir_plot = "Plots" 

  #prefix_path   = "CS_"
  prefix_path   = "CS_NEW_"
  prefix_meta   = "meta_"                 #File where metadata will be written. Don't use preexisting names.
  prefix_in_out = "data_in_out_"
  prefix_skymap_2000_1000 = "skymap_gamma_ray_2000_1000_"
  
  #prefix_path   = "CS_DIRECT_"
  #prefix_meta   = "meta_DIRECT_"
  #prefix_in_out = "data_in_out_DIRECT_"
  #prefix_skymap_2000_1000 = "skymap_gamma_ray_2000_1000_DIRECT_"

  part1 = "09-16Rlc_z"
  part2 = "_OHM_kerr01_pulsar005_deg"
  #part3 = "_r320_sig16_res81_t$t"
  #part3_skymap = "_r320_sig16_res81"

  part3 = "_l6_r320_sig16_res81_t$t"
  part3_skymap = "_l6_r320_sig16_res81"

  frac_modifier = replace(string(FRAC), "." => "")
  comp_modifier = replace(string(COMP), "." => "")
  step_modifier = replace(string(STEP), "." => "")
  step_modifier = string("_int_", step_modifier, "step")
  sigma_modifier = replace(string(sigma), "." => "")

  #filename = string(part1, string(ZLIM), part2, string(CHI), "_comp", comp_modifier, part3)  
  #filename_skymap = string(part1, string(ZLIM), part2, string(CHI), "_comp", comp_modifier, part3_skymap)
  
  filename = string(part1, string(ZLIM), part2, string(CHI), part3)
  filename_skymap = string(part1, string(ZLIM), part2, string(CHI), part3_skymap)

  folder_name = string("/C", comp_modifier, "_deg", string(CHI), "/")
  #folder_name = string("/C", comp_modifier, "_deg", string(CHI), "_direct/")
  
  mkpath(dir_meta * folder_name)
  mkpath(dir_plot * folder_name)
  mkpath(dir_IO * folder_name)

  if interpolate 
    suffix_ONION             = ".csv"
    suffix_data_interpolated = string(step_modifier, ".csv")
    suffix_meta              = string(step_modifier, "_Nvpp$(N_vpp)_sigma$(sigma_modifier).yml")
    suffix_IO                = string(step_modifier, "_Nvpp$(N_vpp)_sigma$(sigma_modifier).h5")
    suffix_skymap_int        = string(step_modifier, "_Nvpp$(N_vpp)_Rcyl", frac_modifier, "_sigma$(sigma_modifier).csv")
    suffix_plot_int          = string(step_modifier, "_Nvpp$(N_vpp)_Rcyl", frac_modifier, "_sigma$(sigma_modifier).png")

    path_data_ONION          = string(dir_IO, folder_name, prefix_path, filename, suffix_ONION)
    path_data_ONION_filt_int = string(dir_IO, folder_name, prefix_path, filename, suffix_data_interpolated)
    
    path_data_save_IO = string(dir_IO, folder_name, prefix_in_out, filename, suffix_IO)  
    path_meta_file    = string(dir_meta, folder_name, prefix_meta, filename, suffix_meta)
    
    path_save_plot_skymap_2000_1000 = string(dir_plot, folder_name, prefix_skymap_2000_1000, filename, suffix_plot_int)
    path_data_save_skymap_2000_1000 = string(dir_IO, folder_name, prefix_skymap_2000_1000, filename, suffix_skymap_int)
    path_data_save_total_skymap_2000_1000 = string(dir_IO, folder_name, prefix_skymap_2000_1000, filename_skymap, suffix_skymap_int)
  else 
    suffix_ONION          = ".csv"
    suffix_ONION_filtered = "_filtered_data_ONION.csv"
    suffix_meta           = "_Nvpp$(N_vpp)_sigma$(sigma_modifier).yml" 
    suffix_IO             = "_Nvpp$(N_vpp)_sigma$(sigma_modifier).h5"
    suffix_skymap         = string("_Nvpp$(N_vpp)_Rcyl", frac_modifier, "_sigma$(sigma_modifier).csv")
    suffix_plot           = string("_Nvpp$(N_vpp)_Rcyl", frac_modifier, "_sigma$(sigma_modifier).png")

    path_data_ONION          = string(dir_IO, folder_name, prefix_path, filename, suffix_ONION)
    path_data_ONION_filt_int = string(dir_IO, folder_name, prefix_path, filename, suffix_ONION_filtered)
    
    path_data_save_IO = string(dir_IO, folder_name, prefix_in_out, filename, suffix_IO)
    path_meta_file    = string(dir_meta, folder_name, prefix_meta, filename, suffix_meta)
    
    path_save_plot_skymap_2000_1000 = string(dir_plot, folder_name, prefix_skymap_2000_1000, filename, suffix_plot)
    path_data_save_skymap_2000_1000 = string(dir_IO, folder_name, prefix_skymap_2000_1000, filename, suffix_skymap)
    path_data_save_total_skymap_2000_1000 = string(dir_IO, folder_name, prefix_skymap_2000_1000, filename_skymap, suffix_skymap)
  end
  return (path_data_ONION, path_data_ONION_filt_int, path_data_save_IO, path_meta_file, path_data_save_skymap_2000_1000, path_save_plot_skymap_2000_1000, path_data_save_total_skymap_2000_1000)  
end

=#