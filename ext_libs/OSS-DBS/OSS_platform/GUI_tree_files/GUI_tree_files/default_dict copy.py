"""@author: trieu,butenko"""
d = {
    'Segm_MRI_processed': 0,
    'DTI_processed': 0,
    'Init_neuron_model_ready': 0,
    'Init_mesh_ready': 0,
    'Adjusted_neuron_model_ready': 0,
    'CSF_mesh_ready': 0,
    'Adapted_mesh_ready': 0,
    'Signal_generated': 0,
    'Parallel_comp_ready': 0,
    'Parallel_comp_interrupted': 0,
    'IFFT_ready': 0,
    'MRI_data_name': 'segmask.nii',
    'MRI_in_m': 0,
    'DTI_data_name': 'Johnson_WS_NormMapping.nii',
    'DTI_in_m': 0,
    'CSF_index': 3.0,
    'WM_index': 2.0,
    'GM_index': 1.0,
    'default_material': 3,
    'Electrode_type': 'AA_rodent_monopolar',
    'Brain_shape_name': '',
    'x_length': 0.0,
    'y_length': 0.0,
    'z_length': 0.0,
    'Aprox_geometry_center': [2.476187675483748, -3.0393966740032874, -0.8115758910846853],
    'Approximating_Dimensions': [80.0, 80.0, 80.0],
    'Implantation_coordinate_X': 2.476187675483748,
    'Implantation_coordinate_Y': -3.0393966740032874,
    'Implantation_coordinate_Z': -0.8115758910846853,
    'Second_coordinate_X': 2.4526094515956562,
    'Second_coordinate_Y': -2.99219413024874,
    'Second_coordinate_Z': -0.195079703163476,
    'Rotation_Z': -0.0,
    'encap_thickness': 0.1,
    'encap_tissue_type': 2,
    'encap_scaling_cond': 1.0,
    'encap_scaling_perm': 1.0,
    'pattern_model_name': '',
    'Axon_Model_Type': 'McIntyre2002',
    'diam_fib': [2.0],
    'n_Ranvier': [5],
    'v_init': -80.0,
    'Neuron_model_array_prepared': 1,
    'Name_prepared_neuron_array': 'Allocated_axons.h5',
    'Global_rot': 1,
    'x_seed': 10.037875525508197,
    'y_seed': -18.294343528964472,
    'z_seed': -6.588918568320674,
    'x_steps': 4,
    'y_steps': 0,
    'z_steps': 4,
    'x_step': 0.5,
    'y_step': 0.5,
    'z_step': 0.5,
    'alpha_array_glob': [0],
    'beta_array_glob': [0],
    'gamma_array_glob': [0],
    'X_coord_old': 0,
    'Y_coord_old': 0,
    'Z_coord_old': 0,
    'YZ_angles': [0],
    'ZX_angles': [0],
    'XY_angles': [0],
    'EQS_core': 'QS',
    'Skip_mesh_refinement': 1,
    'refinement_frequency': '',
    'num_ref_freqs': -1,
    'rel_div_CSF': -1.0,
    'Adaptive_frac_div': 1.0,
    'Min_Scaling': 1.0,
    'CSF_ref_reg': 0.0,
    'rel_div': 10.0,
    'rel_div_current': 1.0,
    'el_order': 3,
    'number_of_processors': 0,
    'FEniCS_MPI': 0,
    'current_control': 1,
    'Pulse_amp': [-0.002, 0.0],  # can be negative (i.e. cathode in VC)
    'Solver_Type': 'GMRES',
    'freq': 130.0,
    'T': 60.0,
    't_step': 1.0,
    'phi': 0.0,
    'Signal_type': 'Rectangle',
    'Ampl_scale': 1.0,
    'CPE_activ': 0,
    'beta': 0.0,
    'K_A': 0.0,
    'beta_ground': 0.0,
    'K_A_ground': 0.0,
    'VTA_approx': 0,  # !!! Here it refers to |E|-field based VTA approximation !!!
    't_step_end': 1200,
    'VTA_from_divE': 0,
    'VTA_from_NEURON': 0,
    'VTA_from_E': 1,
    'Activation_threshold_VTA': 0.343,
    'spectrum_trunc_method': 'Octave Band Method',
    'trunc_param': 130.0,
    'Truncate_the_obtained_full_solution': 0,
    'Show_paraview_screenshots': 0,
    'external_grounding': True,
    'Stim_side': 0,
    'patient_folder' : 'John_Doe',   
}