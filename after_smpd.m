 [resultst] = getting_var_smpd(results);
 %[MTLt] = getting_var_smpd(MTLent);
 [parameterst] = getting_var_smpd(parameters);
  [ocupationnumbert] = getting_var_smpd(ocupationnumber);
   [MTarryocupationt] = getting_var_smpd(MTarryocupation);
   density_estimation;
   Getting_results_from_spreading;
   saving_to_txtfile;
    getting_parameters_final_simul;
    density(density==0)=NaN;
    nanmean(density,2)./(2.*pi.*1000000*(vectorparameters(:,end)))