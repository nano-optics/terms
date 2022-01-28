%% Testing CD of gold dimer
clearvars
close
wavelength = 450:2:850;
% wavelength = 450:100:850;

material.wavelength = wavelength;
material.epsilon = epsilon_Au(wavelength);
material.epsilon = feval(@epsilon_Au, wavelength);
%material.fun_polarisability = @alpha_majic; % more precise than kuwata
  
material.medium=1.33;

a=20;b=a;c=50;
d=100;

cl =  cluster_fingers(d, a, c, pi/4);

xsec = spectrum_oa(cl, material, 'gl',300);
    
clref =  cluster_single(a, b, c , 0, 0, 0);
ref = spectrum_oa(clref, material,'gl',300);

% 
% 
% xsec = simulate_fingers(a, c, d, phi, wavelength, epsilon, medium, N_inc, N_sca, majic);
%                               
% cl =  cluster_dimer(d, a, b, c, pi/4, 0, 0);
% clref =  cluster_single(a, b, c , 0, 0, 0);
% 
% xsec = spectrum_oa(cl, material,'gl',300);
% Incidence = [0,0;0,pi/2;0,0];
% xsec = spectrum_dispersion(cl, material,Incidence);
% % xsec = spectrum_oa(cl, material,'gl',300);
% xsec2 = spectrum_oa(clref, material,'gl',300);


set(groot,'defaultAxesColorOrder', visual_colorscale(2), 'defaultLineLineWidth', 2)
plot(wavelength, xsec(:,1:2),wavelength, ref(:,1:2),':')

save('fingers.mat', 'wavelength', 'xsec', 'ref')
%  plot(wavelength, xsec(:,5),wavelength, xsec2(:,5),':')
% plot(wavelength, xsec.cext,wavelength, xsec.csca,':')
% plot(wavelength, xsec.cext-xsec.cabs,wavelength, xsec.csca,':')
