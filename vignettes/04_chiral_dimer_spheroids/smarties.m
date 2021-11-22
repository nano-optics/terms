clearvars

wavelength = 400:2:800;
epsilon=epsAu(wavelength);
medium=1.33;

% stParams.a=20; stParams.c=50;
stParams.a=10; stParams.c=20;
stParams.N=4; stParams.nNbTheta=40;

stOptions.bGetR = false;
stOptions.Delta = 0;
stOptions.NB = 0;
stOptions.bGetSymmetricT = false;

out = ['tmat_Au',num2glimpse(stParams.a),'x',num2glimpse(stParams.c), ...
       '_Nmax',num2glimpse(stParams.N),'.tmat'];
format = '%3d %3d %3d %3d %3d %3d %.15e %.15e\n';

fileID = fopen(out, 'w');
fprintf(fileID, '# s sp n np m mp Tr Ti | a= %g c= %g\n', stParams.a, stParams.c);
for ii=1:length(wavelength)
    stParams.k1=medium*2*pi/wavelength(ii); 
    stParams.s=sqrt(epsilon(ii)) / medium;
    [~, stT] = slvForT(stParams,stOptions);
    T = exportTmatrix( stT, true, [], format );
    fprintf(fileID, '# lambda= %g nelements= %d epsIn= %f%+fj\n', ...
        wavelength(ii),size(T, 1),real(epsilon(ii)), imag(epsilon(ii)));
    fprintf(fileID, format, T.');
end
fclose(fileID);
