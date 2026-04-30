path(pathdef); % clear previous path changes
addpath(genpath('~/Documents/nano-optics/smarties/'));
addpath(genpath('~/Documents/nano-optics/easyh5/'));
clearvars;

%% example

% prolate Au spheroid in water
% semi-axes a=b=20nm, c=40nm

% wavelengths and epsilon loaded from file (interpolated JJ data)
load('epsAuJC.mat')
epsAuJC
wavelength = epsAuJC.wavelength(1:1:end)';
Nl = length(wavelength);
epsilon = epsAuJC.epsilon(1:1:end).';

medium=1.33;

% constant simulation parameters
stParams.a=20;
stParams.c=40;

% internal options
stOptions.bGetR = false;
stOptions.Delta = 0;
stOptions.NB = 0; % NB will be estimated automatically
stOptions.bGetSymmetricT = false;
stOptions.bOutput = false; % verbosity


globalN = 3;
globalnNbTheta = 300;

% allocate 3D array for all results
qmax = 2*(globalN*(globalN + 1) + globalN );
tmatrix = zeros(Nl, qmax, qmax);

out = ['smarties.txt'];
fileID = fopen(out, 'w');
format = '%3d %3d %3d %3d %3d %3d %.15e %.15e\n';
fprintf(fileID, '# s sp n np m mp Tr Ti | a= %g c= %g\n', stParams.a, stParams.c);

for i=1:Nl
    stParams.k1=medium*2*pi/wavelength(i);
    stParams.s=sqrt(epsilon(i)) / medium;

    stParams.N=globalN; stParams.nNbTheta=globalnNbTheta;

    [stCoa, stT] = slvForT(stParams,stOptions);

    [T, q, qp] = exportTmatrix( stT, true, [], format);
    % plain text version for comparison
    fprintf(fileID, '# lambda= %g nelements= %d epsIn= %f%+fj\n', ...
        wavelength(i),size(T, 1),real(epsilon(i)), imag(epsilon(i)));
    fprintf(fileID, format, T.');

    % convert these indices to (u,up) for treams convention
    [u] = treams_indexing(q, qmax);
    [up] = treams_indexing(qp, qmax);
    % and find linear indices in 3D array to insert these elements
    ind = sub2ind(size(tmatrix),q*0+i, u, up);
    tmatrix(ind) =  T(:,7) + 1i*T(:,8);

end
fclose(fileID);

analytical_zeros = true(qmax,qmax);
analytical_zeros(sub2ind(size(analytical_zeros), u, up)) = false;


%% data to export

epsilon = struct(...
    'embedding', medium^2,...
    'particle', epsilon, ...
    'embedding_name', 'H2O, Water', ...
    'embedding_keywords', 'non-dispersive',...
    'embedding_reference', 'constant', ...
    'material_name', 'Au, Gold', ...
    'material_reference', 'Au from Raschke et al 10.1103/PhysRevB.86.235147',...
    'material_keywords', 'dispersive, plasmonic');

geometry = struct('radiusxy', stParams.a, 'radiusz', stParams.c);


computation = struct('description', 'Computation using SMARTIES, a numerically robust EBCM implementation for spheroids', ...
    'accuracy', 0.0, ... % not determined in this example
    'Lmax', globalN, ...
    'Ntheta', globalnNbTheta, ...
    'analytical_zeros', analytical_zeros);

comments = struct('name', 'Au prolate spheroid in water',...
    'description', 'Computation using SMARTIES, a numerically robust EBCM implementation for spheroids',...
    'keywords', 'gold, spheroid, ebcm', ...
    'script', [mfilename '.m']);


[f] = tmatrix_hdf5('smarties.tmat.h5', tmatrix, wavelength, epsilon, geometry, computation, comments)
