function id_shift = ifftshiftId1D(Nfft,id)
% Compute the new indices of the elements of a vector after a 1D ifftshift.
%-------------------------------------------------------------------------%
% Input:
% > Nfft : vector size
% > id   : vector containing the indices of interest
%
% Output:
% < id_shift : indices corresponding to the elements of id after a 1D
%              ifftshift
%
%-------------------------------------------------------------------------%
%% 
% [26/10/2017], tests OK
% [15/03/2018], P.-A. Thouvenin.
%-------------------------------------------------------------------------%

id1 = (id <= round(Nfft/2));
id_shift = zeros(size(id));
id_shift(id1) = id(id1) + floor(Nfft/2); % Nfft/2 if even, (Nfft+1)/2 otherwise 
id_shift(~id1) = id(~id1) - round(Nfft/2);

end

