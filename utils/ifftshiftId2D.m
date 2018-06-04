function id_shift = ifftshiftId2D(I,J,id)
% Compute the new indices of the elements of a matrix after a 2D ifftshift.
%-------------------------------------------------------------------------%
% Input:
% > I  : number of rows of the matrix
% > J  : number of columns
% > id : vector containing the indices of interest
%
% Output:
% < id_shift : indices corresponding to the elements of id after a 2D
%              ifftshift
%
%-------------------------------------------------------------------------%
%% 
% [26/10/2017], tests OK
% [15/03/2018], P.-A. Thouvenin.
%-------------------------------------------------------------------------%

id_2D_shift = zeros(numel(id),2); 
[idi,idj] = ind2sub([I,J], id(:));

idi1 = (idi <= round(I/2));
id_2D_shift(idi1,1) = idi(idi1) + floor(I/2); % Nfft/2 if even, (Nfft+1)/2 otherwise 
id_2D_shift(~idi1,1) = idi(~idi1) - round(I/2);

idj1 = (idj <= round(J/2));
id_2D_shift(idj1,2) = idj(idj1) + floor(J/2); % Nfft/2 if even, (Nfft+1)/2 otherwise 
id_2D_shift(~idj1,2) = idj(~idj1) - round(J/2);

% conversion to linear indices
id_shift = reshape((id_2D_shift(:,2) - 1)*I + id_2D_shift(:,1), size(id));

end

