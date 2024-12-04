function varargout = mis_read_fortran(fid,n,precision)
% Function to read an unformatted fortran data file.
% ---
% Set up for Tim Houghton's version of Mises, compiled using gfortran.

% Read data
head=fread(fid,2,'int16'); % fortran u/f header
data=fread(fid,n,precision);
foot=fread(fid,2,'int16'); % fortran u/f footer

% Check header and footer
if any(head ~= foot); warning('Fortran unformatted header does not equal footer.'); end

% Manage output arguments
if nargout == 1
    varargout{1} = data;
else
    for ii = 1:nargout; varargout{ii}=data(ii); end
end

return