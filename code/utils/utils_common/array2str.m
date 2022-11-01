function S = array2str(A,Delimiter)
%% Function Description:
% This function creates a delimiter-separated list (string) of the rows in the array A. 
%
% AUTHOR: Sugato Ray | Created on: 25-APR-2017 | ray.sugato[at]gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           PLEASE ACKNOWLEDGE THE AUTHOR IF YOU USE THIS CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
%   A = array 
%   Delimiter = ',' is default. => This creates a CSV String.
%                   The user could choose any characher as the delimiter.
%
% OUTPUT:
%   S = delimiter separated list (string) of the rows of A. 
%
% EXAMPLE:
%   S = array2str(A,Delimiter);
%
%--------------------------------------------------------------------------    
    if nargin<2
        Delimiter = ',' ; % Comma is default delimiter
    end    
    S = strjoin(arrayfun(@(x) num2str(x),A,'UniformOutput',false),Delimiter);

end