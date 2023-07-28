function [loc1, loc2, total] = cip_dofs(connect1, connect2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: cip_dofs
%
% Created by : M. Trezzi
%
%---------------------------------------------------------------------------------------------------
% Purpose
% =======
% This function determinates the DOFs that have to be used in the local CIP assembly 
%
% Input
% =====
% map      : Indexes of the vertices of the edge
% connect1 : Indexes of the vertices of the first element
% connect1 : Indexes of the vertices of the second element
%
% Output
% ======
% loc1       : Position of the DOFs of the first element in the DOFs array
% loc2       : Position of the DOFs of the second element in the DOFs array
% totalindex : The DOFs array
%
%---------------------------------------------------------------------------------------------------
% Function's updates history
% ==========================
% Mai  7, 2022: first realease (by M. Trezzi)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

[~,pos_com,~] = intersect(connect2,connect1);
%pos_com2 = find(ismember(connect2,connect1));                                                        %Position of the common indexes
n_com   = numel(pos_com);                                                                            %Number of common indexes

n_tot = numel(connect1) + numel(connect2) - n_com;                                                   %Number of unique indexes
total = zeros(n_tot,1);                                                                              %Memory allocation

loc1 = (1:numel(connect1))';                                                                         %Construction of loc1
loc2 = zeros(numel(connect2),1);                                                                     %Memory allocation

total(1:numel(connect1)) = connect1;

new             = setdiff(1:numel(connect2), pos_com);
%new2 = setdiff(1:numel(connect2), pos_com2);
[~, ~, pos_old] = intersect(connect2(pos_com),total(1:numel(connect1)));
%pos_old2 = find(ismember(total(1:numel(connect1)),connect2(pos_com2)));

%loc22 = loc2;
loc2(new) = (numel(connect1)+1:n_tot);
%loc22(new2) = (numel(connect1)+1:n_tot);
loc2(pos_com) = pos_old;
%loc22(pos_com2) = pos_old2;

%total2 = total;
total(numel(connect1)+1:end) = connect2(new);
%total2(numel(connect1)+1:end) = connect2(new2);

end