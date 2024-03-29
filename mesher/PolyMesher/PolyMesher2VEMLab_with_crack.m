%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                    xVEMLab
%                          Source code  : not released
%              (See Copyright and License notice in "license.txt")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% FUNCTION:              PolyMesher2VEMLab_with_crack
%
% Created by : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
% Updated by :
% Updated by :
%
%-------------------------------------------------------------------------------
% Purpose
% =======
% Read a PolyMesher [1] mesh and write it into a VEMLab mesh format.
%
% Usage
% =====
% PolyMesher2VEMLab_with_crack(Node,Element,NElem,BoundaryNodes,MeshFile,...
%                              DomainType,BdBox,EnrichedNodes)
%
% Input
% =====
% Node    : PolyMesher array containing the nodal coordinates
% Element : PolyMesher cell array containing the element connectivity 
% NElem   : number of polygonal elements
% BoundaryNodes : Structure containing the boundary nodes
% MeshFile : Location and name of the file for writing the mesh
% DomainType: type of domain
% BdBox: bounding box coordinates of the domain
% EnrichedNodes: [...; node_number enrichment_type; ...]; (enrichment_type: 1 (Heaviside) or 2 (crack-tip))
%
% Output
% ======
%
%-------------------------------------------------------------------------------
% References 
% ==========
% [1] C Talischi, GH Paulino, A Pereira, IFM Menezes, 
%     "PolyMesher: A general-purpose mesh generator for polygonal elements 
%     written in Matlab", Struct Multidisc Optim, 2012,
%     DOI 10.1007/s00158-011-0706-z                                      
%
%-------------------------------------------------------------------------------
% Function's updates history
% ==========================
% Dec 31, 2019: first release (by A. Ortiz-Bernardin)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PolyMesher2VEMLab_with_crack(Node,Element,NElem,BoundaryNodes,MeshFile,...
                                      DomainType,BdBox,EnrichedNodes)
  % for now, only rectangular domain is available
    PolyMesher2VEMLab_with_crack_rectangular_domain(Node,Element,NElem,BoundaryNodes,...
                                                    MeshFile,DomainType,EnrichedNodes);

end

function PolyMesher2VEMLab_with_crack_rectangular_domain(Node,Element,NElem,...
                                                         BoundaryNodes,MeshFile,...
                                                         DomainType,EnrichedNodes)
  fprintf('Printing mesh to a VEMLab mesh format...\n'); 
  fid = fopen(MeshFile,'w');
  % print domain type
  fprintf(fid,'# domain type\n');  
  fprintf(fid,'%s\n',DomainType);  
  % print nodal coordinates
  fprintf(fid,'# nodal coordinates: number of nodes followed by the coordinates\n');
  nnode = size(Node,1);
  fprintf(fid,'%d\n',nnode);                                    
  for node_i = 1:nnode
    fprintf(fid,'%.16f %.16f\n', Node(node_i,1), Node(node_i,2));  
  end
  % print element connectivity
  fprintf(fid,'# element connectivity: number of elements followed by the connectivity (each line: nodes_per_element(nel) node1 node 2 ... node_nel)\n');  
  fprintf(fid,'%d\n',NElem);                                 
  for el = 1:NElem
    NVertex = length(Element{el});
    fprintf(fid,'%d ', NVertex);
    for vertex = 1:(NVertex-1)
      fprintf(fid,'%d ', Element{el}(vertex));
    end
    fprintf(fid,'%d\n', Element{el}(NVertex));
  end
  % print bottom boundary  
  fprintf(fid,'# indices of nodes located on the bottom boundary\n');
  fprintf(fid,'%d ',BoundaryNodes.bottom);  
  fprintf(fid,'\n');
  % print top boundary  
  fprintf(fid,'# indices of nodes located on the top boundary\n');  
  fprintf(fid,'%d ',BoundaryNodes.top);
  fprintf(fid,'\n');  
  % print left boundary
  fprintf(fid,'# indices of nodes located on the left boundary\n');
  fprintf(fid,'%d ',BoundaryNodes.left);  
  fprintf(fid,'\n');  
  % print right boundary  
  fprintf(fid,'# indices of nodes located on the right boundary\n');  
  fprintf(fid,'%d ',BoundaryNodes.right);    
  fprintf(fid,'\n');
  % print xmin, xmax, ymin, ymax for the rectangular domain
  fprintf(fid,'# xmin, xmax, ymin, ymax of the bounding box\n');  
  xmin=Node(BoundaryNodes.blcorner,1);
  xmax=Node(BoundaryNodes.trcorner,1);
  ymin=Node(BoundaryNodes.blcorner,2);
  ymax=Node(BoundaryNodes.trcorner,2);  
  fprintf(fid,'%.16f %.16f %.16f %.16f\n', xmin, xmax, ymin, ymax);  
  % print enriched nodes
  fprintf(fid,'# number of enriched nodes followed by: node_number enrichment_type; (enrichment_type: 1 (Heaviside) or 2 (crack-tip))\n'); 
  nenriched=size(EnrichedNodes,1);
  fprintf(fid,'%d\n',nenriched);
  for node_i=1:nenriched
    fprintf(fid,'%d %d\n',EnrichedNodes(node_i,1), EnrichedNodes(node_i,2));  
  end  
  fclose(fid);
end


