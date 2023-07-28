filew = fopen("quadrati-distorti-0.5-10x10.txt","w");
fprintf(filew,"# domain type\n");
fprintf(filew,"RectangularDomain\n");
fprintf(filew,"# nodal coordinates: number of nodes followed by the coordinates\n");
fprintf(filew,"%d\n",mesh.NV);

arrx={mesh.vertex.x};
arry={mesh.vertex.y};
xx=[];
yy=[];
for x=arrx
    xx=[xx; x{1}];
end
for y=arry
    yy=[yy; y{1}];
end

%Inserisci i nodi
for i=1:mesh.NV

    fprintf(filew, "%.16f %.16f\n",xx(i),yy(i));

end

fprintf(filew,"# element connectivity: number of elements followed by the connectivity (each line: nodes_per_element(nel) node1 node 2 ... node_nel)\n");
fprintf(filew,"%d\n",mesh.NP);

conn={mesh.polygon.vertices};
connv=[];

for x = conn
    fprintf(filew, "%d %d ", numel(cell2mat(x)), cell2mat(x));
    fprintf(filew,"\n");
end

fprintf(filew,"# indices of nodes located on the bottom boundary\n");

for i=1:mesh.NV

    if (yy(i) == 0.0)
        fprintf(filew,"%d ",i);
    end

end

fprintf(filew,"\n# indices of nodes located on the top boundary\n");

for i=1:mesh.NV

    if (yy(i) == 1.0)
        fprintf(filew,"%d ",i);
    end

end

fprintf(filew,"\n# indices of nodes located on the left boundary\n");

for i=1:mesh.NV

    if (xx(i) == 0.0)
        fprintf(filew,"%d ",i);
    end

end

fprintf(filew,"\n# indices of nodes located on the right boundary\n");

for i=1:mesh.NV

    if (xx(i) == 1.0)
        fprintf(filew,"%d ",i);
    end

end

fprintf(filew,"\n# xmin, xmax, ymin, ymax of the bounding box\n");
fprintf(filew,"0 1 0 1");
fprintf(filew,"\n");
fclose(filew);

%close all
%clear all