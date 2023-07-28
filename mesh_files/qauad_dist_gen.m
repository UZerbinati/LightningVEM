meshr = fopen("C:\Users\Manuel.LAPTOP-LBO3KFQI\Desktop\VEM _v102\mesh_files\square_16.txt");       %Mesh di quadrati
meshw = fopen("C:\Users\Manuel.LAPTOP-LBO3KFQI\Desktop\VEM _v102\mesh_files\square_dist_16.txt");

sscanf(fgets(meshr),'%d');
sscanf(fgets(meshr),'%d');
sscanf(fgets(meshr),'%d');
nodes_number = sscanf(fgets(meshr),'%d');

fprintf(meshw,"# domain type\n");
fprintf(meshw,"RectangularDomain\n");
fprintf(meshw,"# nodal coordinates: number of nodes followed by the coordinates\n");
fprintf(meshw,"%d\n",nodes_number(1));

coords = zeros(nodes_number(1),2);

pert = 0.5*(1/4);

for i=1:nodes_number

    pertx = (rand() * 2 - 1) * pert;
    perty = (rand() * 2 - 1) * pert;
    coords(i,:) = sscanf(fgets(mesh_file),'%f');
    fprintf(meshw, "%.16f %.16f\n",coords(i,:) + [pertx, perty]);

end

sscanf(fgets(mesh_file),'%d');
fprintf(meshw,"# element connectivity: number of elements followed by the connectivity (each line: nodes_per_element(nel) node1 node 2 ... node_nel)\n");
n_elem = sscanf(fgets(mesh_file),'%d');
fprintf(meshw,"%d\n",n_elem);

for i=1:n_elem

    vec = sscanf(fgets(mesh_file),'%d');
    fprintf(meshw, "%d \n", vec');
  

end


sscanf(fgets(mesh_file),'%d');
fprintf(meshw,"# indices of nodes located on the bottom boundary\n");
vec = sscanf(fgets(mesh_file),'%d');
fprintf(meshw, "%d ", vec');

sscanf(fgets(mesh_file),'%d');
fprintf(meshw,"\n# indices of nodes located on the top boundary\n");
vec = sscanf(fgets(mesh_file),'%d');
fprintf(meshw, "%d ", vec');

sscanf(fgets(mesh_file),'%d');
fprintf(meshw,"\n# indices of nodes located on the left boundary\n");
vec = sscanf(fgets(mesh_file),'%d');
fprintf(meshw, "%d ", vec');


sscanf(fgets(mesh_file),'%d');
fprintf(meshw,"\n# indices of nodes located on the right boundary\n");
vec = sscanf(fgets(mesh_file),'%d');
fprintf(meshw, "%d ", vec');

fprintf(meshw,"\n# xmin, xmax, ymin, ymax of the bounding box\n");
fprintf(meshw,"0 1 0 1");
fprintf(meshw,"\n");
fclose(meshw);

close all
clear all