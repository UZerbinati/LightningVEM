mesh_file = fopen("C:\Users\Manuel.LAPTOP-LBO3KFQI\Desktop\VEM _v102\mesh_files\square_4096.txt");

sscanf(fgets(mesh_file),'%d');
sscanf(fgets(mesh_file),'%d');
sscanf(fgets(mesh_file),'%d');
nodes_number = sscanf(fgets(mesh_file),'%d');

filew = fopen("trisquare40955556.txt","w");
fprintf(filew,"# domain type\n");
fprintf(filew,"RectangularDomain\n");
fprintf(filew,"# nodal coordinates: number of nodes followed by the coordinates\n");
fprintf(filew,"%d\n",nodes_number(1));


coords = zeros(nodes_number(1),2);

for i=1:nodes_number

    coords(i,:) = sscanf(fgets(mesh_file),'%f');
    fprintf(filew, "%.16f %.16f\n",coords(i,:));

end

sscanf(fgets(mesh_file),'%d');
fprintf(filew,"# element connectivity: number of elements followed by the connectivity (each line: nodes_per_element(nel) node1 node 2 ... node_nel)\n");
n_elem = sscanf(fgets(mesh_file),'%d');
fprintf(filew,"%d\n",2*n_elem);

for i=1:n_elem

    vec = sscanf(fgets(mesh_file),'%d');
    fprintf(filew, "3 %d %d %d \n", vec(2:4)');
    fprintf(filew, "3 %d %d %d \n", vec(4:5)', vec(2)');

end


sscanf(fgets(mesh_file),'%d');
fprintf(filew,"# indices of nodes located on the bottom boundary\n");
vec = sscanf(fgets(mesh_file),'%d');
fprintf(filew, "%d ", vec');

sscanf(fgets(mesh_file),'%d');
fprintf(filew,"\n# indices of nodes located on the top boundary\n");
vec = sscanf(fgets(mesh_file),'%d');
fprintf(filew, "%d ", vec');

sscanf(fgets(mesh_file),'%d');
fprintf(filew,"\n# indices of nodes located on the left boundary\n");
vec = sscanf(fgets(mesh_file),'%d');
fprintf(filew, "%d ", vec');


sscanf(fgets(mesh_file),'%d');
fprintf(filew,"\n# indices of nodes located on the right boundary\n");
vec = sscanf(fgets(mesh_file),'%d');
fprintf(filew, "%d ", vec');

fprintf(filew,"\n# xmin, xmax, ymin, ymax of the bounding box\n");
fprintf(filew,"0 1 0 1");
fprintf(filew,"\n");
fclose(filew);

close all
clear all