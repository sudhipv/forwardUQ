% ExtractGmsh
file = fullfile('square.msh');

% 4th line gives number of nodes in the mesh

nodenumber = dlmread(file,' ',[4 0, 4 0]);

% Lines from 0-4 in .msh file gives descriptions

nodes = dlmread(file,' ',[5 0 nodenumber+4 3]);

% 3 Lines after node number are descriptions
% 3rd line after nodes provides number of elements

elenumber = dlmread(file,' ',[nodenumber+7 0, nodenumber+7 0]);

elements = dlmread(file,' ',[nodenumber+8 0 elenumber+nodenumber+7 7]);

% Here line elements will have last digit assigned 0 by code which is not
% required

% extracting points/nodes edges and traingles from the extracted arrays
% from gmsh

% For 2D only 2 and 3rd columns which gives x and y co-ordinates are taken

points = nodes(1:nodenumber,2:3); % <x-co-ordinate, y co-ordinate>

% Loop to find the position of starting of triangles from element array
for i=1:elenumber
    
    if(elements(i,2) == 2)
        tristart = i;
        break
    end
end

edges = elements(1:tristart-1,5:7); % <edge , node 1, node 2 > 

triangles = elements(tristart:elenumber, 6:8); % <node1, node 2, node3>

triangles(:,4) = 1;

p = points';
e =edges';
t = triangles';

%Save Mesh data
dlmwrite('points.txt', p',' ')
dlmwrite('edges.txt', e',' ')
dlmwrite('triangles.txt', t',' ')






