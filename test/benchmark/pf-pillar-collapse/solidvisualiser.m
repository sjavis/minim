filename_s = 'solid_output.txt';
filename_l = 'liquid_output.txt';
nx = 120;
ny = 100;
nz = 120;

%%  
solid_data = importdata(filename_s);
solid_data = reshape(solid_data,[nx,ny,nz]);

liquid_data = importdata(filename_l);
liquid_data = reshape(liquid_data,[nx,ny,nz]);

%%
[x,y,z] = meshgrid(1:1:ny,1:1:nx,1:1:nz);

isosurface(x,y,z,solid_data,'0.5')
isosurface(x,y,z,liquid_data,'0.5')
view(70,-40)
camlight
lighting gourand
colormap autumn   
