filename_s = 'solid_output.txt';
filename_l = 'liquid_output.txt';
nx = 100;
ny = 100;
nz = 100;

%%  
solid_data = importdata(filename_s);
solid_data = reshape(solid_data,[nx,ny,nz]);

liquid_data = importdata(filename_l);
liquid_data = reshape(liquid_data,[nx,ny,nz]);

%%
[x,y,z] = meshgrid(1:1:ny,1:1:nx,1:1:nz);

isosurface(x,y,z,solid_data,'0.5')
colormap("autumn")
%xlim([0,ny])
ylim([0,nx])
zlim([0,nz])
%isosurface(x,y,z,liquid_data,'0.5')
%colormap("cool")
view(70,-40)



