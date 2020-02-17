function Plot3DBrain(Jen,vert,face,QG)

figure
MyGray = [0.9 0.9 0.9];
cmap = colormap('jet');
cmap_out = [MyGray; cmap];
colormap(cmap_out);
h = patch('Faces', face, 'Vertices', vert,'FaceVertexCData',Jen,'FaceColor','interp');
axis equal;
axis off;
set(h,'edgecolor','none');
set(h,'AmbientStrength',1,'DiffuseStrength',1.0,'SpecularColorReflectance',0.0)
material dull;
camlight headlight; 
lighting phong;
