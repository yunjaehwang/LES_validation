clear;

x_offset = 0.55;
y_offset = 0.0;
z_offset = 0.55;

x0 = -0.05;    x1 = 0.05;   dx = 0.002;
y0 =  0.0;     y1 = 0.08;   dy = 0.002;
z0 = -0.05;    z1 = 0.05;   dz = 0.002;

XX = (x0+dx):dx:(x1-dx);
YY = (y0+dy):dy:(y1-dy);
ZZ = (z0+dz):dz:(z1-dz);

length(XX)*length(YY)*length(ZZ)

% writing file
% write profile number, total number of points
fname_vertical = 'PC_concentration.txt';
fid_vertical = fopen(fname_vertical,'w');
fprintf(fid_vertical, 'x, y, z \n');

for i_x = 1:length(XX)
    for i_y = 1:length(YY)
        for i_z = 1:length(ZZ)
            fprintf(fid_vertical, '%5.4f\t',x_offset + XX(i_x));
            fprintf(fid_vertical, '%5.4f\t',y_offset + YY(i_y));
            fprintf(fid_vertical, '%5.4f\n',z_offset + ZZ(i_z));
        end
    end
end
fclose(fid_vertical);




%%
theta = -45;
Q = [cosd(theta), sind(theta); -sind(theta), cosd(theta)];

x_offset = 0.55;
y_offset = 0.0;
z_offset = 0.55;

x0 = -0.05;    x1 = 0.05;   dx = 0.002;
y0 =  0.0;     y1 = 0.08;   dy = 0.002;
z0 = -0.05;    z1 = 0.05;   dz = 0.002;

XX = (x0+dx):dx:(x1-dx);
YY = (y0+dy):dy:(y1-dy);
ZZ = (z0+dz):dz:(z1-dz);

length(XX)*length(YY)*length(ZZ)

k=1;
for i_x = 1:length(XX)
    for i_y = 1:length(YY)
        for i_z = 1:length(ZZ)
%             fprintf(fid_vertical, '%5.3f\t',x_offset + XX(i_x));
%             fprintf(fid_vertical, '%5.3f\t',y_offset + YY(i_y));
%             fprintf(fid_vertical, '%5.3f\n',z_offset + ZZ(i_z));
            xyz(k,:) = [XX(i_x), YY(i_y), ZZ(i_z)];
            k = k+1;
        end
    end
end
xyz2 = xyz; xyz2(:,[1,3])=(Q*xyz(:,[1,3])')';
figure();
subplot(1,2,1); scatter3(xyz(:,1), xyz(:,3), xyz(:,2),'o');
subplot(1,2,2); scatter3(xyz2(:,1), xyz2(:,3), xyz2(:,2),'o');

% writing file
% write profile number, total number of points
fname_vertical = 'PC_concentration_45.txt';
fid_vertical = fopen(fname_vertical,'w');
fprintf(fid_vertical, 'x, y, z \n');
for k=1:length(xyz2)
    fprintf(fid_vertical, '%5.4f\t',x_offset + xyz2(k,1));
    fprintf(fid_vertical, '%5.4f\t',y_offset + xyz2(k,2));
    fprintf(fid_vertical, '%5.4f\n',z_offset + xyz2(k,3));
end
fclose(fid_vertical);




