function panelHandle = drawMeshList(mesh_list, options)
%DRAWMESHLIST: draw all meshes in mesh_list

GroupSize = length(mesh_list);
DisplayLayout = getoptions(options, 'DisplayLayout', [1,GroupSize]);
linkCamera = getoptions(options, 'linkCamera', 'off');
showNames = getoptions(options, 'showNames', 'on');
DisplayOrient = getoptions(options, 'DisplayOrient', 'Horizontal');
Shading = getoptions(options, 'Shading', 'Smooth');

switch DisplayOrient
    case 'Vertical'
        DisplayOrder = reshape(1:DisplayLayout(1)*DisplayLayout(2), DisplayLayout(2), DisplayLayout(1));
        DisplayOrder = DisplayOrder';
        DisplayOrder = DisplayOrder(:);
    case 'Horizontal'
        DisplayOrder = 1:DisplayLayout(1)*DisplayLayout(2);
end

figurePosition = [10, 10, 800, 400];

panelHandle = figure('Unit', 'pixel', 'Position', figurePosition);
set(gcf, 'ToolBar', 'none');
h = zeros(size(mesh_list));

for i=1:GroupSize
    color_data = repmat([0.9, 0.9, 0.8], mesh_list{i}.nV, 1);
    h(i) = subplot(DisplayLayout(1), DisplayLayout(2), DisplayOrder(i));
    if strcmpi(Shading, 'Smooth')
        mesh_list{i}.draw(struct('FaceColor', 'interp',...
            'FaceVertexCData', color_data, 'CDataMapping','scaled',...
            'EdgeColor', 'none', 'FaceAlpha', 1,...
            'AmbientStrength',0.3,'SpecularStrength',0.0));
        hold on;
        colormap jet(256);
        material shiny
        set(gca, 'CameraUpVector',  [-0.2666,-0.9468, 0.1804]);
        set(gca, 'CameraPosition',  [ 5.4366,-1.4825, 0.2070]);
        set(gca, 'CameraTarget',    [ 0.0059, 0.0039,-0.0166]);
        set(gca, 'CameraViewAngle', 10.0659);
        camlight('headlight');
        camlight(180,0);
    else
        mesh_list{i}.draw();
    end
    
    if (isfield(mesh_list{i}.Aux, 'Aname') && strcmpi(showNames, 'on'))
        title(mesh_list{i}.Aux.name);
    end
end

if strcmpi(linkCamera, 'on')
    Link = linkprop(h, {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
    setappdata(gcf, 'StoreTheLink', Link);
end

% if (exist('camUpVector', 'var'))
%     set(gca, 'CameraUpVector', camUpVector);
%     set(gca, 'CameraPosition', camPosition);
%     set(gca, 'CameraTarget', camTarget);
%     set(gca, 'CameraViewAngle', camViewAngle);
%     set(gcf, 'Tag', 'BundleFunc');
% end

end

% function draw_landmarks(V,IndsOnSource,Type)
% 
% if nargin<3
%     Type = 'Inds';
% end
% 
% if strcmpi(Type,'Coords')
%     Landmarks = IndsOnSource;
% else
%     if (size(V,1)==3)
%         V = V';
%     end
%     Landmarks = V(IndsOnSource,:)';
% end
% 
% NumLandmarks = max(size(Landmarks));
%     
% colmap =  [1,0,0;0,1,0;0,0,1;1,1,0;1,0,1;0,1,1];
% colmap = [colmap;colmap*0.5;colmap*0.25];
% 
% R = 0.015; % 0.025
% CORR_draw_spheres(Landmarks',R,colmap(1:NumLandmarks,:));
% 
% end
% 
