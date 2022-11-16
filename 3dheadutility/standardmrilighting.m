function myLight = standardMRIlighting(option,handles)
% function myLight = standardMRIlighting(option,handles)
% option = 'full' - all lighting adjustments
%          'reflectance' - ambient strength and reflectance only
% handles = [isosurfaceHandle isocapsHandle]            

if strcmp(option,'full')
    
    view(135,30) 

    myLight = lightangle(45,30); 
    %set(gcf,'Renderer','zbuffer'); 
    lighting phong

end

set(handles(2),'AmbientStrength',.6)
set(handles(1),'SpecularColorReflectance',0,'SpecularStrength',.2,'SpecularExponent',200)

if strcmp(option,'full')
    % this is the key to avoiding dark surface head.
    [az,el] = view;
    myLight = lightangle(az,el);
else myLight = [];
end

axis image

drawnow

return