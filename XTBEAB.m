%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Splits spots based on their distance along the z axis to a surface object
% Uses minimum euclidean distance in 3d, and considers the distance along z
% Developed by DP-LAB - www.dp-lab.info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Open imaris first, make sure all the objects are visible and checked

javaaddpath ImarisLib10.jar
vImarisLib = ImarisLib;
aImarisApplicationID = round(str2double(0));
vImarisApplication = vImarisLib.GetApplication(aImarisApplicationID);
vFactory = vImarisApplication.GetFactory;

dataset = vImarisApplication.GetDataSet();
dataset_planes = dataset.GetSizeZ;

%% Choose surface
vScene = vImarisApplication.GetSurpassScene;
if isequal(vScene, [])
    msgbox 'Please create a surpass scene!';
    return
end

vSpotsList = {};
vSpotsNames = {};
vSpotsNumber = 0;

vSurfacesList = {};
vSurfacesNames = {};

vSurfacesNumber = 0;

% get all the spots and surfaces (including sub directories)
  for vChild = 1:vScene.GetNumberOfChildren
    vItem = vScene.GetChild(vChild-1);
    if isempty(vItem) || ~vItem.GetVisible
      continue
    end
    if vFactory.IsSpots(vItem)
      vSpotsNumber = vSpotsNumber + 1;
      vSpotsList{vSpotsNumber} = vFactory.ToSpots(vItem);
      n = vSpotsList{vSpotsNumber}.GetName();
      vSpotsNames{vSpotsNumber} = char(n);
    elseif vFactory.IsSurfaces(vItem)
      vSurfacesNumber = vSurfacesNumber + 1;
      vSurfacesList{vSurfacesNumber} = vFactory.ToSurfaces(vItem);
      n = vSurfacesList{vSurfacesNumber}.GetName();
      vSurfacesNames{vSurfacesNumber} = char(n);
    end
  end

  %% select
  vSurfaceIdx = [];
  while(numel(vSurfaceIdx) ~= 1)
    [vSurfaceIdx, vOk] = listdlg('ListString',vSurfacesNames,...
        'ListSize',[250 150],'Name','Colocalize spots','InitialValue',[1], ...
        'PromptString',{'Please select the SURFACE to use as reference:'});
  end

  vSpotsIdx = [];
  while(numel(vSpotsIdx) ~= 1)
    [vSpotsIdx, vOk] = listdlg('ListString',vSpotsNames,...
        'ListSize',[250 150],'Name','Colocalize spots','InitialValue',[1], ...
        'PromptString',{'Please select the SPOTS to process:'});
  end

%% Params
vQuestion = {sprintf(['Consider the spots  as above / below, if their distance to the surface along the z-axis\n', ...
        'Is more than / less than (default 0, suggested 1 or 2. Values in um):'])};
    vAnswer = inputdlg(vQuestion,'Above/Below Threshold',1,{'1'});
    TH_AB = str2double(vAnswer{1});

    vQuestion = {sprintf(['Exclude the spots whose xy distance is greater than \n', ...
        '(default 100, values in um):'])};
    vAnswer = inputdlg(vQuestion,'Above/Below Threshold',1,{'8'});
    TH_XY = str2double(vAnswer{1});


%% Process
%% spots - surface distance

vXYZ = vSpotsList{vSpotsIdx}.GetPositionsXYZ;
vSurfaces = vSurfacesList{vSurfaceIdx};

surfaceddata = vSurfaces.GetSurfaceData(0);
surf_mask = squeeze(surfaceddata.GetDataFloats());
sx = surfaceddata.GetExtendMinX;
sy = surfaceddata.GetExtendMinY;
sz = surfaceddata.GetExtendMinZ;
ex = surfaceddata.GetExtendMaxX;
ey = surfaceddata.GetExtendMaxY;
ez = surfaceddata.GetExtendMaxZ;

psx = max(1,floor(sx/0.5));
psy = max(1,floor(sy/0.5));
psz = max(1,floor(sz/1.5));
swe = size(surf_mask,1);
she = size(surf_mask,2);
sde = size(surf_mask,3);

vstack_sizes = [swe, she, max(sde, dataset_planes+2)];
stack = zeros(vstack_sizes(1), vstack_sizes(2), vstack_sizes(3), 'single');
stack(psx:psx+swe-1, psy:psy+she-1, psz:psz+sde-1) = squeeze(surf_mask);
%stack(psx:psx+swe-1, psy:psy+she-1, 120-sde+1:end) = squeeze(surf_mask);


[X,Y,Z] = meshgrid(1:vstack_sizes(1), 1:vstack_sizes(2), 1:vstack_sizes(3));
idx = stack > 32768;
xx = X(idx)*0.5;%TODO interface
yy = Y(idx)*0.5; %(in um)
zz = Z(idx)*1.5;


spots_dz = zeros(size(vXYZ,1),1);
spots_dxy = zeros(size(vXYZ,1),1);

TH_Z = 20;
DS_FACTOR = 10;
xx = xx(1:DS_FACTOR:end);
yy = yy(1:DS_FACTOR:end);
zz = zz(1:DS_FACTOR:end);
sXYZ = single([yy,xx,zz]);
h = waitbar(0, 'Computing distances');
for ss = 1:size(vXYZ,1)
    zdx = (vXYZ(ss,2) - xx);
    zdy = (vXYZ(ss,1) - yy);
    zidx = sqrt((zdx.^2) + (zdy.^2))<TH_Z;
    matched_zs = zz(zidx);
    dzs = -(matched_zs - vXYZ(ss,3));
    spots_dz(ss) = mean(dzs);
    d = pdist2(vXYZ(ss,:), sXYZ);
    [md, mdidx] = min(d);
    dz = vXYZ(ss,3) - zz(mdidx);
    dx = vXYZ(ss,2) - xx(mdidx);
    dy = vXYZ(ss,1) - yy(mdidx);
    %spots_dz(ss) = dz;
    spots_dxy(ss) = sqrt((dy.^2)+(dx.^2));
    if(mod(ss,100))
        waitbar(ss/size(vXYZ,1), h);
    end
end
close(h);

%% add spots
vSpots1 = vSpotsList{vSpotsIdx};
vTime1 = vSpots1.GetIndicesT;
vScene = vImarisApplication.GetSurpassScene;
vSpotsGroup = vImarisApplication.GetFactory.CreateDataContainer;
vSpotsGroup.SetName('Above/Below');
%TH_ABOVE_BELOW = 0;
idx_spots_above = (spots_dz < -TH_AB) & (spots_dxy < TH_XY);
idx_spots_below = (spots_dz >= TH_AB) & (spots_dxy < TH_XY);
vNewSpots1 = vImarisApplication.GetFactory.CreateSpots;
vNewSpots1.Set(vXYZ(idx_spots_below, :), vTime1(idx_spots_below), zeros(nnz(idx_spots_below),1)+5);
vNewSpots1.SetName('below');
vNewSpots1.SetColorRGBA(6886); %58906
vSpotsGroup.AddChild(vNewSpots1, -1);

vNewSpots2 = vImarisApplication.GetFactory.CreateSpots;
vNewSpots2.Set(vXYZ(idx_spots_above, :), vTime1(idx_spots_above), zeros(nnz(idx_spots_above),1)+5);
vNewSpots2.SetName('above');
vNewSpots2.SetColorRGBA(58906);
vSpotsGroup.AddChild(vNewSpots2, -1);

vNewSpots1.SetVisible(0);
vNewSpots2.SetVisible(0);
vScene.AddChild(vSpotsGroup, -1);
