function [finalpatch, finalvolfraction,xsep] = CreateBeadedAxonGeometries(rad,amp,targetvolfraction,initialg)
% function to create ply geometry files of 'beaded' axonal geometries
% originally reported in Budde & Frank PNAS 2010.
%
% rad - radius given in microns
% amp - beading amplitude 0 (cylinder) to 1 .
% targetvolfraction - intracellular volume fraction; may be constrained by other 
%                      parameters (packing, g) if unable to fit desired geometries 
%                       max is approximately 0.78 for hexagonal packed
%                       cylinders and high beading amplitudes.
% initialg - initial separation between beads (connecting 'neck')
%             will be modified if surface area is constrained.  default is
%             0.51 and is appropriate for most situtions.
%
%
%  Many other options are allowed below for various situations that have
%  been published (Budde PNAS 2010, Baron Stroke 2014, Skinner NMR Biomed
%  2016) or simply for testing and not reported.
%
%
%  Matt Budde, PhD, Medical College of Wisconsin, 2009-2017
%   mdbudde@mcw.edu
%
% Includes code from other sources, documented and credited where
% appropriate.  
%

%
% Copyright 2017 Matthew Budde
%Permission is hereby granted, free of charge, to any person obtaining a copy of 
% this software and associated documentation files (the "Software"), to use the 
% Software, including without limitation the rights to use, copy, 
% modify, merge the Software. Permission is NOT granted to publish, distribute, sublicense, 
% and/or sell copies of the Software, except with written permission from the copyright holder.
%
% The software is subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all copies 
% or substantial portions of the Software.
%
%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING 
% BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%

%---- Other options that can be modified, but the below values are used for
%most cases.  Some options are legacy from old implementations and others
%are for special cases (cosine and sanorm) that are related to a constrin
%the other input parameters (volumefraction, spacing, amplitude).

nbeads = 3;  % number of beads along the length; not used for hexCropped

layers = 3;  % layers is used for hexblock; not used for hexCropped

packing = 'hexCropped';  % shape of the packing.  Hexagonal (or square) array with an outer perimeter of hexagonal, square, or cropped
                       %'hexhex'  - hexagonal packing with a hexagonal outer border
                       %'hexblock'  - hexagonal packing with a square outer border
                       %'hexCropped' - 'unit' repeating shape cropped.  reduces the number of faces to improve speed
                       % square' - square packing with square outer border
                       
cosinecountour = 1; % 1 yes or 0 no.  Cosine contouer allows for better packing density 
                    % without intersecting adjacent surfaces. No uses the
                    % original shapes in Markin et al 1999 and intersecting
                    % geometries cannot be gauranteed. 
                    
sanorm = 0;  % 1 yes or 0 no; normalize surface area of beaded shapes to equivilant cylinder, only for non-cosine contours

closedsurface = 0;  %includes a commment in ply file as closed surface; used in camino for intra/extra calculations
                    % if not included, camino uses a ray-tracing method, but very slow.
                    % not including it greatly improves speed.
                    % In my edits (local) to the camino code, the normals
                    % of each face point to the extracellular space which
                    % camino can use to determine if spins are intra/extra
                    
vertexDiv = 32; %number of vertices, passed to makeUnd, 25 = smooth; 9= not smooth, but faster simulations.
                %higher number gives smoother surfaces, at the cost of increased simulation time.
                % with the hexCropped option, this can be higher (32) since
                % the total number of faces is considerably reduced.

global pointTol;
pointTol = 1e-12; %tolerance to find similar vertex points, same units as radius 

cutends = 0;  % if using hexblock or other, only clip the ends.

dosurface = 0;  % with clipped vertices, surface or patch is broken due to (potential) different #s of vertices with each face
docontour = 0;
doplyfile = 1;



%--------------end of modifiable parameters ---------------%
   
if nargin < 3
    targetvolfraction = .6;
end
if nargin < 1
    rad = 1.0;
end
if nargin < 2
    amp = 0;%sqrt(3)-1;% for hex: sqrt(3)-1;
end
if nargin < 4
    initialg = 0.51; 
end

% if surface area is to be normalized, need to calculate the spacing (g) that will account for the beading shape
%  generally, this is not necessary since the differences are small.  It
%  was important to show that the beaded shape can encompass a larger
%  volume within the same surface area.
if exist('fsolve','file') && sanorm == 1
    [surfa,vol,len,newrad,intrad,newg] = beadParams(rad,amp,initialg,1)
    fprintf('fsolve found.  Initial: Rad= %.3f  G= %.3f   Final: Rad= %.3f  G=%.3f\n',rad,initialg,newrad,newg);
    fprintf('Final integrated rad=%.3f\n',intrad);
else
    disp('fsolve not on path or sanorm = false, using initial radius and g.n');
    newrad = rad;
    newg=initialg;
end

%  Create the surface of a single bead with the given parameters.
[xx,yy,zz,dd,contour,yym] = makeUnduloid(rad,newrad,amp,newg,nbeads,0,cosinecountour,vertexDiv);


% setup ply structure
finalply.vertex.x = [];
finalply.vertex.y = [];
finalply.vertex.z = [];
finalply.face.vertex_indices = {};

[aa,vv,ll,rav1,rav2] = beadParams(xx(1,:),zz(1,:));

leng=ll/nbeads;
centerpt = floor(size(xx,1)/2);
if dosurface
    sf = figure;
    axis vis3d;
end
if docontour
    cf = figure;
end

switch packing
    case 'hexhex'
        rows = layers*2 - 1;
        ptsperrow = rows - abs(-(layers-1):(layers-1));
        xoffsetmat = mod(ptsperrow-layers,3)+1;

        %normrad(ampindex,gindex)
        fullvol = (vv/ll/2);
        
        xsep = sqrt(fullvol/(targetvolfraction*sqrt(3))); % equal to .80 volume fraction for cylinders, hexagonal spacing. side of equilateral triangle=2x.
        	rmax = max(zz(1,:))
        rmin = min(zz(1,:));
        if xsep < (rmax + rmin)/2
        	xsep=(rmax + rmin)/2
        end
        % don't allow adjacent cylinders to overlap at max radius
        if xsep*sqrt(3) < rmax
        	xsep = rmax/sqrt(3)
        end  
        
        finalvolfraction=(vv/ll/2)/(xsep^2*sqrt(3));
        fprintf('Calculated Volumefraction=%.2f\n',(vv/ll/2)/(xsep^2*sqrt(3)));
        
    case 'hexblock'
        rows = layers;
        ptsperrow = rows*ones(1,rows);
        xoffsetmat = mod(repmat([0 1],[1 layers]),2)*2+1;

        %normrad(ampindex,gindex)
        fullvol = (vv/ll/2);
        
        %xsep = sqrt(pi*rad^2/2/(targetvolfraction*sqrt(3))); % equal to .80 volume fraction for cylinders, hexagonal spacing. side of equilateral triangle=2x.
        xsep = sqrt(fullvol/(targetvolfraction*sqrt(3))); % side of equilateral triangle=2x.
        % try to increase volfraction with increasing amplitude
        %xsep = xsep + (amp)*(newrad-xsep)*2        
        % but don't go with a vol fraction smaller than target.
        %finalvolfraction=(vv/ll/2)/(xsep^2*sqrt(3));
        %if finalvolfraction<targetvolfraction
           % 	xsep = sqrt(fullvol/(targetvolfraction*sqrt(3))) % side of equilateral triangle=2x.
        %end
        % don't allow adjacent cylinders to intersect at min and max radius points, still could intersect along length at abutting beads. (much harder problem).
        rmax = max(zz(1,:));
        rmin = min(zz(1,:));
        if xsep < (rmax + rmin)/2
            xsep=(rmax + rmin)/2;
        end
        % don't allow adjacent cylinders to overlap at max radius
        if xsep*sqrt(3) < rmax
            xsep = rmax/sqrt(3);
        end       

        finalvolfraction=(vv/ll/2)/(xsep^2*sqrt(3));
        fprintf('Volumefraction=%.2f\n',finalvolfraction);    
        
    case 'hexCropped' %same as hexblock but cropped to a unit geometry.  3 rows and 5 layers are needed
        rows = 3; %layers
        yrows = 5; %layers
        ptsperrow = yrows*ones(1,rows);
        xoffsetmat = mod(repmat([0 1],[1 rows]),2)*2+1;

        %normrad(ampindex,gindex)
        fullvol = (vv/ll/2);
        
        %xsep = sqrt(pi*rad^2/2/(targetvolfraction*sqrt(3))); % equal to .80 volume fraction for cylinders, hexagonal spacing. side of equilateral triangle=2x.
        xsep = sqrt(fullvol/(targetvolfraction*sqrt(3))); % side of equilateral triangle=2x.
        % try to increase volfraction with increasing amplitude
        %xsep = xsep + (amp)*(newrad-xsep)*2        
        % but don't go with a vol fraction smaller than target.
        %finalvolfraction=(vv/ll/2)/(xsep^2*sqrt(3));
        %if finalvolfraction<targetvolfraction
        %    	xsep = sqrt(fullvol/(targetvolfraction*sqrt(3))) % side of equilateral triangle=2x.
        %end
        % don't allow adjacent cylinders to intersect at min and max radius points, still could intersect along length at abutting beads. (much harder problem).
        rmax = max(zz(1,:));
        rmin = min(zz(1,:));
        if xsep < (rmax + rmin)/2
            xsep=(rmax + rmin)/2;
        end
        % don't allow adjacent cylinders to overlap at max radius
        if xsep*sqrt(3) < rmax
            xsep = rmax/sqrt(3);
        end       

        finalvolfraction=(vv/ll/2)/(xsep^2*sqrt(3));
        fprintf('Volumefraction=%.2f\n',finalvolfraction);
        
    case 'square'
        rows = layers;
        ptsperrow = ones(1,rows)*rows;
        xoffsetmat = mod(abs(-(layers-1):(layers-1))-layers,2)

        xsep = sqrt(pi*newrad^2/targetvolfraction);
    	
        rmax = max(zz(1,:));
        rmin = min(zz(1,:));
        
        if xsep*sqrt(2) < (rmax)*2 
            xsep=(rmax)*2/sqrt(2);
        end
        if xsep < (rmax + rmin)
            xsep = rmax + rmin;
        end
        
        finalvolfraction=(vv/ll)/(xsep^2);
        fprintf('Volumefraction=%.2f\n',finalvolfraction);
end

for jj=1:rows
    for ii=1:ptsperrow(jj)
       
        switch packing
            case 'hexhex'
                xoffset = (leng/3)*mod((xoffsetmat(jj)+mod(ii-1,3)),3);
                yoffset = (ptsperrow(jj)/2-ii+1/2)*xsep*2;
                zoffset = (xsep*sqrt(3))*(rows/2-jj+1/2);
            case 'hexblock'
                xoffset = (leng/3)*mod((xoffsetmat(jj)+mod(ii-1,3)),3);
                yoffset = (ptsperrow(jj)/2-ii+mod(jj,2)*1/2)*xsep*2;
                zoffset = (xsep*sqrt(3))*(rows/2-jj+1/2);            
            case 'hexCropped' %same as hexblock but cropped to a 'unit' repeating geometry
                xoffset = (leng/3)*mod((xoffsetmat(jj)+mod(ii-1,3)),3);
                yoffset = (ptsperrow(jj)/2-ii+mod(jj,2)*1/2)*xsep*2;
                zoffset = (xsep*sqrt(3))*(rows/2-jj+1/2);
            case 'square'
                xoffset = (leng/2)*mod((xoffsetmat(jj)+mod(ii-1,2)),2);
                %xoffset = (leng/3)*mod(mod(jj,3)*mod(ii,3),3);
                yoffset = ((ii))*xsep; % to center along y add this: - xsep*(rows-1)
                zoffset = ((jj))*xsep;
        end
            solidcolor = mod(xoffsetmat(jj)+mod(ii-1,3),3);
            %set(sf,'NextPlot','add'); %hold(sf,'on');
        
            xxkeep = xx+xoffset;
            yykeep = yy+yoffset;
            zzkeep = zz+zoffset;
            yymkeep = yym;
            
            
        if dosurface
            figure(sf);
            hold on
          
            surf(xxkeep,yykeep,zzkeep,yymkeep,'LineStyle','none','FaceColor','interp','DiffuseStrength',.6,'AmbientStrength',.6,'SpecularStrength',.9);
        end
        
        %set(cf,'NextPlot','add'); %hold(cf,'on');
        if docontour
            figure(cf);
            hold on
            contour(yy+yoffset,zz+zoffset,xx+xoffset,[0 0]);
        end
        
        if doplyfile
            if strcmp(packing,'hexCropped')
                            
                xmin = -1*leng/3;
                xmax = leng/3*2;
                zmax = (xsep*sqrt(3))*(3/2-1+1/2);
                zmin = (xsep*sqrt(3))*(3/2-3+1/2);
                ymax = (ptsperrow(2)/2-1+mod(2,2)*1/2)*xsep*2;
                ymin = (ptsperrow(2)/2-4+mod(2,2)*1/2)*xsep*2;

                tempcoords = vertices2plycoords(xxkeep,yykeep,zzkeep,[xmin xmax],[ymin ymax],[zmin zmax]);
            else
                if cutends == 1  %only clip ends 
                    xmin = -1*leng/3;
                    xmax = leng/3*2;
                    tempcoords = vertices2plycoords(xxkeep,yykeep,zzkeep,[xmin xmax]);
                else %no clipping
                    tempcoords = vertices2plycoords(xxkeep,yykeep,zzkeep);
                end
            end
            
            if isfield(tempcoords.face,'vertex_indices')
            vertices = size(finalply.vertex.x,2);
            
            % ply files are zero indexed.
            tempindices = tempcoords.face.vertex_indices;
            for tind=1:size(tempindices,2)
               % if isempty(finalply.face.vertex_indices)
                %    tempindices{tind} = tempindices{tind} - 1;
                %else
                
                % flip vertices indices here to have normals point to the
                % extracellular side.  These are used in my modified camino
                % code to indicate intra/extra space instead of the ray
                % tracing algorithm in the conventional camino version.
                    tempindices{tind} = tempindices{tind}(end:-1:1) + vertices - 1;
                %end
            end
            finalply.face.vertex_indices = cat(2, finalply.face.vertex_indices, tempindices);

            % put ply files in meters for camino standard units.
            finalply.vertex.x = cat(2, finalply.vertex.x, (tempcoords.vertex.x)./1e6);
            finalply.vertex.y = cat(2, finalply.vertex.y, (tempcoords.vertex.y)./1e6);
            finalply.vertex.z = cat(2, finalply.vertex.z, (tempcoords.vertex.z)./1e6);
            end
        end
    end
end

%for matlab patch command
finalpatch(1).Vertices = [finalply.vertex.x', finalply.vertex.y', finalply.vertex.z'];
%add 1 for matlab indexing
%finalpatch.Faces = cell2mat(finalply.face.vertex_indices') + 1;


    if dosurface
        figure(sf); 
        axis equal
        %axis vis3d
        %    view(3)
        light
        lighting gouraud
    end
    
if doplyfile
    fprintf('Writing Ply file');
    if strcmp(packing,'hexCropped')
        plyname = sprintf('RectCropped_%1.1fRad_%1.2fAmp_%1.2fG_%1.2fVolfract.ply',rad,amp,initialg,finalvolfraction);
        plywrite(finalply,plyname,'ascii','double',closedsurface);
    else
        plyname = sprintf('hexpackInMeters_%1.1fRad_%1.2fAmp_%1.2fG_%1.0flayers_%1.2fVolfract.ply',rad,amp,initialg,layers,finalvolfraction);
        plywrite(finalply,plyname,'ascii','double',closedsurface);
    end
    fprintf('Ply file written: %s',plyname);
end


function Data = vertices2plycoords(x,y,z,xminmax,yminmax,zminmax)
    
    global pointTol;
% generate ploygon faces from the vertices.
    Data = struct('vertex',[],'face',[]);
    
    Data.vertex.x = x(:)';
    Data.vertex.y = y(:)';
    Data.vertex.z = z(:)';
    
    if ~exist('xminmax','var')
        xminmax = [min(Data.vertex.x), max(Data.vertex.x)];
        %xminmax = xminmax + [-xminmax(1)*1.2, xminmax(2)*1.2] % just to be sure nothing is clipped.
    end
    if ~exist('yminmax','var')
        yminmax = [min(Data.vertex.y), max(Data.vertex.y)];
        %yminmax = yminmax + [-yminmax(1)*1.2, yminmax(2)*1.2]
    end
    if ~exist('zminmax','var')
        zminmax = [min(Data.vertex.z), max(Data.vertex.z)];
        %zminmax = zminmax + [-zminmax(1)*1.2, zminmax(2)*1.2]
    end
    
    
    % add all faces to struct
    faceind = 1;
    for si=1:size(x,1)-1
        for sj=1:size(x,2)-1
            sx = size(x);
            facedata = [sub2ind(sx,si,sj), sub2ind(sx,si+1,sj), sub2ind(sx,si+1,sj+1), sub2ind(sx,si,sj+1)];
            Data.face.vertex_indices{faceind} = facedata;
            faceind = faceind + 1;
        end
    end
    
    % clipping
    keeplist = []; rmlist = []; modlist = [];
    for kk=1:length(Data.face.vertex_indices)
        currFace = Data.face.vertex_indices{kk};
        xverts = x(currFace);        
        yverts = y(currFace);        
        zverts = z(currFace);
        
        xin = xverts >= xminmax(1) & xverts <= xminmax(2);
        yin = yverts >= yminmax(1) & yverts <= yminmax(2);
        zin = zverts >= zminmax(1) & zverts <= zminmax(2);
        
        %
        if all(xin) && all(yin) && all(zin)
        % keep if all vertices in clipping area
            keeplist = [keeplist kk];
        elseif ~any(xin) || ~any(yin) || ~ any(zin)
        % remove if no vertices in clipping area
            rmlist = [rmlist kk];
        else  
        % modify if face crosses clipping plane 
            modlist = [modlist kk];
        end       
    end
    NewData.face.vertex_indices = Data.face.vertex_indices(keeplist);
    

    if 1  %modify clipped vertices to ensure points are consistent on each side 
    disp('Modifying clipped vertices');
    for mm=[modlist]
        currFace = Data.face.vertex_indices{mm};
        xverts = x(currFace);        
        yverts = y(currFace);        
        zverts = z(currFace);
        
        verts = [xverts; yverts; zverts]';
        xyzminmax = [xminmax; yminmax; zminmax];
        
        [clippedyz, newvertinds] = sutherlandHodgman3d(verts,xyzminmax,currFace);
        
        if ~isempty(clippedyz)
        clipx = clippedyz(:,1);
        clipy = clippedyz(:,2);
        clipz = clippedyz(:,3);
        
        newfacevertices = [];
        for nv=1:length(clipy)
            if newvertinds(nv) == 0
                vind = length(Data.vertex.x);
                % find any vertex that might already exist, within a
                % tolerance
                closestDist = [];
                closestVertex = [];
                for np=1:vind
                    vertexdist = sqrt((clipx(nv) - Data.vertex.x(np)).^2 + (clipy(nv) - Data.vertex.y(np)).^2 + (clipz(nv) - Data.vertex.z(np)).^2);
                    if vertexdist < pointTol & vertexdist < closestDist
                        closestVertex = np;
                        closestDist = vertexdist;
                    end
                end
                if isempty(closestVertex) %point is new
                    Data.vertex.x(vind + 1) = clipx(nv);
                    Data.vertex.y(vind + 1) = clipy(nv);
                    Data.vertex.z(vind + 1) = clipz(nv);
                    newfacevertices(end + 1) = vind + 1;
                else  % point already exists
                    if any(newfacevertices == closestVertex)
                        % face already contains this point
                        disp('Duplicate vertex in face. removed.');
                        continue;
                    else
                        newfacevertices(end+1) = closestVertex;
                    end
                end
            else %not modified vertex
                newfacevertices(end+1) = newvertinds(nv); 
            end
        end
        if length(newfacevertices) < 3 
            disp('Skipped face with < 3 verts.');
        else
            NewData.face.vertex_indices{end+1} = newfacevertices;
        end
        end
           
    end %for modified faces
    end %if
    
    Data.face.vertex_indices = NewData.face.vertex_indices;
    
    if 0  % not completed yet or tested
    disp('Removing unused vertices');
    % remove unused vertices
    newind = 1;
    for ii=1:length(Data.vertex.x)
        curridx = [Data.vertex.x(ii) Data.vertex.y(ii) Data.vertex.z(ii)];
        for jj=1:length(Data.face.vertex_indices)
            for kk=1:4
                if Data.face.vertex_indices{jj}(kk) == ii
                    NewData.vertex.x(newind) = curridx(1);
                    NewData.vertex.y(newind) = curridx(2);
                    NewData.vertex.z(newind) = curridx(3);
                    NewData.face.vertex_indices{jj}(kk) = newind;
                    newinds = newinds + 1;
                end
            end
        end
    end
    end % if 0
    
%end %vertices2plycoords

    
function [clippedPolygon, clippedVerts] = sutherlandHodgman3d(subjectPolygon,xyzminmax,vertinds)
 % Sutherland-Hodgman Algorithm;  implements cropping of meshes for the
 % hexCropped 
 
    clippedPolygon = subjectPolygon;
    clippedVerts = vertinds;
    
    xyzplane = [1 1 2 2 3 3];
    minmaxplane = [1 2 1 2 1 2];
    for pp=1:6 % number of planes (6) clipVertex = (1:numVerticies)
      
       
       clipPlaneVertex = zeros(1,3);
       clipPlaneVertex(xyzplane(pp)) = xyzminmax(xyzplane(pp),minmaxplane(pp));
       clipPlaneNormal = zeros(1,3);
       clipPlaneNormal(xyzplane(pp)) = -xyzminmax(xyzplane(pp),minmaxplane(pp));
        
        inputList = clippedPolygon;
        inputverts = clippedVerts;
    
        clippedPolygon = []; clippedVerts = [];
        if ~isempty(inputList),
            previousVertex = inputList(end,:);
        end
 
        for subjectVertex = (1:size(inputList,1))
 
            if  insidePlane(inputList(subjectVertex,:),clipPlaneNormal,clipPlaneVertex) 
                
                if not(insidePlane(previousVertex,clipPlaneNormal,clipPlaneVertex))   
                    [I, check] = plane_line_intersect(clipPlaneNormal, clipPlaneVertex, previousVertex, inputList(subjectVertex,:));
                    if check == 1
                        clippedPolygon(end+1,1:3) = I;
                        clippedVerts(end+1) = 0;
                    end
                end
                clippedPolygon(end+1,1:3) = inputList(subjectVertex,:);
                clippedVerts(end+1) = inputverts(subjectVertex);
 
            elseif( insidePlane(previousVertex,clipPlaneNormal, clipPlaneVertex) )
               [I, check] = plane_line_intersect(clipPlaneNormal, clipPlaneVertex, previousVertex, inputList(subjectVertex,:));
                if check == 1 
                    clippedPolygon(end+1,1:3) = I;
                    clippedVerts(end+1) = 0;
                end                
            end
 
            previousVertex = inputList(subjectVertex,:);
           
        end %for subject verticies                
    end %for boundary verticies
%end %sutherlandHodgman

    

function in = insidePlane(point,planenormal, planepoint)

    dotVector = dot(point-planepoint,planenormal);

    if ( dotVector >= 0 )
        in = true;
    else
        in = false;
    end
 

function [I,check]=plane_line_intersect(n,V0,P0,P1)
%plane_line_intersect computes the intersection of a plane and a segment(or
%a straight line)
% Inputs: 
%       n: normal vector of the Plane 
%       V0: any point that belongs to the Plane 
%       P0: end point 1 of the segment P0P1
%       P1:  end point 2 of the segment P0P1
%
%Outputs:
%      I    is the point of interection 
%     Check is an indicator:
%      0 => disjoint (no intersection)
%      1 => the plane intersects P0P1 in the unique point I
%      2 => the segment lies in the plane
%      3=>the intersection lies outside the segment P0P1
%
% Example:
% Determine the intersection of following the plane x+y+z+3=0 with the segment P0P1:
% The plane is represented by the normal vector n=[1 1 1]
% and an arbitrary point that lies on the plane, ex: V0=[1 1 -5]
% The segment is represented by the following two points
% P0=[-5 1 -1]
%P1=[1 2 3]   
% [I,check]=plane_line_intersect([1 1 1],[1 1 -5],[-5 1 -1],[1 2 3]);

%This function is written by :
%                             Nassim Khaled
%                             Wayne State University
%                             Research Assistant and Phd candidate
%If you have any comments or face any problems, please feel free to leave
%your comments and i will try to reply to you as fast as possible.
I=[0 0 0];
u = P1-P0;
w = P0 - V0;
D = dot(n,u);
N = -dot(n,w);
check=0;
if abs(D) < 1e-12        % The segment is parallel to plane
        if N == 0           % The segment lies in plane
            check=2;
            return
        else
            check=0;       %no intersection
            return
        end
end

%compute the intersection parameter
sI = N / D;
I = P0+ sI.*u;

if (sI < 0 || sI > 1)
    check= 3;          %The intersection point  lies outside the segment, so there is no intersection
else
    check=1;
end



function [xm,ym,zm, Data, contour,yym] = makeUnduloid(r,rav,Amp,g,nbeads,doplot,coscontour,vertexDiv)
% run the calculation and plotting of 3d unduloids.
%  r - initial radius, in microns (for neck)
% rav - initial bead radius
%  Amp - amplitude of oscillation
%  g   - bead length/cylinder length
% nbeads  - number of repeatings beads
% doplot - show plot (yes/no)
% coscontour - use a cosine beaded contour instead of the shape in Budde PNAS 2010
% vertexDiv - number of vertices/faces per bead


    if ~exist('doplot','var')
        doplot = 0;
    end
    if ~exist('coscontour','var')
        coscontour = 0;
    end
    if ~exist('doplycoords','var')
        doplycoords = 0;
    else
        doplycoords = 1;
    end
    if ~exist('vertexDiv','var')
        vertexDiv = 25;  
    end
    % compute the curve once and just repeat it if additional length is
    % desired. to be implemented.
    %clear xx yy
    
    jj = -pi:2*pi/vertexDiv:pi;
    xx = zeros(size(jj));
    yy = zeros(size(jj));


    for j=1:size(jj,2)
        if ~coscontour
            [xx(j), yy(j)] = unduloidCoords(rav,Amp,jj(j));
        else
            length = unduloidCoords(rav,0,pi);
        	xx(j) = jj(j)*(length/pi);
        	yy(j) = rav*Amp*cos(jj(j))+rav;
        end
    end

    if g>0
        %generate compacted cylinder, 
        %length is proportinal to unbeaded cylinder length, not of beaded cylinder.
        [minx, miny] = unduloidCoords(r,0,-pi);
        [maxx, dumy] = unduloidCoords(r,0,pi);
        [dumx, maxy] = unduloidCoords(r,Amp,0);

        % or proportional to beaded length
        %[minx, dumy] = unduloid(rav,Amp,-pi);
        %[maxx, dumy] = unduloid(rav,Amp,pi);

        
        
     %   if ~coscontour
            cylLen = (maxx - minx).*g;
            cylDiv = (cylLen/2) / floor((cylLen/2)/(2*pi/vertexDiv*rav));

            if cylDiv >= cylLen/2
                cylDiv = cylLen/2;
            end
            xx2 = 0:cylDiv:cylLen/2;
            yy2 = xx2*0 + yy(1);
            %build one bead-cyl pair and duplicate them
            minusx = min(xx) - cylLen/2 + xx2;% - cylDiv;
            %minusx = min(xx)-(cylLen)/2 + xx2;
            plusx = max(xx)+xx2;% + cylDiv;
            newx = [minusx,xx,plusx];
            newy = [yy2, yy, yy2];
%         else
%             
%             %%NEED TO CHECK THIS 4/2014
%             
%             % include part of the g into the bead maximum.
%             cylLen = (maxx - minx).*g;
%             cylDiv = (cylLen/3) / floor((cylLen/3)/(2*pi/vertexDiv*rav));
% 
%             if cylDiv >= cylLen/3
%                cylDiv = cylLen/3;
%             end
%             xx2 = 0:cylDiv:cylLen/3;
%             yy2 = xx2*0 + yy(1);
%             
%             halfx = floor(size(xx,2)/2);
%             xx3 = 0:cylDiv:cylLen/3 + xx(halfx);
%             yy3 = xx3*0 + max(yy);
%             
%             if ~isempty(xx3)
%             % rebuild bead portion first
%             xx = [xx(1:halfx), xx3, xx(halfx+1:end)+xx3(end)];
%             yy = [yy(1:halfx), yy3, yy(halfx+1:end)];
%             %minusx = min(xx)-(cylLen)/2 + xx2;
%             end
%             
%             %build one bead-cyl pair and duplicate them
%             minusx = min(xx) - cylLen/3 + xx2;% - cylDiv;
%             %minusx = min(xx)-(cylLen)/2 + xx2;
%             plusx = max(xx)+xx2; % + cylDiv;
%             newx = [minusx,xx,plusx];
%             newy = [yy2, yy, yy2];
%         end
        %plot(xx,yy,'ko');
        %hold on;
        %plot(minusx,yy2,'ro');
        %plot(plusx,yy2,'bo');
    else
        cylLen = 0;
        if nbeads>1
            newx = xx(1:end-1);
            newy = yy(1:end-1);
        else
            newx = xx;
            newy = yy;
        end
    end

    %concatenate the shapes together
    a = repmat(0:nbeads-1,[size(newx,2) 1]);
    a = a(:)';
    pairLen = max(xx) - min(xx) + cylLen;
    newx=repmat(newx,[1 nbeads]) + a.*pairLen;
    newy=repmat(newy,[1 nbeads]);
    
    % adding a cylinder to each end
    %newx = [newx, xx2 + (nbeads-1)*pairLen + max(xx)];
    %newy = [newy, yy2];
    
    % center the final structure
xx = newx - min(newx) - (max(newx) - min(newx))/2;
yy = newy;


%just in case we have some overlapping points, which could cause
%overlapping triangles
apts = unique([xx;yy]','rows');
xx = apts(:,1)';
yy = apts(:,2)';
contour.xx = xx;
contour.yy = yy;

% generate surface vertices.
radialDiv = 2*pi/vertexDiv;
pp = 0:radialDiv:(2*pi);
xm = repmat(xx,[size(pp,2),1]);
yym = repmat(yy,[size(pp,2),1]);
pm = repmat(pp,[size(xm,2),1]);
ym = yym.*sin(pm');
zm = yym.*cos(pm');

if doplot
    surf(xm,ym,zm,yym);
    view(0,0);
    caxis([0 rav*2]);
    axis equal
    %axis vis3d
    %max length from cylinder
    %xlim([min(xm(1,:))*1.25 max(xm(1,:))*1.25]);
    %ylim([-rav*2 rav*2]);
    %zlim([-rav*2 rav*2]);
    %set(get(gca,'Children'),'LineStyle','none');
    figure(gcf);
    pause(.01);
end

% generate ploygon faces from the vertices.
    Data = struct('vertex',[],'face',[]);
if doplycoords
    %
    Data.vertex.x = xm(:)';%;
    Data.vertex.y = ym(:)';%;
    Data.vertex.z = zm(:)';%;
    faceind = 1;
    for ii=1:size(xm,1)-1
        for jj=1:size(xm,2)-1
            sx = size(xm);
            
            % ply files are 0 indexed, not 1.
            
            %4 sided polygon
            facedata = [sub2ind(sx,ii,jj), sub2ind(sx,ii+1,jj), sub2ind(sx,ii+1,jj+1), sub2ind(sx,ii,jj+1)];
            Data.face.vertex_indices{faceind} = facedata - 1;
            faceind = faceind + 1;
            
            %triangles
            %facedata1 = [sub2ind(sx,ii,jj), sub2ind(sx,ii+1,jj), sub2ind(sx,ii+1,jj+1)];
            %facedata2 = [sub2ind(sx,ii,jj), sub2ind(sx,ii+1,jj+1), sub2ind(sx,ii,jj+1)];
            %Data.face.vertex_indices{faceind} = facedata1 - 1;
            %faceind = faceind + 1;
            %Data.face.vertex_indices{faceind} = facedata2 - 1;
            %faceind = faceind + 1;
            
        end
    end
end

function [x,y] = unduloidCoords(rav,ao,theta)
% function [x,y] = unduloid(rav,a,theta)
%  Computes coordinates of a beaded cylinder
%   with constant mean curvature (Delauney unduloid)
%   Inputs:
%    rav    -   average radius of cylinder
%    ao      -   amplitude of oscillation (0-1, 0=no beading)
%    theta  -   rotation of ellipsoid
%
%   Outputs:
%    x,y    -   cartesian coordinates for single point
%
%   1/13/09, Matt Budde

if nargin < 3
    error('3 arguments needed');
end

a = rav; 
b = sqrt(a^2*(1-ao^2));

q = quadgk(@(x)integ(x,a,b),0,theta);
x = q + ((a + sqrt(a^2 - b^2)*cos(theta))*sqrt(a^2 - b^2)*sin(theta))./sqrt(a^2*sin(theta)^2 + b^2*cos(theta)^2);
y = b*(a + sqrt(a^2 - b^2)*cos(theta))/sqrt(a^2*sin(theta)^2 + b^2*cos(theta)^2);

function gg = integ(x,a,b)
gg = sqrt(a.^2.*sin(x).^2 + b.^2.*cos(x).^2);



function plywrite(Elements,Path,Format,Str,closedsurface)
%PLYWRITE  Write 3D data as a PLY file.
%   PLYWRITE(DATA,FILENAME) writes the structure DATA as a binary 
%   PLY file.  Every field of DATA is interpreted as an element
%   and every subfield as an element property.  Each subfield of
%   property data must either be an array or a cell array of 
%   arrays.  All property data in an element must have the same
%   length.
%
%   A common PLY data structure has the following fields:
%      DATA.vertex.x = x coordinates, [Nx1] real array
%      DATA.vertex.y = y coordinates, [Nx1] real array
%      DATA.vertex.z = z coordinates, [Nx1] real array
%
%      DATA.face.vertex_indices = vertex index lists, 
%         an {Mx1} cell array where each cell holds a one-
%         dimesional array (of any length) of vertex indices.
%   Some other common data fields:
%      DATA.vertex.nx = x coordinate of normal, [Nx1] real array
%      DATA.vertex.ny = y coordinate of normal, [Nx1] real array
%      DATA.vertex.nz = z coordinate of normal, [Nx1] real array
%
%      DATA.edge.vertex1 = index to a vertex, [Px1] integer array
%      DATA.edge.vertex2 = second vertex index, [Px1] integer array
%   Many other fields and properties can be added.  The PLY format 
%   is not limited to the naming in the examples above -- they are
%   only the conventional naming.
%
%   PLYWRITE(DATA,FILENAME,FORMAT) write the PLY with a specified 
%   data format, where FORMAT is
%      'ascii'                  ASCII text data
%      'binary_little_endian'   binary data, little endian
%      'binary_big_endian'      binary data, big endian (default)
%
%   PLYWRITE(DATA,FILENAME,FORMAT,'double') or
%   PLYWRITE(DATA,FILENAME,'double') write floating-point data as
%   double precision rather than in the default single precision.
%
%   Example:
%   % make a cube
%   clear Data;
%   Data.vertex.x = [0;0;0;0;1;1;1;1];
%   Data.vertex.y = [0;0;1;1;0;0;1;1];
%   Data.vertex.z = [0;1;1;0;0;1;1;0];
%   Data.face.vertex_indices = {[0,1,2,3],[7,6,5,4], ...
%         [0,4,5,1],[1,5,6,2],[2,6,7,3],[3,7,4,0]};
%   plywrite(Data,'cube.ply','ascii');
%
%   See also: PLYREAD

% Pascal Getreuer 2004
% Matt Budde 2009, added option for closedsurface

if nargin < 4
   Str = '';
   
   if nargin < 3
      Format = 'binary_big_endian';
   elseif strcmpi(Format,'double')
      Str = 'double';
      Format = 'binary_big_endian';
   end
end

if ~exist('closedsurface','var')
    closedsurface = 0;
end

[fid,Msg] = fopen(Path,'wt');

if fid == -1, error(Msg); end

PlyTypeNames = {'char','uchar','short','ushort','int','uint','float','double', ...
   'char8','uchar8','short16','ushort16','int32','uint32','float32','double64'};
FWriteTypeNames = {'schar','uchar','int16','uint16','int32','uint32','single','double'};
MatlabTypeNames = {'int8','uint8','int16','uint16','int32','uint32','single','double'};
PrintfTypeChar = {'%d','%u','%d','%u','%d','%u','%-.6e','%-.14e'};
IntegerDataMin = [-128,0,-2^15,-2^31,0];
IntegerDataMax = [127,255,2^16-1,2^31-1,2^32-1];

%%% write PLY header %%%
fprintf(fid,'ply\nformat %s 1.0\ncomment created by MATLAB plywrite\n',Format);
if closedsurface
fprintf(fid,'comment closed surface\n');
end
ElementNames = fieldnames(Elements);
NumElements = length(ElementNames);
Data = cell(NumElements,1);

for i = 1:NumElements
   eval(['tmp=isa(Elements.',ElementNames{i},',''struct'');']);
   
   if tmp
      eval(['PropertyNames{i}=fieldnames(Elements.',ElementNames{i},');']);
   else
      PropertyNames{i} = [];
   end
   
   if ~isempty(PropertyNames{i})
   	eval(['Data{i}{1}=Elements.',ElementNames{i},'.',PropertyNames{i}{1},';']);
      ElementCount(i) = prod(size(Data{i}{1}));
      Type{i} = zeros(length(PropertyNames{i}),1);
   else
      ElementCount(i) = 0;
   end
   
   fprintf(fid,'element %s %u\n',ElementNames{i},ElementCount(i));
   
   for j = 1:length(PropertyNames{i})
      eval(['Data{i}{j}=Elements.',ElementNames{i},'.',PropertyNames{i}{j},';']);
      
      if ElementCount(i) ~= prod(size(Data{i}{j}))
      	fclose(fid);
         error('All property data in an element must have the same length.');
      end
      
      if iscell(Data{i}{j})
         Type{i}(j) = 9;
         Data{i}{j} = Data{i}{j}{1};
      end
      
      for k = 1:length(MatlabTypeNames)
      	if isa(Data{i}{j},MatlabTypeNames{k})
         	Type{i}(j) = Type{i}(j) + k;
	         break;
         end
      end
      
      if ~rem(Type{i}(j),9)
         fclose(fid);
         error('Unsupported data structure.');
      end
      
      % try to convert float data to integer data
      if Type{i}(j) <= 8 			% array data
         if any(strcmp({'single','double'},MatlabTypeNames{Type{i}(j)}))
            if ~any(floor(Data{i}{j}) ~= Data{i}{j})		% data is integer
               MinValue = min(min(Data{i}{j}));
               MaxValue = max(max(Data{i}{j}));
               
               % choose smallest possible integer data format
               tmp = max(min(find(MinValue >= IntegerDataMin)),min(find(MaxValue <= IntegerDataMax)));
               
               if ~isempty(tmp)
                  Type{i}(j) = tmp;
               end
            end
         end
      else								% cell array data
         eval(['Data{i}{j}=Elements.',ElementNames{i},'.',PropertyNames{i}{j},';']);
         tmp = 1;
         
         for k = 1:prod(size(Data{i}{j}))
            tmp = tmp & all(floor(Data{i}{j}{k}) == Data{i}{j}{k});
	      end
         
         if tmp		% data is integer
	         MinValue = inf;
   	      MaxValue = -inf;
         
      	   for k = 1:prod(size(Data{i}{j}))
         	   MinValue = min(MinValue,min(Data{i}{j}{k}));
            	MaxValue = max(MaxValue,max(Data{i}{j}{k}));
	         end
            
            % choose smallest possible integer data format
            tmp = max(min(find(MinValue >= IntegerDataMin)),min(find(MaxValue <= IntegerDataMax)));
            
            if ~isempty(tmp)
               Type{i}(j) = tmp + 9;
            end
         end
      end
      
      % convert double to single if specified
      if rem(Type{i}(j),9) == 8 & ~strcmpi(Str,'double')
      	Type{i}(j) = Type{i}(j) - 1;
      end
      
      if Type{i}(j) <= 8
         fprintf(fid,'property %s %s\n',PlyTypeNames{Type{i}(j)},PropertyNames{i}{j});
      else
         fprintf(fid,'property list uchar %s %s\n',PlyTypeNames{Type{i}(j)-7},PropertyNames{i}{j});
      end
   end
end

fprintf(fid,'end_header\n');

switch Format
case 'ascii'
   Format = 0;
case 'binary_little_endian'
   fclose(fid);
   fid = fopen(Path,'a','ieee-le');
   Format = 1;
case 'binary_big_endian'
   fclose(fid);
   fid = fopen(Path,'a','ieee-be');
   Format = 2;
end

for i = 1:NumElements
   if ~isempty(PropertyNames{i})
   	if ~Format										% write ASCII data
      	for k = 1:ElementCount(i)
         	for j = 1:length(PropertyNames{i})
            	if Type{i}(j) <= 8
               	fprintf(fid,[PrintfTypeChar{Type{i}(j)},' '],Data{i}{j}(k));
            	else
               	fprintf(fid,'%u%s ',length(Data{i}{j}{k}),sprintf([' ',PrintfTypeChar{Type{i}(j)-9}],Data{i}{j}{k}));
					end
            end
            
         	fprintf(fid,'\n');
      	end
   	else												% write binary data
      	if all(Type{i} <= 8) & all(Type{i} == Type{i}(1))
         	% property data without list types (fast)
         	tmp = zeros(length(PropertyNames{i}),ElementCount(i));
         
         	for j = 1:length(PropertyNames{i})
            	tmp(j,:) = Data{i}{j}(:)';
         	end
         
         	fwrite(fid,tmp,FWriteTypeNames{Type{i}(j)});
     		elseif all(Type{i} > 8)
      		% only list types
         	Type{i} = Type{i} - 9;
            
         	if length(PropertyNames{i}) == 1
         		% only one list property
            	tmp = FWriteTypeNames{Type{i}(1)};
            
            	for k = 1:ElementCount(i)
               	fwrite(fid,length(Data{i}{1}{k}),'uchar');
               	fwrite(fid,Data{i}{1}{k},tmp);
            	end
         	else
         		% multiple list properties
	         	for k = 1:ElementCount(i)
   					for j = 1:length(PropertyNames{i})
      					fwrite(fid,length(Data{i}{j}{k}),'uchar');
                  	fwrite(fid,Data{i}{j}{k},FWriteTypeNames{Type{i}(j)});
						end
            	end
         	end
      	else
      		% mixed type
		      for k = 1:ElementCount(i)
   		      for j = 1:length(PropertyNames{i})
      		      if Type{i}(j) <= 8
         		      fwrite(fid,Data{i}{j}(k),FWriteTypeNames{Type{i}(j)});
	         	   else
   	         	   fwrite(fid,length(Data{i}{j}{k}),'uchar');
                     fwrite(fid,Data{i}{j}{k},FWriteTypeNames{Type{i}(j)-9});
                  end
					end
	         end
         end
      end
   end
end

fclose(fid);



function [area, volume, length, rav, rav2, g] = beadParams(r,A,g,SAnorm)
%function to compute parameters from a beaded unduloid.  In particular, it
%compuates the modified values if surface area is to be normalized.

if isscalar(r) && isscalar(A)
   dotype = 'analytic';
    if nargin < 3
        g = 0;
    end
    if nargin < 4
        SAnorm = 0;
    end
else 
    dotype = 'numerical';
    xm = r;
    ym = A;
end



switch dotype

    case 'numerical'
        % numerical solutions
        length = max(xm) - min(xm);
        
        %from Weisstein, Eric W. "Surface of Revolution." From MathWorld--A
        %Wolfram Web Resource. http://mathworld.wolfram.com/SurfaceofRevolution.html 
        area = 0;
        volume = 0;
        for np=1:size(ym,2)-1
            h = (xm(np+1) - xm(np));
            R1 = ym(np+1);
            R2 = ym(np);
            area = area + pi.*(R1+R2).*sqrt((R1-R2).^2 + (h).^2);
            volume = volume + 1/3 * pi * (h) * (R1.^2 + R1*R2 + R2.^2);
        end
        
        rav = (max(abs(ym))+min(abs(ym)))./2;
        rav2 = sqrt(volume/(pi*length));
        
        
    case 'analytic'
    % this is the analytic solution from the paper
    % "Biomechanics of stretch-induced beading" Markin et al, 1999
    % is NaN for A = 0 (cylinder), so substitute 0.01 
    
        [fb, rmax] = unduloidCoords(r,A,0);
        [fb, rmin] = unduloidCoords(r,A,pi);

        rav = (rmax+rmin)/2;
        
        
        q = quadgk(@(x)subfun(x,A),1-A,1+A);
        qq = quadgk(@(x)subfun2(x,A),1-A,1+A);
        qqq = quadgk(@(x)subfun3(x,A),1-A,1+A);

        % this is equal to pi
        smalle = .001;
        qr0 = pi;%quadgk(@(x)subfun(x,smalle),1-smalle,1+smalle);
        qqr0 = pi;%quadgk(@(x)subfun2(x,smalle),1-smalle,1+smalle);
        qqqr0 = quadgk(@(x)subfun3(x,smalle),1-smalle,1+smalle);
        
        % conservation of certain parameters such as surface area, length, etc...
        if SAnorm
            %rav = sqrt((r^2*(qqr0 + g*qr0))./(qq + (1-A)*g*q));
           
            %use the following to try and keep length constant - testing
            %since this will then change rav...
            %g = (r*qr0 + g*r*qr0 - rav*q)./(r*q);
            
            %[rav,g] = fminsearch(@(x)gnorm(x,r,A,g,q,qq,qr0,qqr0),[rav,g]);
            %[rav,fz,e] = fzero(@(x)gnorm(x,r,A,g,q,qq,qr0,qqr0),rav);
            %g = ((r.*qr0 + g.*r.*qr0 - rav.*q)./(r.*q)); 
            
            % g is Lc/LorigBead, not beaded length.
            %[rav,fz,e] = fzero(@(x)myfun(x,r,g,A,q,qq,qr0,qqr0),r);
            
            g0 = g;
            %use the following to try and keep length constant  
            x0 = [rav; g];
            options=optimset('Display','off','TolFun',1e-12);
            % Surface Area and Length conserved
            [x,fval,exitflag,output] = fsolve(@(x)sa_l_norm(x,r,A,g,q,qq,qr0,qqr0),x0,options);
            rav = x(1) ;
            g = x(2) ;  
            
            % Surface Area and Volume conserved
            %[x,fval,exitflag,output] = fsolve(@(x)sa_vol_norm(x,r,A,g,qq,qqq,qqr0,qqqr0),x0,options);
            %rav = x(1) ;
            %g = x(2) ; 
            
            %x0 = [g];
            %[x,fval,exitflag,output] = fsolve(@(x)sa_rav_vol_norm(x,r,A,g,qq,qqq,qqr0,qqqr0),x0,options);
            %rav = rav;
            %g =  x(1);
            
            %x0 = [rav,g];
            %[x,fval,exitflag,output] = fsolve(@(x)sa_vol_eq_r0over2(x,r,A,g,qq,qqq,qqr0,qqqr0),x0,options);
            %rav = x(1);
            %g =  x(2);
            
            
            %[rav,fz,e] = fzero(@(x)myfun(x,r,g,A,q,qq,qr0,qqr0),r);
            %[rav,fz,e] = fzero(@(x)myfun2(x,r,A,q,qq,qr0,qqr0),r)
            %g = (rav.^2*qq - r.*rav.*q)./(r.^2.*qr0-r.*rav.*(1-A).*qr0);
        end
        
        %length = 2*(1+g)*rav*(q);
        %area = 4*pi*rav^2*(qq + (1-A)*g*q);
        %volume = 2*pi*rav^3*(qqq + (1-A)^2*g*q);
        %rav2 = sqrt(volume/(pi*length));
        
        %if g is Lc/LorigBead, not beaded length.
        % this is also still in testing mode
        length = 2*rav*(q) + g*2*r*qr0;
        area = 4*pi*rav^2*(qq) + 4*pi*r*rav*(1-A)*g*qr0;
        volume = 2*pi*rav^3*(qqq) + 2*pi*rav^2*r*(1-A)^2*g*qr0;
        rav2 = sqrt(volume/(pi*length));
    
    
    otherwise
        error('arguments should be all scalars (analytic) or matrices of equal size (numerical)');
end

function f = subfun(x,A)
f = (x.^2 + 1 - A.^2)./sqrt(4.*x.^2 - (x.^2 + 1 - A.^2).^2);

function ff = subfun2(x,A)
ff = x.*sqrt(1 + ((x.^2 + 1 - A.^2)./sqrt(4.*x.^2 - (x.^2 + 1 - A.^2).^2)).^2);

function fff = subfun3(x,A)
fff = x.^2.*(x.^2 + 1 - A.^2)./sqrt(4.*x.^2 - (x.^2 + 1 - A.^2).^2);

function f4 = myfun(x,r,g,A,q,qq,qr0,qqr0)
f4 = x.^2.*qq + (1-A).*x.*r.*g.*qr0 - qqr0.*r.^2 - qr0.*g.*r.^2;

function f4 = myfun2(x,r,A,q,qq,qr0,qqr0)
g = (rav.^2*qq - r.*rav.*q)./(r.^2.*qr0-r.*rav.*(1-A).*qr0);
f4 = x.^2.*qq + (1-A).*x.*r.*g.*qr0 - qqr0.*r.^2 - qr0.*g.*r.^2;
           
function f = sa_l_norm(x,r,A,g0,q,qq,qr0,qqr0)  
% Surface area and total length held constant, Rav and g free parameters
f = [x(1).^2.*qq + (1-A).*x(1).*r.*x(2).*qr0 - qqr0.*r.^2 - qr0.*g0.*r.^2;  
     ((r.*qr0 + g0.*r.*qr0 - x(1).*q)./(r.*qr0)) - x(2)]; 
 
 function f = sa_vol_norm(x,r,A,g0,qq,qqq,qqr0,qqqr0)  
% Surface area and volume held constant, Rav and g free parameters
f = [x(1).^2.*qq + (1-A).*x(1).*r.*x(2).*qqr0 - qqr0.*r.^2 - qqr0.*g0.*r.^2;  
     x(1).^3.*qqq + x(1).^2.*(1-A).^2.*qqr0.*r.*x(2) - qqqr0.*r.^3 - r.^3.*qqqr0.*g0]; 
 
 function f = sa_rav_vol_norm(x,r,A,g0,qq,qqq,qqr0,qqqr0)  
% Surface area volume and rav held constant, g free parameter
f = [r.^2.*qq + (1-A).*r.*r.*x(1).*qqr0 - qqr0.*r.^2 - qqr0.*g0.*r.^2;  
     r.^3.*qqq + r.^2.*(1-A).^2.*qqr0.*r.*x(1) - qqqr0.*r.^3 - r.^3.*qqqr0.*g0]; 

