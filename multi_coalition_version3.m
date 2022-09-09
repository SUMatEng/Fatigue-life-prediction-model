function [coalated]=multi_coalition_version3(defect)
% Combine and dicretize flaws 
%Extract terms from defect population not input = row vectors 
SZ=defect(1,:)'; %Thorsten functions written for column vectors
TH=defect(2,:)';
R=defect(3,:)';
[X,Y] = pol2cart(TH,R);
A_ext=defect(5,:)';
D=(7.1176/1000); %in um
% Step 1: Find flaw interaction.
interaction = flawinteraction(X,Y,SZ);
% Step 2: Combine interacting flaws. 
Xc = X;
Yc = Y;
SZc = SZ;

A_extc=A_ext;

while ~isempty(interaction)
    interaction = flawinteraction(Xc,Yc,SZc);
    [Xc,Yc,SZc,A_extc] = combineflaws(Xc,Yc,SZc,A_extc,interaction);
end
% Step 3: Determine magnification factor
Ymagc=surface_crack_shaft_geo(SZc',D);

% Step 4: Discretize flaws.

[Yfinal, dc ]= discretizeflaws(Xc,Yc,SZc,D,Ymagc');

% Remove nan values from matrices
SZc=SZc(~isnan(SZc));% Crack sizes
Xc=Xc(~isnan(Xc));% cartesion X location
Yc=Yc(~isnan(Yc));% cartesian Y locations
Yfinal=Yfinal(~isnan(Yfinal)); % Stress magnification factors based on location
A_extc=A_extc(~isnan(A_extc)); 

%Prepare for loop and NASGRO
[Angular_location,Radial_location] = cart2pol(Xc,Yc);
coalated=[SZc';Angular_location';Radial_location';Yfinal';A_extc'];


%% Functions
function interaction = flawinteraction(x,y,r)
% FLAWINTERACTION Computes which flaws are interacting.
%   FLAWINTERACTION(x,y,r) Computes which flaw interaction, assuming
%   circular flaws located at position x, y with radius r. The output
%   provides a cell array using the flaw order as given by x, y and r in
%   which the flaw interaction is given.
%
% Thorsten Becker, 2022

% Organise data.
A = [x(:),y(:)];
r = r(:);
interaction = [];
% Number of flaws.
n = numel(x);
% Distance to cicles.
d2c=real(sqrt((ones(n,1)*sum((A.^2)',1))'+ones(n,1)*sum((A.^2)',1)-2.*(A*(A'))));
% Combined radius.
r2c = repmat(r,1,n)+repmat(r',n,1);
% Find interaction.
[i,j] = ind2sub([n,n],find(triu(d2c-r2c,1)<0));
% Assemble output.
%intflaws = sort(sort([i,j],2),1);
intflaws = unique(sort([i,j],2),'rows','sorted');
% Group multiple interacting flaws.
count = 1;
for ii = 1:n
    if ismember(ii,intflaws)
        flawgroups = ii;
        % Iteratively group interacting flaws.
        while any(ismember(intflaws(:),flawgroups))
            % Find and group interacting flaws.
            rows = any(ismember(intflaws,flawgroups),2);
            flawgroups = unique([flawgroups,reshape(intflaws(rows,:),1,[])]);
            % Remove considered flaws to avoid duplicate grouping.
            intflaws(rows,:) = nan;
        end
        interaction{count,1} = flawgroups;
        count = count+1;
    end
end
end

function [xc,yc,rc,a_extc] = combineflaws(x,y,r,a_ext,interaction)
% COMBINEFLAWS Combines flaws using a weighted average.
%	COMBINEFLAWS(x,y,r,interaction) Combines flaws using interaction.
%   Interaction is computed using FLAWINTERACTION, which provides a cell
%   array of interacting flaws. Flaw positions x and y are combined using a
%   weighted average. The flaw radius is combined using a total equivalent
%   area approach so that the combined area is equal to the sum of all
%   interacting flaws. Interacting flaws are replaced by nan.
%
% Thorsten Becker, 2022

% Organise data.
xc = x;
yc = y;
rc = r;
a_extc=a_ext;
if isempty(interaction), return, end
% Number of flaws.
n = numel(xc);
% Loop through flaws.
for ii = 1:n
    % Check if flaw exists.
    if ~isnan(x(ii))
        % Check if interaction exists.
        idx = find(~cellfun(@isempty,cellfun(@(c) find(c == ii),interaction,'uniform',false))); %Have a question about this
        % If interaction exists.
        if ~isempty(idx)
            % Compute weighted position.
            xc(ii,1) = sum(x(interaction{idx}).*r(interaction{idx}))/sum(r(interaction{idx}));
            yc(ii,1) = sum(y(interaction{idx}).*r(interaction{idx}))/sum(r(interaction{idx}));
            % Compute combined radius
            rc(ii,1) = sqrt(sum(r(interaction{idx}).^2));
            %rc(ii,1) = sum(r(interaction{idx})); %I think it should be this
            
            % If combined the crack extension kept as the largest of the
            % Sets extension to minimum of coalescing cracks wo for crack closure model (Nic)
            % Let us check this later
            a_extc(ii,1)=sum(a_ext(interaction{idx}));
            % Remove remining flaws
            x(interaction{idx}(2:end)) = nan;
            y(interaction{idx}(2:end)) = nan;
            r(interaction{idx}(2:end)) = nan;
            % remove from a_extension vector values 
            a_ext(interaction{idx}(2:end))=nan;
         
            
        end
    else
        xc(ii,1) = x(ii);
        yc(ii,1) = y(ii);
        rc(ii,1) = r(ii);
        a_extc(ii,1)=a_ext(ii);
        
    end
end
end

function Y=surface_crack_shaft_geo(a,D)
% a radius of defect, D diameter. Computes magnification factors 
% Note Only valid for a < D/2
siz=a;
%disp(a);
for ii=1:size(siz,2) %Magnification increases -> D/2 then holds
if ~isnan(siz(ii))
   if siz(ii) >D/2 %Holds Y value at max for formulations if is bigger
    siz(ii)=D/2;
   end
    B=(pi/2)*(siz(ii)/D);
G=0.92*(2/pi)*sec(B)*(((tan(B))/B)^(1/2));
Z=1-sin(B);
Y(ii)=G.*(0.752+1.286.*B+0.37*(Z.^3));
else
    Y(ii)=nan;  
end
end
%Nic 2022
end

function [Yalt,d] = discretizeflaws(x,y,r,sr,Ymagc)
% DISCRETIZEFLAWS Groups flaws into surface, near-surface & internal flaws.
%  DISCRETIZEFLAWS(x,y,r,D) Groups flaws into surface, near-surface &
%  internal flaws based on the sample radius sr. When the outer radius of
%  the flaw touches the sample radius, the flaw is grouped as a surface
%  flaw; if the flaw is within 0.8 of the sample radius, the flaw is
%  grouped as a near-surface flaw; lastly, a flaw is within 0.8 of the
%  sample radius is grouped as an internal flaw. Nan values are ignored.
%
% Thorsten Becker, 2022

% Corrects the magnification factors based off
% surface,near-surface,internal

% Defaults.
Dns = 0.8;
% Number of flaws.
n = numel(x);
% Organise data.
[~,Rad] = cart2pol(x,y);
% Combined radial position + radius.
Rr = Rad+r;
% Loop through flaws. 
for ii = 1:n
        if Rr(ii) >= sr
            d{ii,1} = 'surface';
            Yalt(ii,1)=Ymagc(ii,1);
        elseif r(ii)/(sr-Rr(ii)) >= Dns %Proximity rule a/h (Nic)
            d{ii,1} = 'near';
             Yalt(ii,1)=Ymagc(ii,1);
        elseif  r(ii)/(sr-Rr(ii)) < Dns
            d{ii,1} = 'internal';
             Yalt(ii,1)=(0.5/0.65)*Ymagc(ii,1);
        else
            d{ii,1} = 'nan';
            Yalt(ii,1)=nan;
        end
end

end

end