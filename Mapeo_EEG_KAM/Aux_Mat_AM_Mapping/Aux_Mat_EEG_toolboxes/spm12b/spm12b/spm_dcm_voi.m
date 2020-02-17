function spm_dcm_voi(DCM_filename,voi_filenames)
% Insert new regions into a DCM model
% FORMAT spm_dcm_voi (DCM_filename,voi_filenames)
%
% DCM_filename   - Name of DCM file
% voi_filenames  - Cell array of new VOI filenames 
%                  eg. {'VOI_V1','VOI_V5','VOI_PPC'}
%
% The RT is assumed to be the same as before.
%
% This function can be used, for example, to replace subject X's data by 
% subject Y's. The model can then be re-estimated without having to go 
% through model specification again.
%__________________________________________________________________________
% Copyright (C) 2002-2012 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_dcm_voi.m 5092 2012-12-03 16:35:00Z guillaume $


load(DCM_filename);

% Check we have matching number of regions
%--------------------------------------------------------------------------
n = length(voi_filenames);
if n ~= DCM.n
    error('Mismatching number of regions');
end

% Replace relevant fields in DCM with xY
%--------------------------------------------------------------------------
DCM               = rmfield(DCM,'xY');
DCM.Y.y           = [];
for i=1:n
    load(voi_filenames{i});
    
    DCM.v         = size(xY.u,1);
    DCM.Y.y(:,i)  = xY.u;
    DCM.Y.name{i} = xY.name;
    DCM.Y.X0      = xY.X0;
    DCM.Y.Q       = spm_Ce(ones(1,DCM.n)*DCM.v);
    DCM.xY(i)     = xY;
end

% Save (overwrite) new DCM file
%--------------------------------------------------------------------------
save(DCM_filename, 'DCM', spm_get_defaults('mat.format'));
