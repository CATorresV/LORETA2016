function varargout = spm_atlas(action,varargin)
% Atlas multi-function
% FORMAT xA = spm_atlas('load',atlas)
% FORMAT L = spm_atlas('list',{'installed','available'})
% FORMAT [S,sts] = spm_atlas('select',xA)
% FORMAT Q = spm_atlas('query',A,XYZmm)
% FORMAT spm_atlas('label',xA)
% FORMAT VM = spm_atlas('mask',xA,label)
% FORMAT V = spm_atlas('prob',xA,label)
% FORMAT V = spm_atlas('maxprob',xA,thresh)
% FORMAT sts = spm_atlas('install',A)
%
% FORMAT hC = spm_atlas('menu',F)
% FORMAT D = spm_atlas('dir')
% FORMAT def = spm_atlas('def')
%__________________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_atlas.m 5699 2013-10-17 12:21:28Z guillaume $


fprintf('**** Please do not use spm_atlas yet as syntax WILL change ****\n');

if ~nargin, action = 'load'; end

%==========================================================================
switch lower(action), case 'dir'
%==========================================================================
    % FORMAT D = spm_atlas('dir')
    %-Return directory containing atlas files

    d = fullfile(spm('Dir'),'atlas');
    
    varargout = { d };
    

%==========================================================================
case 'def'
%==========================================================================
    % FORMAT def = spm_atlas('def')
    %-Return link to atlas definition file
    
    def = fullfile(spm_atlas('Dir'),'atlas.xml');
    % def = 'http://www.fil.ion.ucl.ac.uk/spm/ext/atlas.xml'
    
    varargout = { def };
    
    
%==========================================================================
case 'load'
%==========================================================================
    % FORMAT xA = spm_atlas('load',atlas)
    %-Load atlas

    if isempty(varargin) || isempty(varargin{1})
        [atlas,sts] = spm_atlas('select');
        if ~sts, varargout = {[]}; return; end
    else
        atlas = varargin{1};
    end
    
    if isstruct(atlas), varargout = { atlas }; return; end
    
    %-Read Description
    switch spm_file(atlas,'ext')
        case 'xml'
            descfile  = atlas;
            xA.X      = convert(xmltree(descfile));
            atlasfile = fullfile(spm_file(descfile,'path'),xA.X.atlasfile);
            labelfile = fullfile(spm_file(descfile,'path'),xA.X.labelfile);
        case {'nii','img'}
            atlasfile = atlas;
            labelfile = spm_file(atlas,'ext','txt');
            descfile  = spm_file(atlas,'ext','xml');
            if spm_existfile(descfile)
                xA.X  = convert(xmltree(descfile));
            else
                xA.X  = struct([]);
            end
        otherwise
            list      = spm_atlas('List','installed');
            idx       = find(ismember({list.name},atlas));
            if ~isempty(idx)
                xA.X  = list(idx).info;
                atlasfile = fullfile(spm_file(list(idx).file,'path'),xA.X.atlasfile);
                labelfile = fullfile(spm_file(list(idx).file,'path'),xA.X.labelfile);
            else
                error('Unknown atlas "%s".',atlas);
            end
    end
    
    %-Read Labels
    xA.L      = read_labels(labelfile);
    
    %-Read Atlas
    xA.V      = spm_vol(atlasfile); % remove frame number for 4D files?

    varargout = { xA };
    
    
%==========================================================================
case 'list'
%==========================================================================
    % FORMAT L = spm_atlas('list',{'installed','available'})
    %-Return list of installed or available atlases
    
    if isempty(varargin), varargin = {'installed'}; end
    
    switch lower(varargin{1})
        case 'installed'
            L = atlas_list_installed(varargin{2:end});
        case 'available'
            atlaslist = spm_atlas('def');
            if ismember(atlaslist(1:4),{'http','file'})
                [s, sts] = urlread(atlaslist);
                if ~sts, error('Cannot access "%s".',atlaslist); end
                atlaslist = s;
            elseif ~spm_existfile(atlaslist)
                error('Cannot open "%s".',atlaslist);
            end
            L = convert(xmltree(atlaslist));
        otherwise
            error('Unknown option.');
    end
    
    varargout = { L };
    
    
%==========================================================================
case 'label'
%==========================================================================
    % FORMAT spm_atlas('label',atlas)
    %-Use atlas to label suprathreshold features
    
    spm('Pointer','Watch')

    xA = spm_atlas('load',varargin{:});

%     F  = spm_figure('GetWin','Satellite');
%     spm_figure('Focus',F);
%     spm_results_ui('Clear',F);
%     
%     %-Display activation labels
%     %----------------------------------------------------------------------
%     FS    = spm('FontSizes');
%     PF    = spm_platform('fonts');
%     
%     hAx   = axes('Parent',F,...
%                  'Position',[0.025 0.05 0.95 0.9],...
%                  'DefaultTextFontSize',FS(8),...
%                  'DefaultTextInterpreter','Tex',...
%                  'DefaultTextVerticalAlignment','Baseline',...
%                  'Tag','XXXXXXXXXXXXXXX',...
%                  'Units','points',...
%                  'Visible','off');
% 
%     AxPos = get(hAx,'Position'); set(hAx,'YLim',[0,AxPos(4)])
%     dy    = FS(9);
%     y     = floor(AxPos(4)) - dy;
% 
%     text(0,y,['Atlas:  \it\fontsize{',num2str(FS(9)),'}',xA.X.name],...
%               'FontSize',FS(11),'FontWeight','Bold');   y = y - dy/2;
%     line([0 1],[y y],'LineWidth',3,'Color','r'),        y = y - 9*dy/8;
%     
%     set(hAx,'DefaultTextFontName',PF.helvetica,'DefaultTextFontSize',FS(8))
%     
%     text(0.01,y,'mm mm mm','Fontsize',FS(8));
%     text(0.15,y,'label','Fontsize',FS(8));
% 
%     y     = y - dy/2;
%     line([0 1],[y y],'LineWidth',1,'Color','r')
%     y     = y - dy;
%     y0    = y;
%     
%     TabDat = evalin('base','TabDat');
%     
%     for i=1:size(TabDat.dat,1)
%         XYZmm = TabDat.dat{i,12};
%         if  ~isempty(TabDat.dat{i,5}), fw = 'Bold'; else fw = 'Normal'; end
%         h = text(0.01,y,sprintf(TabDat.fmt{12},XYZmm),...
%                      'FontWeight',fw);
%         lab = spm_atlas('query',xA,XYZmm);
%         h = text(0.1,y,strrep(lab,'_','\_'),'FontWeight',fw);
%         y = y - dy;
%     end
    
    %-Add contextual menus to coordinates
    %----------------------------------------------------------------------
    hAx       = findobj('Tag','SPMList');
    if isempty(hAx), spm('Pointer','Arrow'), return; end
    UD        = get(hAx,'UserData');
    HlistXYZ  = UD.HlistXYZ(ishandle(UD.HlistXYZ));

    for i=1:numel(HlistXYZ)
        h     = uicontextmenu('Parent',ancestor(hAx,'figure'));
        XYZmm = get(HlistXYZ(i),'UserData');

        %-Consider peak only
        %------------------------------------------------------------------
        labk  = spm_atlas('query',xA,XYZmm);
        
        hi    = uimenu(h,'Label',['<html><b>' labk '</b></html>']);

        %-Consider a 10mm sphere around the peak
        %------------------------------------------------------------------
        [labk,P] = spm_atlas('query',xA,...
            struct('def','sphere','spec',10,'xyz',XYZmm));
        
        for j=1:numel(labk)
            hj   = uimenu(hi,'Label',sprintf('<html><b>%s</b> (%.1f%%)</html>',labk{j},P(j)),...
                'Callback',['web(''' spm_atlas('weblink',XYZmm) ''',''-notoolbar'');']);
        end

        set(HlistXYZ(i),'UIContextMenu',h);
    end
    
    spm('Pointer','Arrow')


%==========================================================================
case 'menu'
%==========================================================================
    % FORMAT hC = spm_atlas('menu',F)
    %-Create user interface atlas menu
    
    if nargin < 2, Finter = 'Interactive'; else Finter = varargin{1}; end
    Finter = spm_figure('GetWin',Finter);
    
    %hC  = uicontextmenu;
    hC   = uimenu(Finter,'Label','Atlas', 'Tag','AtlasUI');
    
    hC1  = uimenu(hC,'Label','Use Atlas');
    
    list = spm_atlas('List','installed');
    for i=1:numel(list)
        uimenu(hC1,'Label',list(i).name,...
            'Callback',sprintf('spm_atlas(''label'',''%s'');',list(i).name));
    end
    if isempty(list), set(hC1,'Enable','off'); end
    
    
    hC2  = uimenu(hC,'Label','Download Atlas...',...
        'Separator','on',...
        'Callback','spm_atlas(''install'');');
    
    %set(Finter,'uicontextmenu',hC);
    
    varargout = { hC };
    
    
%==========================================================================
case 'select'
%==========================================================================
    % FORMAT [S,sts] = spm_atlas('select',xA)
    %-Select atlas or labels
    
    S = '';
    if isempty(varargin)
        [S,sts] = spm_select(1, {'image','xml'},...
            'Select Atlas...', {}, spm_atlas('Dir'));
        if ~sts, varargout = { S, sts }; return; end
    else
        xA = spm_atlas('load',varargin{1});
        [sel,sts] = listdlg(...
            'ListString',xA.L{2},...
            'SelectionMode','multiple',...
            'ListSize', [400 300],...
            'Name','Select label(s)',...
            'PromptString',sprintf('Labels from %s atlas:',xA.X.name));
        if ~sts, varargout = { S, sts }; return; end
        S  = xA.L{2}(sel);
    end
    varargout = { S, sts };
    
    
%==========================================================================
case 'query'
%==========================================================================
    % FORMAT Q = spm_atlas('query',xA,XYZmm)
    % FORMAT [Q,P] = spm_atlas('query',xA,xY)
    %-Atlas query
    
    xA = spm_atlas('load',varargin{1});
    if nargin > 2, xY = varargin{2}; else xY = struct; end
    
    unknown = '????';
    
    if numel(xA.V) == 1 % or xA.X contains type definition
        if isnumeric(xY)
            %-peak
            XYZmm     = xY;
            XYZ       = xA.V.mat\[XYZmm;1];
            vpeak     = spm_sample_vol(xA.V,XYZ(1),XYZ(2),XYZ(3),0);
            j         = xA.L{3} == vpeak;
            if any(j) == 1
                Q     = xA.L{2}{j};
            else
                Q     = unknown;
            end
            
            varargout = { Q };
        else
            %-cluster
            v         = spm_summarise(xA.V,xY);
            vu        = unique(v);
            vun       = histc(v,vu);
            [vun,is]  = sort(vun(:),1,'descend');
            vu        = vu(is);
            for j=1:numel(vu)
                k     = xA.L{3} == vu(j);
                if any(k) == 1
                    Q{j} = xA.L{2}{k};
                else
                    Q{j} = unknown;
                end
                P(j)  = 100*vun(j)/numel(v);
            end
            
            varargout = { Q, P };
        end
    else
        P             = xA.L{2};
        if isnumeric(xY)
            %-peak
            XYZmm     = xY;
            XYZ       = xA.V(1).mat\[XYZmm;1];
            Q         = spm_get_data(xA.V,XYZ); % which interp to use?
            
            varargout = { Q, P };
        else
            %-cluster
            v         = spm_summarise(xA.V,xY);
            v         = mean(v,2);
            
            varargout = { v, P };
        end
    end
    
     
%==========================================================================
case 'mask'
%==========================================================================
    % FORMAT VM = spm_atlas('mask',xA,label)
    %-Return binary mask for given labels

    if nargin < 2, xA = ''; else xA = varargin{1}; end
    xA = spm_atlas('load',xA);
    if nargin < 3, label = spm_atlas('select',xA);
    else label = varargin{2}; end
    label = filter_labels(xA,label);
    
    if numel(xA.V) == 1 % or xA.X contains type definition
        VM = struct(...
            'fname',   [xA.X.name '_mask' spm_file_ext],...
            'dim',     xA.V(1).dim,...
            'dt',      [spm_type('uint8') spm_platform('bigend')],...
            'mat',     xA.V(1).mat,...
            'n',       1,...
            'pinfo',   [1 0 0]',...
            'descrip', sprintf('%s mask',xA.X.name));
        VM.dat = false(VM.dim);
        
        D = spm_read_vols(xA.V);
        for i=1:numel(label)
            j = find(ismember(xA.L{2},label{i}));
            for k=1:numel(j)
                VM.dat = VM.dat | (D == xA.L{3}(j(k)));
            end
        end
        VM.dat = uint8(VM.dat);
    else
        if nargin < 4, thresh = 0.5; else thresh = varargin{3}; end
        VM       = spm_atlas('prob',xA,label);
        VM.dt(1) = spm_type('uint8');
        VM.dat   = uint8(VM.dat > thresh);
    end
    
    varargout = { VM };
    % The output mask can be saved to disk with:
    % VM = spm_write_vol(VM,VM.dat);
    % VM = rmfield(VM,'dat');
    
	
%==========================================================================
case 'maxprob'
%==========================================================================
    % FORMAT V = spm_atlas('maxprob',xA,thresh)
    
    if nargin < 2, xA = ''; else xA = varargin{1}; end
    xA = spm_atlas('load',xA);
    if nargin < 3, thresh = 0; else thresh = varargin{2}; end
    
    V = struct(...
        'fname',   [xA.X.name '_maxprob_thresh' num2str(thresh) spm_file_ext],...
        'dim',     xA.V(1).dim,...
        'dt',      [spm_type('uint8') spm_platform('bigend')],...
        'mat',     xA.V(1).mat,...
        'n',       1,...
        'pinfo',   [1 0 0]',...
        'descrip', sprintf('%s mask',xA.X.name));
    V.dat = zeros(V.dim);
    
    for i=1:V.dim(3)
        Y = zeros(V.dim(1),V.dim(2),numel(xA.V));
        for j=1:numel(xA.V)
            Y(:,:,j) = spm_slice_vol(xA.V(j),spm_matrix([0 0 i]),V.dim(1:2),0);
        end
        [Y,V.dat(:,:,i)] = max(Y,[],3);
        V.dat(:,:,i) = V.dat(:,:,i) .* (Y > thresh);
    end
    V.dat = uint8(V.dat);
    
    varargout = { V };
    
    
%==========================================================================
case 'prob'
%==========================================================================
    % FORMAT V = spm_atlas('prob',xA,label)
    
    if nargin < 2, xA = ''; else xA = varargin{1}; end
    xA = spm_atlas('load',xA);
    if nargin < 3, label = spm_atlas('select',xA);
    else label = varargin{2}; end
    [label,idx] = filter_labels(xA,label);
    
    if numel(idx) == 1
        descrip = label{1};
    else
        descrip = sprintf('%s prob',xA.X.name);
    end
    V = struct(...
        'fname',   [xA.X.name '_prob' spm_file_ext],...
        'dim',     xA.V(1).dim,...
        'dt',      [spm_type('float32') spm_platform('bigend')],...
        'mat',     xA.V(1).mat,...
        'n',       1,...
        'pinfo',   [1 0 0]',...
        'descrip', descrip);
    V.dat = zeros(V.dim);
    
    for i=1:numel(idx)
        V.dat = V.dat + spm_read_vols(xA.V(idx(i)));
    end
    V.dat = single(V.dat);
    
    varargout = { V };
    
    
%==========================================================================
case 'install'
%==========================================================================
    % FORMAT sts = spm_atlas('install',A)
    %-Install Atlas
    
    if isempty(varargin)
        h = atlas_figure('SPM Atlases');
        H = getappdata(h,'Hbrowser');
        
        if ispc, s = '/'; else s = ''; end
        localurl = @(f) sprintf(['file://' s strrep(spm('Dir'),'\','/') '/help/' f]);
        
        tpl = spm_file_template(fullfile(spm('Dir'),'help'),'keep');
        tpl = tpl.file('TPL_ATLAS','spm_atlas.tpl');
        tpl = tpl.block('TPL_ATLAS','loading','load');
        tpl = tpl.block('TPL_ATLAS','listing','list');
        tpl = tpl.block('listing','atlas','atl');
        tpl = tpl.var('SPM',spm('Ver'));
        tpl = tpl.var('SPM_CSS',localurl('spm.css'));

        tpl = tpl.var('IMG_LOADING',localurl('images/loading.gif'));
        tpl = tpl.var('TXT_LOADING','Loading description of atlases...');
        tpl = tpl.parse('load','loading',0);
        tpl = tpl.var('list','');
        tpl = tpl.parse('OUT','TPL_ATLAS');
        html = get(tpl,'OUT');
        
        spm_browser(html,H);
        
        try
            L = spm_atlas('List','Available');
        catch
            tpl = tpl.var('TXT_LOADING','Cannot access atlases description.');
            tpl = tpl.var('IMG_LOADING','');
            tpl = tpl.parse('load','loading',0);
            tpl = tpl.parse('OUT','TPL_ATLAS');
            html = get(tpl,'OUT');
            spm_browser(html,H);
            varargout = { false };
            return;
        end
        
        tpl = tpl.var('load','');
        tpl = tpl.var('tbx','');
        for i=1:numel(L.atlas)
            tpl = tpl.var('ATLAS_URL',L.atlas{i}.website);
            tpl = tpl.var('ATLAS_NAME',L.atlas{i}.name);
            tpl = tpl.var('ATLAS_ID',L.atlas{i}.name);
            if isstruct(L.atlas{i}.maintainer)
                L.atlas{i}.maintainer = {L.atlas{i}.maintainer};
            end
            for j=1:numel(L.atlas{i}.maintainer)
                tpl = tpl.var('ATLAS_AUTHOR',L.atlas{i}.maintainer{j}.name);
                tpl = tpl.var('ATLAS_EMAIL',L.atlas{i}.maintainer{j}.email);
            end
            tpl = tpl.var('ATLAS_SUMMARY',L.atlas{i}.description);
            tpl = tpl.parse('atl','atlas',1);
        end
        tpl = tpl.parse('list','listing',0);
        tpl = tpl.parse('OUT','TPL_ATLAS');
        html = get(tpl,'OUT');
        
        spm_browser(html,H);
        
        varargout = { true };
        return;
        
    else
        A  = varargin{1};
        
        if ~(nargin == 3 && strcmpi(varargin{2},'-force'))
            AI = spm_atlas('list','installed');
            if ismember(A,{AI.name})
                warning('Atlas "%s" is already installed.',A);
                varargout = { true };
                return;
            end
        end
        
        AA = spm_atlas('list','available');
        for i=1:numel(AA.atlas)
            if strcmp(A,AA.atlas{i}.name)
                url = AA.atlas{i}.download;
                if isempty(url)
                    if desktop('-inuse')
                        str = sprintf('<a href="%s">%s</a>',url,url);
                    else
                        str = url;
                    end
                    fprintf(['Atlas "%s" is only available through website:\n' ...
                        '  ' str '\n'], A);
                    if ~spm('CmdLine'), web(url,'-browser'); end
                    varargout = { false };
                    return;
                else
                    %-Check folder permissions
                    %------------------------------------------------------
                    dest = spm_atlas('Dir');
                    [sts, attrb] = fileattrib(dest);
                    if ~sts, error('"%s"\n%s',dest,attrb); end
                    if ~attrb.UserWrite
                        error('No write access to "%s".\nMaybe use "%s" instead.',...
                            dest, strrep(userpath,pathsep,''));
                    end
                    
                    %-Download atlas archive
                    %------------------------------------------------------
                    tmpfile = [tempname(dest) '.zip'];
                    try
                        F = urlwrite(url,tmpfile);
                    catch
                        l = lasterror;
                        switch l.identifier
                            case 'MATLAB:urlwrite:ConnectionFailed'
                                error('Could not access URL "%s".',url);
                            case 'MATLAB:urlwrite:InvalidOutputLocation'
                                error('Could not create output file "%s".',tmpfile);
                            otherwise
                                rethrow(l);
                        end
                    end
                    
                    %-Unzip archive in destination folder
                    %------------------------------------------------------
                    try
                        FS = unzip(F,dest);
                    catch
                        spm_unlink(F);
                        error('Error when unpackig atlas archive');
                    end
                    
                    %-Delete atlas archive
                    %------------------------------------------------------
                    spm_unlink(F);
                    
                    %-Display <ATLAS>_README.txt if present
                    %------------------------------------------------------
                    fprintf('Atlas "%s" installed.',A);
                    
                    idx = find(~cellfun('isempty',...
                        regexpi(FS,sprintf('%s_README.txt$',A))));
                    if ~isempty(idx)
                        type(FS{idx(1)});
                    end
                    
                    % Refresh list of installed atlases.
                    spm_atlas('list','installed','-refresh');
                    
                    varargout = { true };
                    return;
                end
            end
        end
        
        warning('Cannot find atlas "%s".',A);
        varargout = { false };
        return;
    end
    
    
%==========================================================================
case 'weblink'
%==========================================================================
    % FORMAT url = spm_atlas('weblink',XYZmm,website)
    %-Return URL for coordinates query
    
    XYZmm = varargin{1};
    if nargin < 3, website = 'Brede'; end
    
    switch lower(website)
        case 'brede'
            %-Brede Database - Talairach coordinate search
            url = 'http://neuro.imm.dtu.dk/cgi-bin/brede_loc_query.pl?q=%d+%d+%d';
        case 'neurosynth'
            %-Neurosynth - structure-to-function mappings
            url = 'http://neurosynth.org/locations/%d_%d_%d';
        otherwise
            error('Unknown website "%s".',website);
    end
    
    url = sprintf(url,XYZmm);
    
    varargout = { url };
    
    
%==========================================================================
otherwise
%==========================================================================
    error('Unknown action.');
end


%==========================================================================
% FUNCTION labels = read_labels(labelfile)
%==========================================================================
function labels = read_labels(labelfile)
fid        = fopen(labelfile,'rt');
if fid == -1, error('Cannot find atlas labels: "%s".',labelfile); end
try
    % AAL
    % ShortName Long_Name Key
    L = textscan(fid,'%s %s %d');
    labels = L;
    
    % FreeSurfer
    % Key Long-Name Red Green Blue Alpha
    % L = textscan(fid,'%d %s %d %d %d %d','CommentStyle','#');
    % labels = {{} L{2} L{1}};
    
    % Brainvisa
    % Key, X, Y, Z, Red, Green, Blue, Acronym, Short Name
    % L = textscan(fid,'%d %f %f %f %d %d %d %s %s','Delimiter',',','HeaderLines',1);
    % labels = {L{8} L{9} L{1}};
    
    % Talairach
    % Key <TAB> Long Name
    % L = textscan(fid,'%d %s','Delimiter','\t');
    % labels = {{} L{2} L{1}};
    
    % MRIcron JHU-WhiteMatter
    % Key <TAB> Long_Name
    % L = textscan(fid,'%d %s','Delimiter','\t');
    % labels = {{} L{2} L{1}};
    
    % MRIcron Brodmann
    % labels = {{} cellfun(@(x) sprintf('Brodmann %d',x),num2cell(1:52),'UniformOutput',0) (1:52)};
    
    % WFU PickAtlas
    % Key <TAB> Long Name <TAB> *
    % L = textscan(fid,'%d %s %*[^\n]','Delimiter','\t','CommentStyle',{'[' ']'});
    
    % Hammers_mith
    % Long_Name *
    % L = textscan(fid,'%s','Delimiter',',','HeaderLines',2,'MultipleDelimsAsOne',1);
    % labels = {{} L{1}(1:end-2) (1:numel(L{1})-2)};
    
    % Joern's Cerebellum MNIsegment-MRICroN
    % Key Long_Name ???
    % L = textscan(fid,'%d %s %d','Delimiter',' ');
    % labels = {{} L{2} L{1}};
    
    % FSL
    % XML: <label index="Key" x="X" y="Y" z="Z">Long Name</label>
    % X = xmltree(labelfile);
    % I = find(X,'/atlas/data/label');
    % for i=1:numel(I)
    %   A = attributes(X,'get',I(i)); A = [A{:}];
    %   L{3}(i) = str2double(A(strcmp({A.key},'index')).val); % + 1 ?
    %   L{2}{i} = get(X,children(X,I(i)),'value');
    % end
    % labels = L;
    
    % Colin 27: ITK-SNAP Label Description File
    % Key Red Green Blue Alpha Vis Idx "Long Name"
    % L = textscan(fid,'%d %d %d %d %d %d %d "%[^"]"','CommentStyle','#');
    % labels = {{} L{8} L{1}};
catch
    fclose(fid);
    error('Cannot read atlas labels in: "%s".',labelfile);
end
fclose(fid);


%==========================================================================
% FUNCTION save_labels(labelfile,labels)
%==========================================================================
function save_labels(labelfile,labels)
fid = fopen(labelfile,'wt');
if fid == -1, error('Cannot write file "%s".',labelfile); end
for i=1:numel(labels{3})
    fprintf(fid,'%d\t%s\n',labels{3}(i),labels{2}{i});
end
fclose(fid);


%==========================================================================
% FUNCTION [labels,i] = filter_labels(xA,labels)
%==========================================================================
function [labels,i] = filter_labels(xA,labels)
if ~iscellstr(labels)
    labels = xA.L{2}(~cellfun(@isempty,regexp(xA.L{2},labels)));
end
[labels,i] = intersect(xA.L{2},labels);


%==========================================================================
% FUNCTION L = atlas_list_installed(refresh)
%==========================================================================
function L = atlas_list_installed(refresh)
persistent atlas_list
if nargin && strcmpi(refresh,'-refresh'), atlas_list = []; end
if isempty(atlas_list)
    L = spm_select('FPList',spm_atlas('Dir'),'^.*\.xml$');
    if isempty(L), L = {}; else L = cellstr(L); end
    atlas_list = struct('file',{},'info',{},'name',{});
    for i=1:numel(L)
        try
            A.file = L{i};
            A.info = convert(xmltree(L{i}));
            A.name = A.info.name;
            atlas_list(end+1) = A;
        catch
            fprintf('Atlas "%s" could not be read.\n',L{i});
        end
    end
end
L = atlas_list;


%==========================================================================
% FUNCTION atlas_figure
%==========================================================================
function h = atlas_figure(name)

if ~nargin, name = 'Atlases'; end
h = spm_figure('FindWin','SPMatlas');
if ~isempty(h), set(h,'Name',name); return; end

h = figure(...
    'MenuBar',     'none',...
    'NumberTitle', 'off',...
    'Name',        name,...
    'Resize',      'off',...
    'Toolbar',     'none',...
    'Tag',         'SPMatlas',...
    'WindowStyle', 'Normal',... %'Modal'
    'Color',       [1 1 1],...
    'Visible',     'off');
pos = get(h,'Position');
pos([3 4]) = [350 400];
set(h,'Position',pos);

[H, HC] = spm_browser('<html></html>',h,[2 2 pos(3)-4 pos(4)-4],'html');
setappdata(h,'Hbrowser',H); setappdata(h,'HCbrowser',HC);
set(h,'Resize','on','ResizeFcn',@atlas_figure_resize);

set(h,'Visible','on');


%==========================================================================
% FUNCTION atlas_resize(obj,evt,varargin)
%==========================================================================
function atlas_figure_resize(obj,evt,varargin)
old_units = get(obj,'Units');
set(obj, 'Units','pixels');
figpos    = get(obj,'Position');
H         = getappdata(obj,'HCbrowser');
set(H,   'pos',[2 2 figpos(3)-4 figpos(4)-4]);
set(obj, 'Units',old_units);
