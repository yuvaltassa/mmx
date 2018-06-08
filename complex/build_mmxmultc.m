function build_mmxmultc(verbose)
% build_mmxc - compiles mmxc and mmxc_nothreads
%
%       build_mmxmultc()
%       build_mmxmultc(verbose)
%
% where verbose is a boolean. When it is set to true, all compilation info will be displayed
%
%  build_mmxmultc has been tested on Win64 and Linux 64

if nargin == 0
    verbose = false;
end

% make sure to be in the mmx/complex folder
cd(fileparts(mfilename('fullpath')))

% get the current computer architecture
arch        = computer('arch');

% preallocate the link and define cell arrays
[link, define]  = deal({});
% fill the link and define cell arrays
switch arch
    case {'win64','win32'}
        define   = {'WIN_SYSTEM'};
    case {'glnxa64','glnx86'}
        link     = {'pthread'};
        define   = {'UNIX_SYSTEM'};
    case {'maci64'}
        link     = {'pthread'};
        define   = {'UNIX_SYSTEM'};     
    otherwise
        error unsupported_architecture
end
% add the -l and -D prefix before the includes and defines
prefix   = @(pref,str_array) cellfun(@(x)[pref x],str_array,'UniformOutput',0);
l_link   = prefix('-l',link);
D_define = prefix('-D',define);
% add the verbose flag if verbose is set to true
if verbose
    verb  = {'-v'};
else
    verb  = {};
end
%% compile mmx_complex
% remove previously compiled mex file from the matlab memory
clear('mmxmultc')
% combine everything together to generate the build command
command = {verb{:}, l_link{:}, D_define{:}}; %#ok<*CCAT>
fprintf('==========\nCompiling mmxmultc, using \n');
fprintf('%s, ',command{:})
fprintf('\n')
% compile
mex(command{:}, '-output', 'mmxmultc', 'mmxmultc.cpp');
fprintf('Compilation of mmxmultc succeeded.\n');
%% compile mmx_complex_nothreads
% remove previously compiled mex file from the matlab memory
clear('mmxmultc_nothreads')
% combine everything together to generate the build command
command = {verb{:}, l_link{:}, D_define{:}}; %#ok<*CCAT>
fprintf('==========\nCompiling mmxmultc_nothreads, using \n');
fprintf('%s, ',command{:})
fprintf('\n')
% compile
mex(command{:}, '-output', 'mmxmultc_nothreads', 'mmxmultc_nothreads.cpp');
fprintf('Compilation of mmxmultc_nothreads succeeded.\n');

end



