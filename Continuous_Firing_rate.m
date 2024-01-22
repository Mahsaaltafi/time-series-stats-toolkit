function Hfr=Continuous_Firing_rate(Eth, spikeTimes, prec)

% Calculating continuous firing rate of each neuron during a behaviour. The function computes the continuous firing rate
% of all trials of the behaviour, accounting for different trial lengths,
% and generates a plot displaying the average firing rate across trials.
%
% Input:   
%       Eth           :a cell array containig 2d arrays for each neuron where rows correspond to different trials of
%                      a behaviours and columns refer to [start, end] time (at 20kHz) of each trial,
%                      respectibvely.
%       spikeTimes    :a cell array containing spike times (at 20kHz) of each neuron. Eth and and spikeTimes must have the 
%                      the same length.
%                 
%       prec          :an integer, precison of recorded behaviour in ms (prec=500ms in Chen, Altafi et al. 2024)
%
%
% Output:
%       Hfr           :a matrix of firing rate of a behaviour for each
%                      neuron (rows) aligned. Columns refer to the duration
%                      of behaviours where 1st column is the start and 2nd
%                      column is the end of the behaviour.
%=====================================================================================================================
% Function plots a heatmap of normalized firing rates in case it recieves input with
% more than 10 neurons (sorted by the average firing rate)
%
% Usage: Hfr=Continuous_Firing_rate(Eth, spikeTimes, prec)
%
%
%=====================================================================================================================

% to check arguments 
if ~iscell(Eth) || ~iscell(spikeTimes) 
    error('Input variables Eth and spikeTimes must be in cell format');
end
    

if abs(length(Eth)-length(spikeTimes))>0 
    error('Eth and spiket should be given for same number of cells'); 
end




k=1;
H1=[];
for i=1:length(Eth)
    
    inst_firR=[];
    fprintf(['Processing unit # ', num2str(i), ' out of ', num2str(length(Eth)), '\n']);
    if ~isempty(Eth{i})
       for j=1:size(Eth{i},1)
    
           % computes firing rate for each trial of behaviour
           inst_fir=InstantFirRate(spikeTimes{i}, [Eth{i}(j,1), Eth{i}(j,2)], prec);
   
           
           % segmenting the firing rate to 140 bins for alignment 
           edge=linspace(1,length(inst_fir),141);
           vq=interp1([1:length(inst_fir)], inst_fir, edge);
           for ii=1:140
               inst_firR(j,ii)=nanmean([vq(ii),vq(ii+1)]);
    
           end

       end
       h=nanmean(inst_firR,1);
    else
       h=NaN(1,140);
    end

    H1(i,1:140)=h;   
    k=k+i;        
end

n=size(H1,2);
H1=H1(all(~isnan(H1),2),:); 
iMyPhasewt = smoothn([H1,H1,H1,H1], [1 3],'gaussian',2); % smoothing the firing rate with convolving a Gaussian kernel of size 3 and SD 2
H1=iMyPhasewt(:,n+1:2*n);
H1=H1';

% normalizing the firing firing to (0, 1)
m_1 = H1 - repmat(nanmin(H1), size(H1, 1), 1);
m_2 = nanmax(H1) - nanmin(H1);
mm = m_1./repmat(m_2, size(H1, 1), 1);
Hfr =mm';

if length(Eth)>=10
   HH1=Hfr(all(~isnan(Hfr),2),:); %removing cells with NaN firing rates

   % to sort the cells by the average firing rate
   a1=nanmean(HH1');
   HH1 = [HH1, a1'];

   HHH1 = sortrows(HH1, size(HH1, 2)); 
   HHH1 = HHH1(:,1:size(H1, 1));


   % Plot
   imagescnan([1:n], [1:size(HHH1,1)], HHH1,'CDataMapping', 'scaled','NanColor','w') 
   colormap default
   c = colorbar;
   xlim([1 n])
   set(gca, 'XTick', [1, 140])
      set(gca, 'xTickLabels', {'start', 'end'});
      xlabel('Behaviour')
      ylabel(' # Cells ')
      ylabel(c, 'Firing rate (Norm.)')
   pbaspect([1 1 1])    

end
end

function fir=InstantFirRate(spikeTimes,win, prec)

% Computing instantanous firing rate during a window
% win & spikeTimes at 20KHz

p=prec;
srate=20;    % resolution at 20Hz
t=0;        % the time away from win (It could be adjusted to get also firing rates for +/-t sec away from win.)

spikeTimes=(spikeTimes/20000)*(srate);     % changing resolutions 
win=(win/20000)*srate;

kernel_length= p* 10^ (-3)* srate; % kernel window size, 500ms is the preciosion of behavior
                                   
temp=[win(1)-(t*srate): win(2)+(t*srate)];   % computing for +/- 5 sec around wiondow
x_vals=[temp(1)-kernel_length: temp(end)+kernel_length];


y_vals=zeros(length(x_vals),1);
for j=1:length(x_vals)
    count=sum(spikeTimes>=x_vals(j)-(1/2) & spikeTimes<x_vals(j)+(1/2));
    y_vals(j)=count;
end
cond1=spikeTimes >= temp(1);
cond2=spikeTimes <=temp(end);
spikes= spikeTimes(cond1 & cond2);


    
inst_fir=zeros(length(temp),1);
    
    
parfor k=1:length(temp)
        
       sigma=(p/6)*10^(-3)*srate;  % to calculate kernel for the range [mu-3*sigma:mu+3*sigma] 
       kernel= (1/(sigma*(sqrt(2*pi))))*exp(-(x_vals - temp(k)) .^ 2 / (2*(sigma^2)));%computing kernel for each spike time
           
       kernel(find(kernel<=1*10^(-6)))=0; % to make kernel tobe zero below the given threshold
          
       kernel= kernel/sum(kernel);  % make kernel sum to 1
         
       kernelconv= sum(y_vals .* kernel'); % convolving 
      
       inst_fir(k)= srate*kernelconv;
           
end
fir2=smoothn([inst_fir; inst_fir; inst_fir; inst_fir],[5, 1], 'gaussian',2);  % smoothing the firing rate with convolving a Gaussian kernel of size 5 and SD 2
fir=fir2(length(inst_fir)+1:2*length(inst_fir));

end

% Local functions
function Y = smoothn(X,sz,filt,std)


% SMOOTHN Smooth N-D data
%   Y = SMOOTHN(X, SIZE) smooths input data X. The smoothed data is
%       retuirned in Y. SIZE sets the size of the convolution kernel
%       such that LENGTH(SIZE) = NDIMS(X)
%
%   Y = SMOOTHN(X, SIZE, FILTER) Filter can be 'gaussian' or 'box' (default)
%       and determines the convolution kernel.
%
%   Y = SMOOTHN(X, SIZE, FILTER, STD) STD is a vector of standard deviations 
%       one for each dimension, when filter is 'gaussian' (default is 0.65)

%     $Author: ganil $
%     $Date: 2001/09/17 18:54:39 $
%     $Revision: 1.1 $
%     $State: Exp $

if nargin == 2,
  filt = 'b';
elseif nargin == 3,
  std = 0.65;
elseif nargin>4 | nargin<2
  error('Wrong number of input arguments.');
end

% check the correctness of sz
if ndims(sz) > 2 | min(size(sz)) ~= 1
  error('SIZE must be a vector');
elseif length(sz) == 1
  sz = repmat(sz,ndims(X));
elseif ndims(X) ~= length(sz)
  error('SIZE must be a vector of length equal to the dimensionality of X');
end

% check the correctness of std
if filt(1) == 'g'
  if length(std) == 1
    std = std*ones(ndims(X),1);
  elseif ndims(X) ~= length(std)
    error('STD must be a vector of length equal to the dimensionality of X');
  end
  std = std(:)';
end

sz = sz(:)';

% check for appropriate size
padSize = (sz-1)/2;
if ~isequal(padSize, floor(padSize)) | any(padSize<0)
  error('All elements of SIZE must be odd integers >= 1.');
end

% generate the convolution kernel based on the choice of the filter
filt = lower(filt);
if (filt(1) == 'b')
  smooth = ones(sz)/prod(sz); % box filter in N-D
elseif (filt(1) == 'g')
  smooth = ndgaussian(padSize,std); % a gaussian filter in N-D
else
  error('Unknown filter');
end


% pad the data
X = padreplicate(X,padSize);

% perform the convolution
Y = convn(X,smooth,'valid');

end

function h = ndgaussian(siz,std)

% Calculate a non-symmetric ND gaussian. Note that STD is scaled to the
% sizes in SIZ as STD = STD.*SIZ


ndim = length(siz);
sizd = cell(ndim,1);

for i = 1:ndim
  sizd{i} = -siz(i):siz(i);
end

grid = gridnd(sizd);
std = reshape(std.*siz,[ones(1,ndim) ndim]);
std(find(siz==0)) = 1; % no smoothing along these dimensions as siz = 0
std = repmat(std,2*siz+1);


h = exp(-sum((grid.*grid)./(2*std.*std),ndim+1));
h = h/sum(h(:));

end

function argout = gridnd(argin)

% exactly the same as ndgrid but it accepts only one input argument of 
% type cell and a single output array

nin = length(argin);
nout = nin;

for i=nin:-1:1,
  argin{i} = full(argin{i}); % Make sure everything is full
  siz(i) = prod(size(argin{i}));
end
if length(siz)<nout, siz = [siz ones(1,nout-length(siz))]; end

argout = [];
for i=1:nout,
  x = argin{i}(:); % Extract and reshape as a vector.
  s = siz; s(i) = []; % Remove i-th dimension
  x = reshape(x(:,ones(1,prod(s))),[length(x) s]); % Expand x
  x = permute(x,[2:i 1 i+1:nout]);% Permute to i'th dimension
  argout = cat(nin+1,argout,x);% Concatenate to the output 
end

end

function b=padreplicate(a, padSize)
% Pad an array by replicating values.
numDims = length(padSize);
idx = cell(numDims,1);
for k = 1:numDims
  M = size(a,k);
  onesVector = ones(1,padSize(k));
  idx{k} = [onesVector 1:M M*onesVector];
end

b = a(idx{:});

end

function [H,HNAN] = imagescnan(varargin)
%IMAGESCNAN   Scale data and display as image with uncolored NaNs.
%
%   SYNTAX:
%                imagescnan(U)
%                imagescnan(U,...,'NanColor',CNAN)
%                imagescnan(U,...,'NanMask',MNAN)
%                imagescnan(U,...,IOPT)
%                imagescnan(X,Y,U,...)
%     [H,HNAN] = imagescnan(...);
%
%   INPUT:
%     U    - 2 dimensional N-by-M image or N-by-M-by-3 RGB image.
%     X    - 2 extrema X-axis data; or the M values; or the N-by-M values
%            as obtained from MESHGRID (see DESCRIPTION below). 
%            DEFAULT: [1 N]
%     Y    - 2 extrema X-axis data; or the N values; or the N-by-M values
%            as obtained from MESHGRID (see DESCRIPTION below). 
%            DEFAULT: [1 M]
%     CNAN - Color for the NaNs elements. May be a char specifier or an [R
%            G B] triplet specifying the color.
%            DEFAULT: invisible (axes background color)
%     MNAN - Elements to be ignored besides not finite values. May be an
%            scalar or a logical M-by-N matrix indicating the elements to
%            be ignored.
%            DEFAULT: []
%     IOPT - IMAGE function normal optional pair arguments like
%            ('Parent',H) or/and CLIM like optional last argument as in
%            IMAGESC. 
%            DEFAULT: none
%
%   OUTPUT (all optional):
%     H    - Image handle
%     HNAN - Handle of every ignored (NaN) value colored patch.
%
%   DESCRIPTION:
%     MATLAB function IMAGESC does not work properly with NaNs. This
%     programs deals with this problem by including colored patches over
%     this elements and maybe others specyfied by the user with MNAN. 
%
%     Besides, those functions does not work properly with X,Y values
%     variable interval, but this functions does it by generating a whole
%     new image of several rectangular patches, but whose centers may not
%     lay in the specified coordinate (see NOTE below). This functionality
%     is experimental and not recommended (see ADDITIONAL NOTES inside this
%     program).
%
%     In previous release, 2-dim input images were transformed into a
%     3-dim RGB image. This is not used anymore (see ADDITIONAL NOTES
%     inside this file).
%
%   NOTE:
%     * Optional inputs use its DEFAULT value when not given or [].
%     * Optional outputs may or not be called.
%     * If X is a two element vector, min(X) will be the coordinate of the
%       first column and max(X) of the last column.
%     * If Y is a two element vector, min(Y) will be the coordinate of the
%       first row and max(Y) of the last row.
%     * If vector X-axis is decreasing U=fliplr(U) will be used.
%     * If vector Y-axis is decreasing U=flipud(U) will be used.
%     * When X or Y do not have a constant increasing/decreasing step, the
%       vertices of the color rectangules are set in the middle of each
%       pair of coordinates. For this reason its center may not lay on the
%       specified coordinate, except on the coordinates at the edges where
%       it always lays on the center.
%     * To get a non-scaled image (IMAGE instead of IMAGESC) use:
%         >> H = imagescnan(...);
%         >> set(H,'CDataMapping','direct')
%     * ADDITIONAL NOTES are included inside this file.
%
%   EXAMPLE:
%     % Compares with normal IMAGESC:
%      N     = 100;
%      PNaNs = 0.10;
%      U     = peaks(N);
%      U(round(1 + (N^2-1).*rand(N^2*PNaNs,1))) = NaN;         % Adds NaNs
%      subplot(221), imagesc(U)
%       title('With IMAGESC: ugly NaNs')
%      subplot(222), imagescnan(U) 
%       title('With IMAGESCNAN: uncolored NaNs')
%     % Compares with SPY:
%      subplot(223), spy(isnan(U))
%       title('SPY(isnan(U))')
%      subplot(224), imagescnan(isnan(U),'NaNMask',0), axis equal tight
%       title('SPY with IMAGESCNAN')
%     
%   SEE ALSO:
%     IMAGE, IMAGESC, COLORBAR, IMREAD, IMWRITE
%     and
%     CMAPPING, CBFREEZE by Carlos Vargas
%     at http://www.mathworks.com/matlabcentral/fileexchange
%
%
%   ---
%   MFILE:   imagescnan.m
%   VERSION: 2.1 (Aug 20, 2009) (<a href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/authors/11258')">download</a>) 
%   MATLAB:  7.7.0.471 (R2008b)
%   AUTHOR:  Carlos Adrian Vargas Aguilera (MEXICO)
%   CONTACT: nubeobscura@hotmail.com

%   ADDITIONAL NOTES:
%     * I keep getting a kind of BUG with the edges of the patched NaNs. I
%       added two NOTE inside this program that may fix this problem.
%       Another way is to convert the intensity matrix U into RGB colors by
%       using the CMAPPING function, as used by the first version of this
%       program.
%     * Besides, if the matrix is too large, sometimes there is an
%       undocumented failure while drawing the patch NaNs. Is recommended
%       to use U = cmapping(U,[],'k','discrete') instead, and change the
%       CLIM to [min(U(:)) max(U(:))].
%     * The use of not homogeneous step interval X,Y axes is not
%       recommended because the program tries to put its value in the
%       middle of the colored rectangule (as IMAGESC does) and soetimes the
%       result may not be what the user wants. So this is for experimental
%       use only.

%   REVISIONS:
%   1.0      Released. (Jun 30, 2008)
%   1.1      Fixed bug when CAXIS used. Colorbar freezed colormap. Fixed
%            bug in color vector input (Found by Greg King) and now 
%            accets RGB image as input. (Jul 14, 2008)
%   2.0      Totally rewritten code. Do not converts to RGB anymore. Do not
%            freezes the colormap anymore. Do not output any colorbar. New
%            X and Y variable steps accepted input. Now uses patches. (Jun
%            08, 2009)
%   2.1      Fixed bug with RGB input. Added a NOTE about the use of
%            CMAPPING. (Aug 20, 2009)

%   DISCLAIMER:
%   imagescnan.m is provided "as is" without warranty of any kind, under
%   the revised BSD license.

%   Copyright (c) 2008,2009 Carlos Adrian Vargas Aguilera


% INPUTS CHECK-IN
% -------------------------------------------------------------------------

% Initializes:
X    = [];
Y    = [];
CNAN = [];
MNAN = [];
ha   = [];

% Checks number of inputs:
if     nargin<1
 error('CVARGAS:imagescnan:notEnoughInputs',...
  'At least 1 input is required.')
elseif nargout>2
 error('CVARGAS:imagescnan:tooManyOutputs',...
  'At most 2 outputs are allowed.')
end

% Gets X,Y,U:
if ((nargin==1) || (nargin==2))
 U = varargin{1};
 varargin(1) = [];
else
 if (isnumeric(varargin{1}) && isnumeric(varargin{2}) && ...
   isnumeric(varargin{3}))
  X = varargin{1};
  Y = varargin{2};
  U = varargin{3};
  varargin(1:3) = [];
 else
  U = varargin{1};
  varargin(1) = [];
 end
end

% Check U:
ndim = ndims(U);
if     (ndim==2)
 [M,N]   = size(U);
 O = 1;
elseif (ndim==3)
 [M,N,O] = size(U);
 if (O~=3)
  error('CVARGAS:imagescnan:incorrectRgbImage',...
   'RGB image must be of size M-by-N-by-3.')
 end
else
 error('CVARGAS:imagescnan:incorrectImageSize',...
  'Image must be 2-dimensional or a 3-dim RGB image.')
end

% Check X:
aequal = true;    % Equal intervals on x-axis?
dX     = [];
if isempty(X)
 X = [1 N];
else
 if (ndims(X)>2)
  error('CVARGAS:imagescnan:incorrectXDims',...
   'X must be a vector or a matrix as a result of MESHGRID.')
 end
 if any(~isfinite(X(:)))
  error('CVARGAS:imagescnan:incorrectXValue',...
   'X elements must be numeric and finite.')
 end
 [Mx,Nx] = size(X);
 if ((Mx*Nx)==2)
  if X(2)<X(1)
   X = X([2 1]);
   for k = 1:O % Fixed bug Aug 2009
    U(:,:,k) = fliplr(U(:,:,k));
   end
  end 
 else
  if     ((Mx==M) && (Nx==N))
   % Checks if generated with MESHGRID:
   dX    = abs(X(2:M,:)-repmat(X(1,:),M-1,1));
   if any(abs(dX(:))>(eps*max(abs(dX(:)))*1000))
    error('CVARGAS:imagescnan:incorrectXMatrix',...
     'X matrix must be as generated by MESHGRID.')
   end
   X = X(1,:);
  elseif (~any([Mx Nx]==1) && ~((Mx*Nx)==N))
   error('CVARGAS:imagescnan:incorrectXSize',...
     'X must be an scalar or a matrix.')
  end     
  % Forces ascending x-axis:
  [X,I] = sort(X(:).');
  for k = 1:O % Fixed bug Aug 2009
   U(:,:,k) = U(:,I,k);
  end
  clear I
  % Checks equal intervals:
  dX = diff(X);
  if any(abs(dX(1)-dX(2:end))>(eps*max(dX)*1000))
   if aequal
    aequal = false;
   end
  else
   X  = [X(1) X(end)];
   dX = [];
  end
 end
end

% Check Y:
dY = [];
if isempty(Y)
 Y = [1 M];
else
 if (ndims(Y)>2)
  error('CVARGAS:imagescnan:incorrectYDims',...
   'Y must be a vector or a matrix as a result of MESHGRID.')
 end
 if any(~isfinite(Y(:)))
  error('CVARGAS:imagescnan:incorrectYValue',...
   'Y elements must be numeric and finite.')
 end
 [My,Ny] = size(Y);
 if ((My*Ny)==2)
  if Y(2)<Y(1)
   Y = Y([2 1]);
   for k = 1:O % Fixed bug Aug 2009
    U(:,:,k) = flipud(U(:,:,k));
   end
  end
 else
  if     ((My==M) && (Ny==N))
   % Checks if generated with MESHGRID:
   dY = abs(Y(:,2:N)-repmat(Y(:,1),1,N-1));
   if any(abs(dY(:))>(eps*max(abs(dY(:)))*1000))
    error('CVARGAS:imagescnan:incorrectYMatrix',...
     'Y matrix must be as generated by MESHGRID.')
   end
   Y = Y(:,1);
  elseif (~any([My Ny]==1) && ~((My*Ny)==M))
   error('CVARGAS:imagescnan:incorrectYSize',...
     'Y must be an scalar or a matrix.')
  end     
  % Forces ascending y-axis:
  [Y,I] = sort(Y(:).');
  for k = 1:O % Fixed bug Aug 2009
   U(:,:,k) = U(I,:,k);
  end
  clear I
  % Checks equal intervals:
  dY = diff(Y);
  if any(abs(dY(1)-dY(2:end))>(eps*max(dY)*1000))
   if aequal
    aequal = false;
   end
  else
   Y  = [Y(1) Y(end)];
   dY = [];
  end
 end
end

% Checks varargin:
ind  = [];
Nopt = length(varargin); 
for k = 1:Nopt-1
 if (~isempty(varargin{k}) && ischar(varargin{k}))
  if     strncmpi(varargin{k},'NanColor',4)
   CNAN = varargin{k+1};
   ind  = [ind k k+1];
  elseif strncmpi(varargin{k},'NanMask',4)
   MNAN = varargin{k+1};
   ind  = [ind k k+1];
  elseif (strncmpi(varargin{k},'Parent',2) && isempty(CNAN))
   try
    CNAN = get(varargin{k+1},'Color');
    ha   = varargin{k+1};
   catch
    error('CVARGAS:imagescnan:incorrectParentHandle',...
     '''Parent'' must be a valid axes handle.')
   end
  end
 end
end
varargin(ind) = [];
Nargin = length(varargin);

% Check ha:
if isempty(ha)
 ha = gca;
end

% Check CNAN:
if     isempty(CNAN)
 CNAN = get(ha,'Color');
elseif ischar(CNAN)
 switch lower(CNAN)
  case 'y', CNAN = [1 1 0];
  case 'm', CNAN = [1 0 0];
  case 'c', CNAN = [0 1 1];
  case 'r', CNAN = [1 0 0];
  case 'g', CNAN = [0 1 0];
  case 'b', CNAN = [0 0 1];
  case 'w', CNAN = [1 1 1];
  case 'k', CNAN = [0 0 0];
  otherwise
   error('CVARGAS:imagescnan:incorrectNancString',...
    'Color string must be a valid color identifier. One of ''ymcrgbwk''.')
 end
elseif isnumeric(CNAN) && (length(CNAN)==3)
 CNAN = CNAN(:).'; % Forces row vector.
else
 error('CVARGAS:imagescnan:incorrectNancInput',...
  'Not recognized CNAN input.')
end

% Check MNAN:
if isempty(MNAN)
 MNAN = any(~isfinite(U),3);
else
 if (ndims(MNAN)==2)
  [Mm,Nm] = size(MNAN);
  if     ((Mm*Nm)==1)
   MNAN = (any(~isfinite(U),3) | any(U==MNAN,3));
  elseif ((Mm==M) && (Nm==N) && islogical(MNAN))
   MNAN = (any(~isfinite(U),3) | MNAN);
  else
   error('CVARGAS:imagescnan:incorrectNanmSize',...
   'MNAN must be an scalar or a logical matrix of size M-by-N.')
  end
 else
  error('CVARGAS:imagescnan:incorrectNanmDims',...
   'MNAN must be an scalar or a matrix.')
 end
end

% -------------------------------------------------------------------------
% MAIN
% -------------------------------------------------------------------------

% Generates the image:
if aequal
 % IMAGESC way.
 H = imagesc(X,Y,U,varargin{:});

else
 % PATCH way.

 % Check clim:
 if (rem(Nargin,2)==1)
  clim          = varargin{end};
  varargin(end) = [];
  if ((length(clim)~=2) || (clim(1)>clim(2)))
   error('CVARGAS:imagescnan:incorrectClimInput',...
    'clim must be a 2 element increasing vector.')
  end
 else
  clim = [];
 end

 % Generates vertices between coordinates (coordinates may not be at the
 % center of these vertices): 
 if (length(X)~=N)
  X  = (0:N-1)*((X(2)-X(1))/(N-1)) + X(1);
 end
 if (length(Y)~=M)
  Y  = (0:M-1)*((Y(2)-Y(1))/(M-1)) + Y(1);
 end
 if isempty(dX)
  dX = diff(X);
 end
 if isempty(dY)
  dY = diff(Y);
 end
 [X,Y] = meshgrid([X(1)-dX(1)/2 X+dX([1:N-1 N-1])/2],...
                  [Y(1)-dY(1)/2 Y+dY([1:M-1 M-1])/2]);
 
 % Generates faces:
 ind              = (1:(M+1)*N)';
 ind(M+1:M+1:end) = [];
 
 % Generates patches:
 H = patch(...
  'Vertices'       ,[X(:) Y(:)],...
  'Faces'          ,[ind ind+1 ind+M+2 ind+M+1],...
  'FaceVertexCData',U(:),...
  'FaceColor'      ,'flat',...
  'EdgeColor'      ,'none',... % NOTE: Sometimes this is not required.
  varargin{:});
 set(ha,...
  'YDir' ,'reverse',...
  'View' ,[0 90],...
  'Box'  ,'on',...
  'Layer','top')
 axis(ha,'tight')
 
 % Sets clim:
 if ~isempty(clim)
  set(ha,'CLim',clim)
 else
  set(ha,'CLimMode','auto')
 end
 
end

% Adds NaNs patches:
if any(MNAN(:))
 if aequal
  % dX and dY is constant:
  [MNAN,NNAN] = ind2sub([M,N],find(MNAN));
  Nnan        = length(MNAN);
  dX   = (X(2)-X(1))/(N-1)/2;
  dY   = (Y(2)-Y(1))/(M-1)/2;
  HNAN = patch(repmat((X(1)+(NNAN(:)'-1)*(2*dX)),4,1) + ...
                                       (repmat([-1 1 1 -1]'*dX,1,Nnan)),...
               repmat((Y(1)+(MNAN(:)'-1)*(2*dY)),4,1) + ...                                       
                                       (repmat([1 1 -1 -1]'*dY,1,Nnan)),...
               CNAN,...
               'EdgeColor',CNAN,... 'EdgeColor','none',... 
               varargin{1:Nargin-rem(Nargin,2)});
 else
  % dX and/or dY is not constant:
  MNAN = find(MNAN);
  HNAN = patch(...
  'Vertices'       ,[X(:) Y(:)],...
  'Faces'          ,[ind(MNAN) ind(MNAN)+1 ind(MNAN)+M+2 ind(MNAN)+M+1],...
  'FaceColor'      ,CNAN,...
  'EdgeColor'      ,'none',... 'EdgeColor',CNAN,... % NOTE: may be better?
  varargin{:});
 end
else
 HNAN = [];
end

 
% OUTPUTS CHECK-OUT
% -------------------------------------------------------------------------

% Clears outputs?:
if (nargout==0)
 clear H
end
end