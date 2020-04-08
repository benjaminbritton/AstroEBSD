%% PCA & RTM postanalysis - TPM 25/10/19

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

clear all
home

% Set up and locate plugins
InputUser.MTEX_loc='C:\CommunalMatlab_plugins\mtex-5.2.beta2';
InputUser.Astro_loc='I:\TomMcA\GitHub\AstroEBSD_v2';

run([InputUser.MTEX_loc,'\startup_mtex']);
run([InputUser.Astro_loc,'\start_AstroEBSD']);

% Can provide a list of files, or just specify a directory
list={'I:\TomMcA\PCA_output\SuperalloyExample_EBSD+EDS'};

%folders=["I:\TomMcA\PCA_output\trialCr+Ti_EBSD+EDS"];
%cd(folders)
%list=dir('*EBSD*');
% ind_forbidden=[10,12,19]; %if you want to skip any files in the list, eg. in a big directory

%%
%inds=1;%comment this out if you want to run for all datasets in list, otherwise it will just do this index

namesetting='weighting1_vt0,1_Guo_rad3'; %the name of the analyses you want to do on this list

%analyse single tile? Unless you specifically want to do this keep as 0.
post.singletile=0; %0 by default ; May need to be turned on for older PCA
%analysis datasets

%What to plot in the post-analysis
post.EBSD=1; %plot IPFs and colour keys

% IMPORTANT!! %
% In current implementation RCSpectra and AvSpectra require quantified 'results.xls' table in their
% saved directories (in Bruker export style). Example is provided in the
% repository

post.RCSpectra=1; %plot RC spectra quant
post.RCSpectra_std=1; %include standard deviation lines on the radar plot
post.AvSpectra=1; %plot average spectra quant
post.AvSpectra_std=1; %include standard deviation lines on the radar plot

post.comb=[4,6,7;4,6,9;7,4,9]; %element combinations for combo EDS maps %[r1,g1,b1;r2,g2,b2;r3,g3,b3];

post.chemlims.minV='autoscale'; %if this is a string this will autoscale the radar plots, otherwise specify in a list: eg. [0,0,0,...] for [C, Zr, Ti...]
post.chemlims.maxV='autoscale'; %as above string or list: eg. [10,12,14,...] for [C, Zr, Ti...] 
post.phase=1; %plot a phase map

%set up colourmap
cmap=cbrewer('qual','Dark2',15); %can be useful to change the number in here dep on how many phases you are using

%%
%%%%%
PCA_postanalysis %runs the postanalysis

