%% Matched RC-pattern comparison to simulations - for AstroEBSD v2 - T P McAuliffe 19/02/20

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

clear all
home

%Choose the input patterns (RC-EBSPs) you want to compare to library.
%These will be created from the PCA_deck file
load('I:\TomMcA\PCA_output\SuperalloyExample_EBSD+EDS\weighting1_vt0,1_Guo_rad3\Results.mat');

% Set up and locate plugins
InputUser.MTEX_loc='C:\CommunalMatlab_plugins\mtex-5.2.beta2';
InputUser.Astro_loc='I:\TomMcA\GitHub\AstroEBSD_v2';

run([InputUser.MTEX_loc,'\startup_mtex']);
run([InputUser.Astro_loc,'\start_AstroEBSD']);

%%
InputPats=tile.pat; %where should the patterns come from?

component=1; %what component number should we compare patterns for?
% These can visually be chosen by plotting tile.map_reshaped

patcomp.print=0; %print the matched patterns to the CURRENT directory

subplot_h=1;
subplot_w=4; % h * w must be larger than your input phases length(InputUser.Phases)
PCA_patcomparison %Run the script


