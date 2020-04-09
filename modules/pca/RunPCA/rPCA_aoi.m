% rPCA_aoi - T P McAuliffe 20/02/2020 for AstroEBSD v2

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

% This runs PCA on testArray, whether it is composed of EBSD and/or EDX...


%  run PCA analysis to split data into basis patterns (eigenvectors) and signal strength in map (loadings)
%   this uses a user input choice on number of factors to retain (NumPCA)

%coeff = maps
%scores = patterns

%% perform the PCA
if PCA_Setup.weighted==1;
    weightings=[Output.weighting_ratio*ones(1,EBSD_Info.PatSizeW*EBSD_Info.PatSizeH),ones(1,Settings_Cor.channum)];
    [Output.coeff,Output.score,Output.latent,~,Output.explained]=pca(testArray, 'Centered',false,'Weights',weightings);
else
    [Output.coeff,Output.score,Output.latent,~,Output.explained]=pca(testArray, 'Centered',false);%,'NumComponents',200);
end

clear weightings
%coeff: returns all principal components. Each column contains coefficients
%for each principal component.
%score: princ. components scores are representations of input (testArray2)
%in principal component space. Rows = observations; columns = components. 
%latent: variances of principal components
