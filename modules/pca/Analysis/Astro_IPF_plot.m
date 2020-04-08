% TPM 20/02/2020 for AstroEBSD_V2

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

function []=Astro_IPF_plot(EBSD_template,printing)

% Now plot the MTEX data
setMTEXpref('xAxisDirection','west');
setMTEXpref('zAxisDirection','OutOfPlane');

for a=1:length(EBSD_template.indexedPhasesId);
    phaseno=EBSD_template.indexedPhasesId(a);
    Mapping.oM(a)=ipfHSVKey(EBSD_template(EBSD_template.mineralList(phaseno)));
    clear phaseno
end

% Plot IPFX
figure;
for a=1:length(EBSD_template.indexedPhasesId);
    phaseno=EBSD_template.indexedPhasesId(a);
    Mapping.oM(a).inversePoleFigureDirection = xvector;
    plot(EBSD_template(EBSD_template.mineralList(phaseno)), Mapping.oM(a).orientation2color(EBSD_template(EBSD_template.mineralList(phaseno)).orientations),'MicronBar','off')
    hold on
end
hold off

if printing==1
    print(gcf,'IPFX','-dpng','-r500')
end

% Plot IPFY
figure;
for a=1:length(EBSD_template.indexedPhasesId);
    phaseno=EBSD_template.indexedPhasesId(a);
    Mapping.oM(a).inversePoleFigureDirection = yvector;
    plot(EBSD_template(EBSD_template.mineralList(phaseno)), Mapping.oM(a).orientation2color(EBSD_template(EBSD_template.mineralList(phaseno)).orientations),'MicronBar','off')
    hold on
end
hold off

if printing==1
    print(gcf,'IPFY','-dpng','-r500')
end

% Plot IPFZ
figure;
for a=1:length(EBSD_template.indexedPhasesId);
    phaseno=EBSD_template.indexedPhasesId(a);
    Mapping.oM(a).inversePoleFigureDirection = zvector;
    plot(EBSD_template(EBSD_template.mineralList(phaseno)), Mapping.oM(a).orientation2color(EBSD_template(EBSD_template.mineralList(phaseno)).orientations),'MicronBar','off')
    hold on
end
hold off

if printing==1
    print(gcf,'IPFZ','-dpng','-r500')
end


% Plot colour keys
for a=1:length(EBSD_template.indexedPhasesId);
    figure
    phaseno=EBSD_template.indexedPhasesId(a);
    plot(Mapping.oM(a))
    title(EBSD_template.mineralList(phaseno))
    if printing==1
        print(gcf,['IPFkey-',EBSD_template.mineralList{phaseno}],'-dpng','-r300')
    end
    end
end

