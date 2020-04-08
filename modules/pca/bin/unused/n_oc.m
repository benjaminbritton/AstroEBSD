plot(Output.explained,'r','LineWidth',1.5)
xlim([0,30])
xlabel('Principal Component')
set(gca,'YScale','log')
ylabel('% contribution to dataset variance')

p=200;
l=Output.latent;

%%%% See paper by Raiche et al


nK=min(find(Output.latent<1));

%%
for i=1:200
b=(l(p)-l(i+1))./(p-i-1);
b2=(l(p)-l(i+2))./(p-i);
a=l(i+1)-b2;

L(i)=a+b;
end

for i=1:200
    difference(i)=L(i)-l(i);
end

noc=min(find(difference>-0.001));


%%

for i=2:199;
    af(i)=L(i+1)-2.*L(i)-L(i-1)
end

naf=min(find(af>-0.01));

%%

