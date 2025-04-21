clear all
close all

 for n=1:15
     for m=0:n

close all


%Initialize
k_steps=10^5;
Hist_steps=k_steps/50;
gamma=2.74;
%gamma=1;

%Unit vectors in real space
a_1=[1/2,sqrt(3)/2];
a_2=[-1/2,sqrt(3)/2];

%Unit vectors in k-space
b_1=[2*pi,2*pi/sqrt(3)];
b_2=[-2*pi,2*pi/sqrt(3)];

%Set to 1 if you want to generate .html file
do_save=0;

%Count unit cells
d_r=gcd(2*n+m,n+2*m);
N_CNT=(2*(n^2+m^2+n*m))/d_r;

%Translation vectors
t_1=(2*m+n)/d_r;
t_2=-(2*n+m)/d_r;
T=t_1.*a_1+t_2.*a_2;

%Graphene K-points
K_1=(1/N_CNT).*(-t_2.*b_1+t_1.*b_2);
K_2=(1/N_CNT).*(m.*b_1-n.*b_2);
K_1_norm=norm(K_1);
K_2_norm=norm(K_2);

l_2=linspace(-pi/norm(T),pi/norm(T),k_steps);

E_p=nan(N_CNT,k_steps);
E_n=E_p;

for l_1=0:N_CNT
    
    %Allowed k-values
    k_x=l_1.*K_1(1)+l_2.*K_2(1)./K_2_norm;
    k_y=l_1.*K_1(2)+l_2.*K_2(2)./K_2_norm;
    
    %Calculate Dispersion Relation
    val=gamma*sqrt(1+4*cos(sqrt(3)*k_y/2).*cos(k_x/2)+4*cos(k_x/2).^2);
    E_p(l_1+1,:)=val;
    E_n(l_1+1,:)=-val;
    
end

 %Plot Dispersion Relation
%fig=figure('visible','off');
fig=figure;
subplot(1,2,1)
plot(linspace(-1,1,k_steps),E_p,'-r',linspace(-1,1,k_steps),E_n,'-b')
xlabel('k/k_{max}')
ylabel('E [eV]')
xlim([-1,1])
ylim([-3.2*gamma,3.2*gamma])


% %Plot DOS
E=[E_p;E_n];
[D_E_per_bin_per_a,edges]=histcounts(E,Hist_steps);
bin_lenghts=diff(edges);
x_DOS=edges(1:end-1)+bin_lenghts(1)/2;
D_E_per_eV_per_m=D_E_per_bin_per_a*2.46*10^10*bin_lenghts(1);
subplot(1,2,2)
stairs(x_DOS,D_E_per_eV_per_m)
xlabel('E [eV]');
ylabel('D(E) [ eV^{-1} m^{-1}]')
ylim([0,5*median(D_E_per_eV_per_m)])
xlim([-3*gamma,3*gamma])

set(gcf, 'Position', [300, 300, 1000, 500])



%Make html file
if do_save==1
 %set to the folder you want to save in
 folder = 'C:\Users\Sebastian\Desktop\student project\new_2';
 filename_html=[folder,'\',num2str(n),'_',num2str(m),'.html'];
 filename_fig=[folder,'\',num2str(n),'_',num2str(m),'.png'];
 htmltags_1=['<html> <body> <h2 style="text-align:center">Dispersion relation and density of states of a ' ,num2str(n),',' num2str(m),' carbon nanotube</h2>'];
 htmltags_2=['<p style="text-align:center"><img  src="',num2str(n),'_' num2str(m),'.png" /></p>'];
 htmltags_3=['<p style="text-align:center">E &nbsp;[eV]&emsp;&emsp;<i>D</i>(E)&nbsp;[eV<sup>-1</sup>&nbsp;m<sup>-1</sup>] <br> <textarea rows="20" cols="40" id="output">'];
data_t_n=[x_DOS(1:end/2)',D_E_per_eV_per_m(1:end/2)'];
data_t_p=[x_DOS(end/2:end)',D_E_per_eV_per_m(end/2:end)'];
dlmwrite(filename_html,[htmltags_1,htmltags_2,htmltags_3],'delimiter','')
dlmwrite(filename_html,data_t_n,'delimiter','\t','-append','precision','%.5f')
dlmwrite(filename_html,data_t_p,'delimiter','\t','-append','precision','%.6f')
htmltags_4=['</textarea></p></body></html>'];
dlmwrite(filename_html,htmltags_4,'delimiter','', '-append');
print(filename_fig,'-dpng')
else
end

end

end
