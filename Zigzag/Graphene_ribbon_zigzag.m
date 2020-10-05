%Deleting unimportant data
    clear
    clc
    close all
%-------------------------------------------
%initializing parameters & Bravais vectors
a1 = 0.5*[-sqrt(3),3];
a2 = 0.5*[sqrt(3),3];
a = sqrt(3)*norm(a1); %unit cell width
%%
%Generating Hamiltonian
n = 19; %Number of hexagon cells in unit cell

N = 40; %number of atoms in unit cell
        %Parameters:
            t = 1.;
            L = 0.1;
            V = 0.0;
        %
grid = 200;
K = -2*pi/a:pi/grid/a:2*pi/a;
E = zeros(length(K),N);
figure(1);
% filename = 'ChangingV_L=0,1.gif';
% for V=-0.5:0.05:0.5
    for g = 1:length(K)
        k = K(g);
        u = exp(i*k*a); %u is phase going right
            H_NN = zeros(N,N);
            H_NNN = zeros(N,N);
            
            H_NN(N,N-1) = t*(1+conj(u));  
            H_NNN(N-1,N-1) = -V+2*L*sin(k*a);
            H_NNN(N,N) = V-2*L*sin(k*a);
            
        for ii=1:N-2
           %Diagonal elements
            if(mod(ii,2)==0)
                H_NNN(ii,ii) = V-2*L*sin(k*a); %iL( exp(ika)-exp(-ika) )
            else
                H_NNN(ii,ii) = -V+2*L*sin(k*a); %iL( exp(-ika)-exp(ika) )
            end

            if(mod(ii,4)==1)
                H_NN(ii,ii+1) = t*(1+conj(u)); %1->2
                H_NNN(ii,ii+2) = i*L*(1 - conj(u)); % 1->3
            end
            if(mod(ii,4)==2)
                H_NN(ii,ii+1) = t; %2->3
                H_NNN(ii,ii+2) = i*L*(1 - u); %2->4
            end
            if(mod(ii,4)==3)
                H_NN(ii,ii+1) = t*(1+u); %3->4
                H_NNN(ii,ii+2) = -i*L*(1 - u); %3->5
            end
            if(mod(ii,4)==0)
                H_NN(ii,ii+1) = t; %4->5
                H_NNN(ii,ii+2) = -i*L*(1 - conj(u)); % 4->6
            end
        end
        H = H_NN + H_NNN;
        H = H + H';
        E(g,:) = eig(H);
    end
    %%
    %Plotting energy
    for ii=1:N
        if((ii==N/2) || (ii==N/2+1))
            plot(K,E(:,ii),'r-');
        else
            plot(K,E(:,ii),'k-');
        end
        hold on
    end
    drawnow
    %title(sprintf('Energy spectrum for graphene ribbon with zizzag edge for %0.0f atoms width\n using parameters: t=%0.2f, \\lambda=%0.2f and V=%0.2f \n Topological phase transition',N,t,L,V)); 
    title(sprintf('Energy spectrum for graphene ribbon with zizzag edge for %0.0f atoms width\n using parameters: t=%0.2f, \\lambda=%0.2f and V=%0.2f',N,t,L,V)); drawnow
    xlabel(sprintf('k [\\pi/a]')); drawnow
    ylabel(sprintf('E [eV]')); drawnow
    xticks([-3/2*pi -pi -pi/2 0 pi/2 pi 3/2*pi]/a); drawnow
    xticklabels({'-3/2','-1','-1/2','0','1/2','1','3/2'}); drawnow
    axis([min(K) max(K) min(E(:,1)) max(E(:,N))]); drawnow
    hold off
%       % Creating .gif file for animation
%    frame = getframe(1);
%    im = frame2im(frame);
%    [imind, cm] = rgb2ind(im,256);
%    if V == -0.5;
%          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%    else
%          imwrite(imind,cm,filename,'gif','WriteMode','append');
%    end
% end








