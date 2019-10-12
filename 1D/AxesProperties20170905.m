if Parameter1 == 1
    set(gca,'fontsize',13);
    ylabel('{\bf J_{actin}}/\theta_* (nm/s)','fontsize',13)
elseif Parameter1 == 2
    set(gca,'fontsize',13,'yscale','log');
    ylabel('\eta_{st} (Pa s/{\mu}m^2)','fontsize',13)
elseif Parameter1 == 3
    set(gca,'fontsize',13,'yscale','log');
    ylabel('\xi_{m} (Pa s/{\mu}m)','fontsize',13)
elseif Parameter1 == 4
    set(gca,'fontsize',13);
    ylabel('c_0^f (mM)','fontsize',13)
elseif Parameter1 == 5
    set(gca,'fontsize',13,'yscale','log');
%     xticks(logspace(-1,3,5)*1d9/1d6);
    ylabel('d_g (Pa s/{\mu}m)','fontsize',13)
elseif Parameter1 == 6
    set(gca,'fontsize',13,'yscale','log');
    ylabel('\phi_n^f (Pa)','fontsize',13);
elseif Parameter1 == 7
    set(gca,'fontsize',13);
    ylabel('\Theta_n','fontsize',13);
elseif Parameter1 == 8
    set(gca,'fontsize',13,'yscale','log');
    ylabel('\eta (Pa s/{\mu}m^2)','fontsize',13)
elseif Parameter1 == 9
    set(gca,'fontsize',13,'yscale','log');
    ylabel('\eta_n (Pa s/m^2)','fontsize',13)
elseif Parameter1 == 10
    set(gca,'fontsize',13,'yscale','log');
    ylabel('\zeta_c (Pa s/m^2)','fontsize',13)
elseif Parameter1 == 11
    set(gca,'fontsize',13,'yscale','log');
    ylabel('\zeta_n (Pa s/m^2)','fontsize',13)
elseif Parameter1 == 12
    set(gca,'fontsize',13);
    ylabel('f_{ext}^f (Pa)','fontsize',13)
elseif Parameter1 == 13
    set(gca,'fontsize',13);
    ylabel('f_{ext}^b (Pa)','fontsize',13)
elseif Parameter1 == 14
    set(gca,'fontsize',13);
    ylabel('J_{water} (nm/s)','fontsize',13)
elseif Parameter1 == 18
    set(gca,'fontsize',13);
    ylabel('J_{c,active}^f (mM nm/s)','fontsize',13)
elseif Parameter1 == 19
    set(gca,'fontsize',13);
    ylabel('p_0^b (Pa)','fontsize',13)
elseif Parameter1 == 20
    set(gca,'fontsize',13,'yscale','log');
    ylabel('\gamma_m (Pa s/{\mu}m^2)','fontsize',13)
elseif Parameter1 == 21
    set(gca,'fontsize',13,'yscale','log');
    ylabel('\alpha^{f,b} ({\mu}m/Pa/s)','fontsize',13)
elseif Parameter1 == 22
    set(gca,'fontsize',13,'yscale','log');
    ylabel('-\sigma_{a,0} (Pa)','fontsize',13)
elseif Parameter1 == 23
    set(gca,'fontsize',13,'yscale','log');
    ylabel('\xi_{m,1} (Pa s/{\mu}m^2)','fontsize',13)
elseif Parameter1 == 24
    set(gca,'fontsize',13,'yscale','log');
    ylabel('\eta_\theta (1/s)','fontsize',13)
elseif Parameter1 == 25
    set(gca,'fontsize',13,'yscale','log');
    ylabel('\eta_\theta^f (1/s)','fontsize',13)
end
    

if Parameter2 == 1
    set(gca,'fontsize',13);
    xlabel('{\bf J_{actin}}/\theta_* (nm/s)','fontsize',13)
elseif Parameter2 == 2
    set(gca,'fontsize',13,'xscale','log');
    xlabel('\eta_{st} (Pa s/{\mu}m^2)','fontsize',13)
elseif Parameter2 == 3
    set(gca,'fontsize',13,'xscale','log');
    xlabel('\xi_{m} (Pa s/{\mu}m)','fontsize',13)
elseif Parameter2 == 4
    set(gca,'fontsize',13);
    xlabel('c_0^f (mM)','fontsize',13)
elseif Parameter2 == 5
    set(gca,'fontsize',13,'xscale','log');
%     xticks(logspace(-1,3,5)*1d9/1d6);
    xlabel('d_g (Pa s/{\mu}m)','fontsize',13)
elseif Parameter2 == 6
    set(gca,'fontsize',13,'xscale','log');
    xlabel('\phi_n^f (Pa)','fontsize',13);
elseif Parameter2 == 7
    set(gca,'fontsize',13);
    xlabel('\Theta_n','fontsize',13);
elseif Parameter2 == 8
    set(gca,'fontsize',13,'xscale','log');
    xlabel('\eta (Pa s/{\mu}m^2)','fontsize',13)
elseif Parameter2 == 9
    set(gca,'fontsize',13,'xscale','log');
    xlabel('\eta_n (Pa s/m^2)','fontsize',13)
elseif Parameter2 == 10
    set(gca,'fontsize',13,'xscale','log');
    xlabel('\zeta_c (Pa s/m^2)','fontsize',13)
elseif Parameter2 == 11
    set(gca,'fontsize',13,'xscale','log');
    xlabel('\zeta_n (Pa s/m^2)','fontsize',13)
elseif Parameter2 == 12
    set(gca,'fontsize',13);
    xlabel('f_{ext}^f (Pa)','fontsize',13)
elseif Parameter2 == 13
    set(gca,'fontsize',13);
    xlabel('f_{ext}^b (Pa)','fontsize',13)
elseif Parameter2 == 14
    set(gca,'fontsize',13);
    xlabel('J_{water} (nm/s)','fontsize',13)
elseif Parameter2 == 18
    set(gca,'fontsize',13);
    xlabel('J_{c,active}^f (mM nm/s)','fontsize',13)
elseif Parameter2 == 19
    set(gca,'fontsize',13);
    xlabel('p_0^b (Pa)','fontsize',13)
elseif Parameter2 == 20
    set(gca,'fontsize',13,'xscale','log');
    xlabel('\gamma_m (Pa s/{\mu}m^2)','fontsize',13)
elseif Parameter2 == 21
    set(gca,'fontsize',13,'xscale','log');
    xlabel('\alpha^{f,b} ({\mu}m/Pa/s)','fontsize',13)
elseif Parameter2 == 22
    set(gca,'fontsize',13,'xscale','log');
    xlabel('-\sigma_{a,0} (Pa)','fontsize',13)
elseif Parameter2 == 23
    set(gca,'fontsize',13,'xscale','log');
    xlabel('\xi_{m,1} (Pa s/{\mu}m^2)','fontsize',13)
elseif Parameter2 == 24
    set(gca,'fontsize',13,'xscale','log');
    xlabel('\eta_\theta (1/s)','fontsize',13)
elseif Parameter2 == 25
    set(gca,'fontsize',13,'xscale','log');
    xlabel('\eta_\theta^f (1/s)','fontsize',13)
end


    
    vaxis = axis;
    vaxis(1) = min(min(Xmesh));
    vaxis(2) = max(max(Xmesh));
    vaxis(3) = min(min(Ymesh));
    vaxis(4) = max(max(Ymesh));
    axis(vaxis);
    
