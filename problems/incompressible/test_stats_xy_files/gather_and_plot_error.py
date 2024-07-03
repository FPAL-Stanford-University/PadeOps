import numpy as np
import matplotlib.pyplot as plt
import code
import sys

if __name__=='__main__':

    if len(sys.argv) > 1:
        interscale_only = True if sys.argv[1].lower() == 'true' else False
    else:
        interscale_only = False

    # Domain setup
    nz = 128
    zmin = -4
    zmax = 2
    dz = (zmax-zmin)/nz
    z = zmin+dz/2 + np.arange(0,nz)*dz
    for scl in range(2):
        # First, plot the PadeOps and MATLAB implementations for visual comparison
        datadir = '/home1/06632/ryanhass/codes/PadeOps/problems/incompressible/test_stats_xy_files/data'

        stats = np.genfromtxt(datadir+'/Run78_stats_xy_t000004_scale'+str(scl+1).zfill(2)+'.out')
        stats_sca = np.genfromtxt(datadir+'/Run78_stats_xy_sca01_t000004_scale'+str(scl+1).zfill(2)+'.out')

        turb_trans_X_PO = stats[:,40]
        turb_trans_Y_PO = stats[:,55]
        turb_trans_Z_PO = stats[:,70]
        turb_trans_T_PO = stats_sca[:,3]
        turb_trans_wT_PO = stats_sca[:,23]

        turb_trans_mixed_X_PO = stats[:,50]
        turb_trans_mixed_Y_PO = stats[:,65]
        turb_trans_mixed_Z_PO = stats[:,80]
        turb_trans_mixed_T_PO = stats_sca[:,10]
        turb_trans_mixed_wT_PO = stats_sca[:,36]

        interscale_X_PO = stats[:,48]
        interscale_Y_PO = stats[:,63]
        interscale_Z_PO = stats[:,78]
        interscale_T_PO = stats_sca[:,8]
        interscale_wT_PO = stats_sca[:,34]

        turb_trans_X_MATLAB = np.genfromtxt(datadir+'/turb_trans_X_scl'+str(scl+1).zfill(2)+'.dat')
        turb_trans_Y_MATLAB = np.genfromtxt(datadir+'/turb_trans_Y_scl'+str(scl+1).zfill(2)+'.dat')
        turb_trans_Z_MATLAB = np.genfromtxt(datadir+'/turb_trans_Z_scl'+str(scl+1).zfill(2)+'.dat')
        turb_trans_T_MATLAB = np.genfromtxt(datadir+'/turb_trans_T_scl'+str(scl+1).zfill(2)+'.dat')
        turb_trans_wT_MATLAB = np.genfromtxt(datadir+'/turb_trans_wT_scl'+str(scl+1).zfill(2)+'.dat')

        turb_trans_mixed_X_MATLAB = np.genfromtxt(datadir+'/turb_trans_mixed_X_scl'+str(scl+1).zfill(2)+'.dat')
        turb_trans_mixed_Y_MATLAB = np.genfromtxt(datadir+'/turb_trans_mixed_Y_scl'+str(scl+1).zfill(2)+'.dat')
        turb_trans_mixed_Z_MATLAB = np.genfromtxt(datadir+'/turb_trans_mixed_Z_scl'+str(scl+1).zfill(2)+'.dat')
        turb_trans_mixed_T_MATLAB = np.genfromtxt(datadir+'/turb_trans_mixed_T_scl'+str(scl+1).zfill(2)+'.dat')
        turb_trans_mixed_wT_MATLAB = np.genfromtxt(datadir+'/turb_trans_mixed_wT_scl'+str(scl+1).zfill(2)+'.dat')

        interscale_X_MATLAB = np.genfromtxt(datadir+'/interscale_X.dat')
        interscale_Y_MATLAB = np.genfromtxt(datadir+'/interscale_Y.dat')
        interscale_Z_MATLAB = np.genfromtxt(datadir+'/interscale_Z.dat')
        interscale_T_MATLAB = np.genfromtxt(datadir+'/interscale_T.dat')
        interscale_wT_MATLAB = np.genfromtxt(datadir+'/interscale_wT.dat')

        fig1, axs1 = plt.subplots(2,3,figsize=(15,10))
        fig2, axs2 = plt.subplots(2,3,figsize=(15,10))
        fig3, axs3 = plt.subplots(2,3,figsize=(15,10))

        axs1[0,0].plot(z,turb_trans_X_MATLAB)
        axs1[0,0].plot(z,turb_trans_X_PO,'--')
        axs1[0,1].plot(z,turb_trans_Y_MATLAB)
        axs1[0,1].plot(z,turb_trans_Y_PO,'--')
        axs1[0,2].plot(z,turb_trans_Z_MATLAB)
        axs1[0,2].plot(z,turb_trans_Z_PO,'--')
        axs1[1,0].plot(z,turb_trans_T_MATLAB)
        axs1[1,0].plot(z,turb_trans_T_PO,'--')
        axs1[1,1].plot(z,turb_trans_wT_MATLAB)
        axs1[1,1].plot(z,turb_trans_wT_PO,'--')
        axs1[1,2].plot(z,0.5*(turb_trans_X_MATLAB+turb_trans_Y_MATLAB+turb_trans_Z_MATLAB))
        axs1[1,2].plot(z,stats[:,4],'--')
        axs1[0,0].set_title('X')
        axs1[0,1].set_title('Y')
        axs1[0,2].set_title('Z')
        axs1[1,0].set_title('T')
        axs1[1,1].set_title('wT')
        axs1[1,2].set_title('TKE')

        axs2[0,0].plot(z,turb_trans_mixed_X_MATLAB)
        axs2[0,0].plot(z,turb_trans_mixed_X_PO,'--')
        axs2[0,1].plot(z,turb_trans_mixed_Y_MATLAB)
        axs2[0,1].plot(z,turb_trans_mixed_Y_PO,'--')
        axs2[0,2].plot(z,turb_trans_mixed_Z_MATLAB)
        axs2[0,2].plot(z,turb_trans_mixed_Z_PO,'--')
        axs2[1,0].plot(z,turb_trans_mixed_T_MATLAB)
        axs2[1,0].plot(z,turb_trans_mixed_T_PO,'--')
        axs2[1,1].plot(z,turb_trans_mixed_wT_MATLAB)
        axs2[1,1].plot(z,turb_trans_mixed_wT_PO,'--')
        axs2[1,2].plot(z,0.5*(turb_trans_mixed_X_MATLAB+turb_trans_mixed_Y_MATLAB+turb_trans_mixed_Z_MATLAB))
        axs2[1,2].plot(z,stats[:,15],'--')
        axs2[0,0].set_title('X')
        axs2[0,1].set_title('Y')
        axs2[0,2].set_title('Z')
        axs2[1,0].set_title('T')
        axs2[1,1].set_title('wT')
        axs2[1,2].set_title('TKE')

        axs3[0,0].plot(z,interscale_X_MATLAB)
        axs3[0,0].plot(z,interscale_X_PO,'--')
        axs3[0,1].plot(z,interscale_Y_MATLAB)
        axs3[0,1].plot(z,interscale_Y_PO,'--')
        axs3[0,2].plot(z,interscale_Z_MATLAB)
        axs3[0,2].plot(z,interscale_Z_PO,'--')
        axs3[1,0].plot(z,interscale_T_MATLAB)
        axs3[1,0].plot(z,interscale_T_PO,'--')
        axs3[1,1].plot(z,interscale_wT_MATLAB)
        axs3[1,1].plot(z,interscale_wT_PO,'--')
        axs3[1,2].plot(z,0.5*(interscale_X_MATLAB+interscale_Y_MATLAB+interscale_Z_MATLAB))
        axs3[1,2].plot(z,stats[:,13],'--')
        axs3[0,0].set_title('X')
        axs3[0,1].set_title('Y')
        axs3[0,2].set_title('Z')
        axs3[1,0].set_title('T')
        axs3[1,1].set_title('wT')
        axs3[1,2].set_title('TKE')

        fig1.savefig('figures/turb_trans_scl'+str(scl+1).zfill(2)+'.png',dpi=600)
        fig2.savefig('figures/turb_trans_mixed_scl'+str(scl+1).zfill(2)+'.png',dpi=600)
        fig3.savefig('figures/interscale_scl'+str(scl+1).zfill(2)+'.png',dpi=600)

    if interscale_only:
        exit()
        
    # Now plot the error for the terms we have analytical fields for
    Nvec = np.array([64,128,256])
    terms = ['SGS tensor variance','SGS heat flux moments','Convective transport',\
            'Pressure transport','Molecular transport',\
            'SGS transport','SGS dissipation',\
            'Force production',\
            'Unsteady terms']
    
    for term in terms:
        errorR11r = []
        errorR22r = []
        errorR33r = []
        errorR11s = []
        errorR22s = []
        errorR33s = []
        errorTKEr = []
        errorTKEs = []
        errorTKEbr = []
        errorTKEbs = []
        errorwT1r = []
        errorwT1s = []
        errorwT2r = []
        errorwT2s = []
        errorTTr = []
        errorTTs = []
        print("Processing error for "+term,flush=True)
        for N in Nvec:
            readError = False
            scale = 'r'
            fname = 'test{}.log'.format(N)
            print("Reading "+fname)
            with open(fname,'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith(term):
                        #code.interact(local=locals())
                        readError = True
                        print("    > readError is True",flush=True)
                    if readError:
                        if line.startswith('> max absolute difference TKE_b'):
                            print("    > Reading TKE_b error",flush=True)
                            lsplt = line.split()
                            if scale == 'r':
                                errorTKEbr.append(float(lsplt[-1]))
                            elif scale == 's':
                                errorTKEbs.append(float(lsplt[-1]))
                        elif line.startswith('> max absolute difference TKE'):
                            print("    > Reading TKE error",flush=True)
                            lsplt = line.split()
                            if scale == 'r':
                                errorTKEr.append(float(lsplt[-1]))
                            elif scale == 's':
                                errorTKEs.append(float(lsplt[-1]))
                        elif line.startswith('> max absolute difference R11'):
                            lsplt = line.split()
                            print("    > Reading R11 error",flush=True)
                            if scale == 'r':
                                errorR11r.append(float(lsplt[-1]))
                            elif scale == 's':
                                errorR11s.append(float(lsplt[-1]))
                        elif line.startswith('> max absolute difference R22'):
                            print("    > Reading R22 error",flush=True)
                            lsplt = line.split()
                            if scale == 'r':
                                errorR22r.append(float(lsplt[-1]))
                            elif scale == 's':
                                errorR22s.append(float(lsplt[-1]))
                        elif line.startswith('> max absolute difference R33'):
                            print("    > Reading R33 error",flush=True)
                            lsplt = line.split()
                            if scale == 'r':
                                errorR33r.append(float(lsplt[-1]))
                            elif scale == 's':
                                errorR33s.append(float(lsplt[-1]))
                        elif line.startswith('> max absolute difference TT'):
                            print("    > Reading TT error",flush=True)
                            lsplt = line.split()
                            if scale == 'r':
                                errorTTr.append(float(lsplt[-1]))
                            elif scale == 's':
                                errorTTs.append(float(lsplt[-1]))
                        elif line.startswith('> max absolute difference wT1'):
                            print("    > Reading wT1 error",flush=True)
                            lsplt = line.split()
                            if scale == 'r':
                                errorwT1r.append(float(lsplt[-1]))
                            elif scale == 's':
                                errorwT1s.append(float(lsplt[-1]))
                        elif line.startswith('> max absolute difference wT2'):
                            print("    > Reading wT2 error",flush=True)
                            lsplt = line.split()
                            if scale == 'r':
                                errorwT2r.append(float(lsplt[-1]))
                            elif scale == 's':
                                errorwT2s.append(float(lsplt[-1]))
                        elif line.startswith('> max absolute difference Var(tau13s)'):
                            print("    > Reading Var(tau13s) error",flush=True)
                            lsplt = line.split()
                            errorR11s.append(float(lsplt[-1]))
                        elif line.startswith('> max absolute difference Var(tau23s'):
                            print("    > Reading Var(tau23s) error",flush=True)
                            lsplt = line.split()
                            errorR22s.append(float(lsplt[-1]))
                        elif line.startswith('> max absolute difference Var(q3s)'):
                            print("    > Reading Var(q3s) error",flush=True)
                            lsplt = line.split()
                            errorR11s.append(float(lsplt[-1]))
                        elif line.startswith('> max absolute difference for <fifi>,r'):
                            print("    > Reading <fifi>,r error:",flush=True)
                            lsplt = line.split()
                            errorR11r.append(float(lsplt[-1]))
                        elif line.startswith('> max absolute difference for <fifi>,s'):
                            print("    > Reading <fifi>,s erro:",flush=True)
                            lsplt = line.split()
                            errorR11s.append(float(lsplt[-1]))
                        elif line.startswith('> Scale 2'):
                            scale = 's'
                        elif line.startswith('*** *** ***'):
                            readError = False


        fig, axs = plt.subplots(1,2,figsize=(10,5))
        if len(errorR11r) == len(Nvec):
            axs[0].plot(Nvec,errorR11r,'o-',label='R11')
            print("Error for R11r: {}".format(errorR11r),flush=True)
        if len(errorR22r) == len(Nvec):
            axs[0].plot(Nvec,errorR22r,'x-',label='R22')
            print("Error for R22r: {}".format(errorR22r),flush=True)
        if len(errorR33r) == len(Nvec):
            axs[0].plot(Nvec,errorR33r,'d--',label='R33')
            print("Error for R33r: {}".format(errorR33r),flush=True)
        if len(errorTKEr) == len(Nvec):
            axs[0].plot(Nvec,errorTKEr,'s-',label='TKE')
            print("Error for TKEr: {}".format(errorTKEr),flush=True)
        if len(errorTKEbr) == len(Nvec):
            axs[0].plot(Nvec,errorTKEbr,'^-',label='TKEb')
            print("Error for TKEbr: {}".format(errorTKEbr),flush=True)
        if len(errorTTr) == len(Nvec):
            axs[0].plot(Nvec,errorTTr,'^-',label='TT')
            print("Error for TTr: {}".format(errorTTr),flush=True)
        if len(errorwT1r) == len(Nvec):
            axs[0].plot(Nvec,errorwT1r,'^-',label='wT1')
            print("Error for wT1r: {}".format(errorwT1r),flush=True)
        if len(errorwT2r) == len(Nvec):
            axs[0].plot(Nvec,errorwT2r,'^-',label='wT2')
            print("Error for wT2r: {}".format(errorwT2r),flush=True)
        print("Finished plotting Scale 1",flush=True)
        axs[0].set_title('Scale 1')
        if len(errorR11s) == len(Nvec):
            axs[1].plot(Nvec,errorR11s,'o-',label='R11')
            print("Error for R11s: {}".format(errorR11s),flush=True)
        if len(errorR22s) == len(Nvec):
            axs[1].plot(Nvec,errorR22s,'x-',label='R22')
            print("Error for R22s: {}".format(errorR22s),flush=True)
        if len(errorR33s) == len(Nvec):
            axs[1].plot(Nvec,errorR33s,'d--',label='R33')
            print("Error for R33s: {}".format(errorR33s),flush=True)
        if len(errorTKEs) == len(Nvec):
            axs[1].plot(Nvec,errorTKEs,'s-',label='TKE')
            print("Error for TKEs: {}".format(errorTKEs),flush=True)
        if len(errorTKEbs) == len(Nvec):
            axs[1].plot(Nvec,errorTKEbs,'^-',label='TKEb')
            print("Error for TKEbs: {}".format(errorTKEbs),flush=True)
        if len(errorTTs) == len(Nvec):
            axs[1].plot(Nvec,errorTTs,'^-',label='TT')
            print("Error for TTs: {}".format(errorTTs),flush=True)
        if len(errorwT1s) == len(Nvec):
            axs[1].plot(Nvec,errorwT1s,'^-',label='wT1')
            print("Error for wT1s: {}".format(errorwT1s),flush=True)
        if len(errorwT2s) == len(Nvec):
            axs[1].plot(Nvec,errorwT2s,'^-',label='wT2')
            print("Error for wT2s: {}".format(errorwT2s),flush=True)
        print("Finished plotting Scale 2",flush=True)
        axs[1].set_title('Scale 2')
        axs[1].legend()
        for i in range(2):
            if term == 'Unsteady terms':
                axs[i].plot(Nvec,1.0e3*Nvec**(-2.0),'k')
            else:
                axs[i].plot(Nvec,1.0e5*Nvec**(-6.0),'k')
            axs[i].set_yscale('log')
            axs[i].set_xscale('log')
            axs[i].grid()
        fig.tight_layout()
        fig.suptitle(term)
        np.savez('data/'+term.split()[0]+'_'+term.split()[1]+'.npz',errorR11r=errorR11r,\
                Nvec=Nvec)
        fig.savefig('figures/'+term.split()[0]+'_'+term.split()[1]+'.png',dpi=600)
    #plt.show()
