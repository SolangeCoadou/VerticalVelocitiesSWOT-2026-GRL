from scipy import signal, ndimage, interpolate, stats
import pyinterp
import pyinterp.fill
import seawater as sw
from gsw import f
import xarray as xr
import numpy as np
from scipy.interpolate import griddata

#######################################################################
### Functions for processing the data (interpolating, filtering...) ###
#######################################################################
# Some of the following code has been created by Elisa Carli, 
# Lia Siegelman and Patrice Klein.
# https://doi.org/10.5281/zenodo.15088480

def interp_to_prof(data, z1, z2):
    """
    1D interpolation.

    Parameters:
    ----------
    data : array_like (M)
        1D array of the data to be interpolated.
    z1 : array_like (M)
        1D array of the original coordinates.
    z2 : array_like (M)
        1D array of the coordinates at which to interpolate
        the data.
    Returns:
    --------
       : array_like (M)
        1D array of the value of data along z2
    """
    y_interp = interpolate.interp1d(z1, data)
    return y_interp(z2)

def interp_to_grid(data, x1, y1, x2, y2):
    """
    2D interpolation.

    Parameters:
    ----------
    data : array_like (M x N)
        2D array of the data to be interpolated.
    x1 : array_like (M x N)
        2D array of the original coordinates in the longitude direction.
    y1 : array_like (M x N)
         2D array of the original coordinates in the latitude direction.
    x2 : array_like (M x N)
        2D array of the coordinates in the longitude direction at which to interpolate the data.
    y2 : array_like (M x N)
        2D array of the coordinates in the latitude direction at which to interpolate the data.
    Returns:
    --------
       : array_like (M x N)
        2D array of the value of data at grid points (x2,y2)
    """
    return griddata((x1.ravel(), y1.ravel()), data.ravel(), (x2, y2), method='linear')

def LanczosKernel(cutoff_frequency, hws):
    """Lanczos 2D kernel.
    
    Parameters
    ----------
    hws : int
        Half window size in indexes
    cutoff_frequency : float
        Cutoff frequency in inverse abscissas steps
    Returns:
    --------
    kernel : array_like (hws*2+1 x hws*2+1)
        Lanczos 2D kernel
    """
    x = np.arange(-hws, hws + 1)
    sinc_filter = 2 * cutoff_frequency * np.sinc(2 * cutoff_frequency * x)
    sigma_factor = np.sinc(x / hws)
    
    kernel_x = sinc_filter * sigma_factor
    kernel_y = sinc_filter * sigma_factor
    kernel = np.dot(
        np.expand_dims(kernel_x, 1),
        np.expand_dims(kernel_y, 1).T)
    return kernel
    
def filter_convolution2d(data, cutoff_length,mode='reflect', cval=0.0):
    """Lanczos 2D filter.
    
    Parameters
    ----------
    data : array_like (M x N)
        2D array of the data to be filtered.
    cutoff_length : float
        Cutoff length scale in km
    Returns:
    --------
    result : array_like (M x N)
        Filtered data using a 2D Lanczos filter
    """
    
    dx=2
    hws=2*(cutoff_length/dx)-1
    cutoff_frequency = 1/(cutoff_length/dx)
    kernel = LanczosKernel(cutoff_frequency, hws)
    data=data.copy()
    mask = np.isnan(data)
    data[mask] = 0
    z_filtered = ndimage.convolve(data, kernel, mode=mode, cval=cval)

    # Weight correction for ignoring the nan
    w = np.ones_like(data)
    w[mask] = 0
    correction = ndimage.convolve(w, kernel, mode=mode, cval=cval)

    # Weight correction application and masking
    result = z_filtered / correction
    result[mask] = np.nan
    
    return result

def filtering_lanczos(dataset, cutoff_length):
    """Apply a Lanczos 2D filter to a DataArray.
    
    Parameters
    ----------
    dataset : DataArray_like
        DataArray to be filtered along dimensions x and y.
    cutoff_length : float
        Cutoff length scale in km
    Returns:
    --------
    data_final : DataArray_like
        Filtered DataArray reshape in the x and y dimensions to a square
    """
    data_filter = xr.apply_ufunc(
        filter_convolution2d,# first the function
        dataset,cutoff_length,
        input_core_dims=[["y","x"],[]],
        output_core_dims=[["y","x"]],
        vectorize=True,
        dask="parallelized",
        output_dtypes=['float64']
    )

    #Resize DataArray to be a square
    size_x=len(data_filter.x)
    size_y=len(data_filter.y)

    if size_y>size_x:
        data_final = data_filter.isel(y=slice(0,size_x))
    else:
        data_final = data_filter.isel(x=slice(0,size_y))
    return data_final


##############################
### eSQG general functions ###
##############################
# The following code has been created by Elisa Carli, 
# Lia Siegelman and Patrice Klein.
# https://doi.org/10.5281/zenodo.15088480

def doubly_per(data):
    """ Make data doubly periodic by mirror symmetry
    
    Parameters
    ----------
    data : array_like (M x N)
        2D array of the data to be transformed.
    Returns:
    --------
    data_per : array_like (2*M x 2*N)
        Data which has been doubly periodized
    """
    data = np.hstack((data,np.fliplr(data)))
    data_per = np.vstack((np.flipud(data),data))
    return data_per

def ssh_preprocessing_sqg(ssh):
    """ Apply doubly_per and remove the mean
    
    Parameters
    ----------
    ssh : array_like (M x N)
        2D array of the data to be transformed.
    Returns:
    --------
    ssh_det1_nomean : array_like (2*M x 2*N)
        Data which has been doubly periodized and from which
        the mean has been removed
    """
    
    ssh_doub_periodic = doubly_per(ssh)
    ssh_det1_nomean = ssh_doub_periodic - np.mean(ssh_doub_periodic)

    return ssh_det1_nomean

def get_kxky(L,M,dx,dy):
    """ Compute the horizontal wavenumbers vectors, matrix and norm
    
    Parameters
    ----------
    L : int
      number of grid points along the x axis
    M : int
      number of grid points along the y axis
    dx : float
      grid spacing in the x direction in m
    dy : int
      grid spacing in the y direction in m
    Returns:
    --------
    kx : array_like (L)
      1D array of the wavenumber vector in the x direction
    ky : array_like (M)
      1D array of the wavenumber vector in the y direction
    kxx : array_like (L x M)
      2D array of the wavenumber vector in the x direction
    kyy : array_like (L x M)
      2D array of the wavenumber vector in the y direction
    kk : array_like (L x M)
      2D array of the norm of the wavenumber
    """
    Lx = L*dx
    Ly = M*dy
    coefx = Lx / (2*np.pi)
    coefy = Ly / (2*np.pi)

    kx = np.hstack((np.arange(L/2),0,np.arange(-L/2+1,0)))/coefx
    ky = np.hstack((np.arange(M/2),0,np.arange(-M/2+1,0)))/coefy

    kkx,kky = np.meshgrid(kx,ky)
    kk = np.sqrt(kkx**2 + kky**2)

    return kx,ky,kkx,kky,kk

def sqg_rel_vort(ssh, M, L, N, kk, zz, gof0):
    
    '''
    M and L are the size of the double periodic field
    '''

    spec_vortssh = - gof0*(kk**2)*np.fft.fft2(ssh)

    spec_vortssh_z=np.zeros((M,L))
    spec_vortssh1=spec_vortssh

    kz = 0
    for k in range(0,N):
        if kz == 0:
            vort_sqg_matr = np.zeros([len(np.arange(0,N)),int(M/2),int(L/2)]) 

        #print(r'k = '+str(k))

        func_exp = np.zeros((M, L))
        func_exp = np.exp(kk*zz[k])
        spec_vortssh_z=spec_vortssh1*func_exp
        xis_sqg=np.flipud(np.fft.ifft2(spec_vortssh_z).real[:int(np.ceil((M/2))), :int(np.ceil((L/2)))])
        vort_sqg_matr[k,:,:] = xis_sqg
        kz = kz+1
        
    return vort_sqg_matr

def sqg_w(ssh, cc, rho_cst, M, L, N, kk, kkx, kky, zz, goN02orho0, gof0):
    
    '''
    M and L are the size of the double periodic field
    '''
    
    spec_u_ssh=-1j*gof0*kky*np.fft.fft2(ssh)
    spec_v_ssh= 1j*gof0*kkx*np.fft.fft2(ssh)
    spec_rho_ssh=-(kk*np.fft.fft2(ssh))*rho_cst

    u_ssh=np.fft.ifft2(spec_u_ssh).real
    v_ssh=np.fft.ifft2(spec_v_ssh).real
    rho_ssh=np.fft.ifft2(spec_rho_ssh).real

    jac1s=u_ssh*rho_ssh
    jac2s=v_ssh*rho_ssh

    spec_jac1= 1j*kkx*np.fft.fft2(jac1s)
    spec_jac2= 1j*kky*np.fft.fft2(jac2s)
    spec_jacs= spec_jac1+spec_jac2


    kz=0
    for k in range(0,N):  

        if kz == 0:
            w_sqg_matr = np.zeros([len(np.arange(0,N)),int(M/2),int(L/2)])

        func_exp = np.zeros((M, L))
        func_exp = np.exp(kk*zz[k])

        v_z= np.fft.ifft2(spec_v_ssh*func_exp).real
        u_z= np.fft.ifft2(spec_u_ssh*func_exp).real
        rho_z= np.fft.ifft2(spec_rho_ssh*func_exp).real

        jac1z=u_z*rho_z
        jac2z=v_z*rho_z
        spec_jac1= 1j*kkx*np.fft.fft2(jac1z)
        spec_jac2= 1j*kky*np.fft.fft2(jac2z)
        spec_jacz= spec_jac1+spec_jac2
        jacz=np.fft.ifft2(spec_jacz).real

        spec_jacsz= spec_jacs*func_exp
        jacs=np.fft.ifft2(spec_jacsz).real

        w_sqg=-jacs + jacz
        w_sqg= goN02orho0*w_sqg*cc
        w_sqg=86400*np.flipud(w_sqg[:int(np.ceil((M/2))), :int(np.ceil((L/2)))])
        w_sqg_matr[kz,:,:] = w_sqg

        kz = kz+1
        
    return(w_sqg_matr)

def sqg_strain(ssh, M, L, N, kk, kkx, kky, zz, gof0):
    
    '''
    M and L are the size of the double periodic field
    '''
    
    spec_sn = 2*kkx*kky*gof0*np.fft.fft2(ssh)
    spec_ss = (kky**2 - kkx**2)*gof0*np.fft.fft2(ssh)

    kz = 0
    for k in range(0,N):
        if kz == 0:
            vort_sqg_matr = np.zeros([len(np.arange(0,N)),int(M/2),int(L/2)]) 

        print(r'k = '+str(k))

        func_exp = np.zeros((M, L))
        func_exp = np.exp(kk*zz[k])
        
        spec_sn_z = spec_sn * func_exp
        spec_ss_z = spec_ss * func_exp
        
        sn = np.fft.ifft2(spec_sn_z).real
        ss = np.fft.ifft2(spec_ss_z).real
        
        strain = np.sqrt(sn**2 + ss**2)

        xis_sqg=np.flipud(strain[:int(np.ceil((M/2))), :int(np.ceil((L/2)))])
        vort_sqg_matr[k,:,:] = xis_sqg

        kz = kz+1
        
    return vort_sqg_matr

####################
### eSQG INATL60 ###
####################

def find_scientific_notation(number):
    base10 = np.log10(abs(number))
    return abs(np.floor(base10))

def compute_SQG_zeta(data,N02,f0,z_model, dx):
    ssh_sqg = ssh_preprocessing_sqg(data) 
    N0=np.sqrt(N02)        
    g = 9.81
    N0of0=N0/f0
    gof0=g/f0
    zz= z_model*N0of0
    N=len(zz)
    M=len(ssh_sqg[:,0])
    L=len(ssh_sqg[0,:])
            
    kx,ky,kkx,kky,kk = get_kxky(L,M,dx*1000,dx*1000)
    zeta = sqg_rel_vort(ssh_sqg, M, L, N, kk, zz, gof0)
    return zeta

def compute_SQG_w(data,N02,cc,f0,z_model, dx):
    ssh_sqg = ssh_preprocessing_sqg(data) 
    N0=np.sqrt(N02)        
    ## calcul de zz_rmean= z_rmean*N/f ##
    g = 9.81
    N0of0=N0/f0
    rho0=1.e+3
    gof0=g/f0
    goN02orho0=g/N02/rho0
    rho_cst=N0of0*rho0
    zz= z_model*N0of0
    N=len(zz)
    M=len(ssh_sqg[:,0])
    L=len(ssh_sqg[0,:])
            
    kx,ky,kkx,kky,kk = get_kxky(L,M,dx*1000,dx*1000)
    wref = sqg_w(ssh_sqg, cc, rho_cst, M, L, N, kk, kkx, kky, zz, goN02orho0, gof0)
    return wref

####################
### eSQG on SWOT ###
####################

def process_SSH_before_computing_derivatives(file_id, reduce=True):
    
    latmin, latmax = np.min(latBIG)-.2, np.max(latBIG)+.2
    lonmin, lonmax = np.min(lonBIG)-.2, np.max(lonBIG)+.2
    pos_west,pos_east,pos_south,pos_north  = lonmin,lonmax,latmin,latmax

    #Cut data around the region of interest
    swot = xr.open_dataset(file_id)
    id_select = np.where((swot.longitude>pos_west) & (swot.longitude<pos_east) & (swot.latitude>pos_south) & (swot.latitude<pos_north))[0]
    swot = swot.isel(num_lines=slice(id_select[0], id_select[-1]))

    nadir_pos=np.arange(32,37)

    ssha = swot['ssha_filtered'].data
    for n in nadir_pos:
        if len(np.where(np.isnan(ssha[:,n])==0)[0])>0:
            ssha[:,n]=np.nan
    swot['ssha_filtered'] = xr.DataArray(ssha,dims=('num_lines', 'num_pixels'),coords={'latitude':swot.latitude, 'longitude':swot.longitude})

    ssh = np.squeeze(swot.ssha_filtered.data.reshape(-1,1)+swot.mdt.data.reshape(-1,1))

    #Mask data around NaNs
    thresh_nan=.5/9

    ssha_complex = np.zeros_like(ssha, dtype=np.complex64)
    ssha_complex[np.isnan(ssha)] = np.array((1j))
    ssha_complex[np.bitwise_not(np.isnan(ssha))] = ssha[np.bitwise_not(np.isnan(ssha))]

    k = np.array([
        [1, 1, 1],
        [1, 1, 1],
        [1, 1, 1]])
    k = k/k.sum() 

    convolution = signal.convolve2d(ssha_complex, k, 'same', 'wrap')

    col_all_nan = [np.isnan(swot['ssha_filtered'].data[:,c]).all() for c in range(len(swot['ssha_filtered'].data[0,:]))]
    count_m=0
    for cc in range(len(col_all_nan)-1):
        if np.diff(col_all_nan)[cc]:
            if count_m%2==0:   
                idx=cc+1
            else:
                idx=cc
            count_m+=1
            img_value,count_v = np.unique(np.imag(convolution[:,idx]),return_counts=True)
            most_common = img_value[np.argmax(count_v)]
            convolution[np.where(np.imag(convolution[:,idx])==most_common)[0],idx]=1

    idnotnan_conv = np.where(np.imag(convolution.reshape(-1,1))>thresh_nan)[0]
    ssh[idnotnan_conv]=np.nan

    #Filter
    idnotnan = np.where(np.isnan(ssh)==0)[0]
    if len(idnotnan)>0:
        y_axis = pyinterp.Axis(swot.num_lines.values)
        x_axis = pyinterp.Axis(swot.num_pixels.values, is_circle=False)
        grid = pyinterp.Grid2D(x_axis, y_axis, ssh.reshape(swot.longitude.data.shape).T)
        has_converged, filled = pyinterp.fill.gauss_seidel(grid)
        swot['ssh_interp'] = (('num_lines', 'num_pixels'), filled.T)
 
    else:
        ssh_interp=np.zeros_like(swot.longitude.data)+np.nan
        swot['ssh_interp'] = xr.DataArray(ssh_interp,dims=('num_lines', 'num_pixels'),coords={'latitude':swot.latitude, 'longitude':swot.longitude})
    swot['mask_nan'] = xr.DataArray(np.isnan(swot.ssha_filtered),dims=('num_lines', 'num_pixels'),coords={'latitude':swot.latitude, 'longitude':swot.longitude})

    return swot

def compute_sqg_swot(depth, file_id, field):
    swot_processed = process_SSH_before_computing_derivatives(file_id)
    
    #Remove lines and columns where some nans remain
    #First columns
    id_nan_col = np.where(np.isnan(swot_processed.ssh_interp)==1)[1]
    col, nb_time = np.unique(id_nan_col, return_counts=True)
    idsel = np.where(nb_time>1)[0]
    col = col[idsel]
    
    col_all = np.arange(0,len(swot_processed.ssh_interp[0,:]), dtype='float64')
    for c in col:
       col_all[c]=np.nan
    col_all = col_all[np.where(np.isnan(col_all)==0)[0]]
    swot_processed = swot_processed.isel(num_pixels=col_all.astype(int))
    if len(swot_processed.ssh_interp[0,:])==0:
         wref, mask_ref, lon_ref, lat_ref = [np.nan], [np.nan], [np.nan], [np.nan]
    else:
    
        #Then lines
        id_nan_col = np.where(np.isnan(swot_processed.ssh_interp)==1)[0]
        col, nb_time = np.unique(id_nan_col, return_counts=True)
        idsel = np.where(nb_time>1)[0]
        col = col[idsel]
            
        col_all = np.arange(0,len(swot_processed.ssh_interp[:,0]), dtype='float64')
        for c in col:
            col_all[c]=np.nan
        col_all = col_all[np.where(np.isnan(col_all)==0)[0]]
        swot_processed = swot_processed.isel(num_lines=col_all.astype(int))
        
        if len(swot_processed.ssh_interp[:,0])==0:
            wref, mask_ref, lon_ref, lat_ref = [np.nan], [np.nan], [np.nan], [np.nan]
        else:
            #Check that we didn't remove a column or a line in the middle of the swath, thus breaking continuity
            dist_lines = sw.dist(swot_processed.latitude[0,:],swot_processed.longitude[0,:],units='km')[0]
            dist_column = sw.dist(swot_processed.latitude[:,0],swot_processed.longitude[:,0],units='km')[0]
            
            if len(np.where(np.concatenate((dist_lines,dist_column))>3)[0])>0:
                wref, mask_ref, lon_ref, lat_ref = [np.nan], [np.nan], [np.nan], [np.nan]
            else:
                #Start with the filtering
                dx = 2
                cutoff_lenght = 40
                hws=int(2*(cutoff_lenght/dx)-1)
                cutoff_frequency = 1/(cutoff_lenght/dx)
        
                lon = swot_processed.longitude.data
                lat = swot_processed.latitude.data
                
                field_filter_ref = filter_convolution2d(LanczosKernel(cutoff_frequency, hws),swot_processed.ssh_interp.data)#0 - 60
                lon_ref=lon[6:-6,6:-6]
                lat_ref=lat[6:-6,6:-6]
                mask_ref = swot_processed.mask_nan.data[6:-6,6:-6]
                    
                #Compute w
                ssh_det2 = ssh_preprocessing_sqg(field_filter_ref)
                
                #Set eSQG parameters
                N02 = 3.4e-6
                N0=np.sqrt(N02)
                cc=2
        
                g = 9.81
                f0 = abs(np.mean(f(swot_processed.latitude.data.flatten())))
                N0of0=N0/f0
                rho0=1.e+3
                gof0=g/f0
                Rho0oN02og=rho0/N02/g
                goN02orho0=g/N02/rho0
                rho_cst=N0of0*rho0
                M=len(ssh_det2[:,0])
                L=len(ssh_det2[0,:])
                kx,ky,kkx,kky,kk = get_kxky(L,M,2000,2000)
                if len(depth_list)==1:
                    zz = np.array([-depth])*N0of0
                else:
                    zz = np.array(np.squeeze([-depth]))*N0of0
                N=len(zz)
                
                if field=='vorticity':
                    wref = np.squeeze(sqg_rel_vort(ssh_det2, M, L, N, kk, zz, gof0))/f0
                elif field=='velocity':
                    wref = np.squeeze(sqg_w(ssh_det2, cc, rho_cst, M, L, N, kk, kkx, kky, zz, goN02orho0, gof0))
                elif field=='strain':
                    wref = np.squeeze(sqg_strain(ssh_det2, M, L, N, kk, kkx, kky, zz, gof0))/f0
                if len(wref.shape)==2:
                    wref=wref[6:-6,6:-6]
                else:
                    wref=wref[:,6:-6,6:-6]
    return wref, mask_ref, lon_ref, lat_ref

#########################################
### Functions for spectra computation ###
#########################################
# Coded by Hector Torres (NASA-JPL)
# Adapted by Felix Vivant (Scripps-UCSD and ENS Paris-Saclay)
# https://doi.org/10.5281/zenodo.13923050

def spectra_w(A,B,d1,d2):
    ## Wavenumber spectral analysis
    ### Outputs:
    ## A,B,cs,coh,f1,f2,df1,df2,Atofft,Btofft
    ##
    ## where
    ## A is the power spectrum of A
    ## B is the power spectrum of B
    ## cs is the cospectrum
    ## coh is the spectral coherence
    ## f1,f2 are the spectral coordinates
    ## df1,df2 are the spectral resolution
    ## Atofft,Btofft are the variable in the physical space 
    ## used for the FFt, i.e. after detrending and hanning windowing
    
    import numpy as np
    l1,l2 = A.shape
    df1 = 1./(l1*d1)
    df2 = 1./(l2*d2)
    f1Ny = 1./(2*d1)
    
    f1 = np.arange(-f1Ny,f1Ny,df1)
    #f2 = np.arange(0,l2/2+1)*df2
    if l2%2==1:
        f2 = np.arange(0,l2/2)*df2 
    elif l2%2==0:
        f2 = np.arange(0,l2/2+1)*df2 
    # Spectral window
    # the spatial window
    w1 = np.hanning(l1)
    wx = np.matrix(w1)
    w2 = np.hanning(l2)
    wy = np.matrix(w2)
    window_s = np.array(wx.T*wy).reshape(l1,l2)

    # ===== Spectral space ======

    corr_fact=d1*d2
    Atofft=window_s*A
    Btofft=window_s*B
    Ahat = np.fft.rfftn(Atofft)*corr_fact
    Bhat = np.fft.rfftn(Btofft)*corr_fact
    
    # Cospectrum of A and B
    cs = (Ahat*Bhat.conjugate()).real
    # cs = (Ahat*Ahat.conjugate()).real
    # Power spectrum of A
    A = (Ahat*Ahat.conjugate()).real
    # Power spectrum of B
    B = (Bhat*Bhat.conjugate()).real
    #:::::
    ### Coherence #####
    coh = cs/(np.sqrt(A)*np.sqrt(B))
  
    # zero-padding
    coh = np.fft.fftshift(coh,axes=(0))
    # cospectrum density
    cs = np.fft.fftshift(cs,axes=(0))

    # power spectrum density
    A = np.fft.fftshift(A,axes=(0))
    B = np.fft.fftshift(B,axes=(0))
    return A,B,cs,coh,f1,f2,df1,df2,Atofft,Btofft

def calc_ispec(k,l,E):
    """ calculates isotropic spectrum from 2D spectrum """

    dk,dl = k[1,]-k[0],l[1]-l[0]
    l,k = np.meshgrid(l,k)
    wv = np.sqrt(k**2 + l**2)

    if k.max()>l.max():
        kmax = l.max()
    else:
        kmax = k.max()

    # create radial wavenumber
    dkr = np.sqrt(dk**2 + dl**2)
    kr =  np.arange(dkr/2.,kmax+dkr,dkr)
    ispec = np.zeros(kr.size)

    for i in range(kr.size):
        fkr =  (wv>=kr[i]-dkr/2) & (wv<=kr[i]+dkr/2)
        # infinitesimal phase space angle in polar coordinates
        # dth = 2pi/fkr with fkr the discrete number of points 
        # in the disk bounded by two circle of radius kr-dkr/2 qnd kr+dkr/2
        # 2D case (x,y) of a w-omega spectrum with kx or ky not truncated for fft algorithm
        # fkr is the number over all the phase space
        # i.e. dth = pi / (fkr.sum())
        dth = np.pi / (fkr.sum())
        ispec[i] = E[fkr].sum() * kr[i] * dth *dkr

    return kr, ispec, dkr

def compute_coherence(data1, data2, dx):
    if len(np.where(np.isnan(data2.reshape(-1,1))==1)[0])>0:
        x=np.arange(len(data2[0,:]))
        y=np.arange(len(data2[0,:]))
        x, y=np.meshgrid(x, y)
        x=np.squeeze(x.reshape(-1,1))
        y=np.squeeze(y.reshape(-1,1))
        idnotnan=np.where(np.isnan(data2.reshape(-1,1))==0)[0]
        data_interp = griddata((x[idnotnan], y[idnotnan]), np.squeeze(data2.reshape(-1,1))[idnotnan], (x, y), method='cubic').reshape(data2.shape)
        B=data_interp
    else:
        B=data2
    A = signal.detrend(data1,axis=0,type='linear')
    A = signal.detrend(A,axis=1,type='linear')
    B = signal.detrend(B,axis=0,type='linear')
    B = signal.detrend(B,axis=1,type='linear')
    EuA,EuB,cs,coh,k,l,d1,d2,Atofft,Btofft = spectra_w(A,B,dx,dx)
    #Isospectra
    kiso,EiA,dkiso = calc_ispec(k,l,EuA[:,:])
    kiso,EiB,dkiso = calc_ispec(k,l,EuB[:,:])
    kiso,Ei_cs,dkiso = calc_ispec(k,l,cs[:,:])
    return kiso, np.abs(Ei_cs)**2 / EiA / EiB

def compute_correlation(data1, data2):
    print(f"data: {data1.shape}")
    a=data1.ravel()
    b=data2.ravel()
    idnotnan=np.where((np.isnan(a+b)==0) & (np.isinf(a+b)==0))[0]
    if len(a[idnotnan])>2:
        return stats.pearsonr(a[idnotnan], b[idnotnan])[0]
    else:
        return np.nan
