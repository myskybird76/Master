import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.colors import LogNorm
from scipy import ndimage
from scipy import stats
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from astropy.wcs import WCS
import cv2
from numpy import unravel_index
import time
import os
import glob
from openpyxl import Workbook
from astropy.visualization import hist
import seaborn as sns
import matplotlib.patches as mpatches
import sys

write_wb = Workbook()
write_ws = write_wb.active
write_ws.title = 'Result1'
write_ws.append(
    ['Name', 'w_T', 'w_+', 'w_-', 'c_T', 'c_+', 'c_-', 'G_T', 'G_+', 'G_-',
     'M20_T', 'M20_+', 'M20_-', 'A180_T', 'A180_+', 'A180_-', 'A90_T', 'A90_+', 'A90_-',
     'P2_T', 'P2_+', 'P2_-', 'P3_T', 'P3_+', 'P3_-', 'P4_T', 'P4_+', 'P4_-'])

R500 = np.genfromtxt('D:/Master thesis/ForSangwoo(final)/r500(0.2).csv', dtype = None, delimiter=',', skip_header = 1, encoding = 'UTF8', names = ('Name', 'R500'))

os.chdir("D:/Master thesis/ForSangwoo(final)/0.2/")

rawimage2 = glob.glob("*_b8_ghostmos.fits")
backgrd2 = glob.glob('*_b8_ghostbgdmos.fits')
backgrd_sm2 = glob.glob('*_b8_ghostbgdmos_sm.fits')
expmap2 = glob.glob('*_b8_ghostexpmos.fits')
ClModel2 = glob.glob('*_b8_ghostmosCR_sm.fits')

for i in range(len(rawimage2)) :

    r500 = R500['R500'][i]
    Name = R500['Name'][i]
    print(Name)
    rawimage, header = fits.getdata(rawimage2[i], header = True)
    backgrd = fits.getdata(backgrd2[i])
    backgrd_sm = fits.getdata(backgrd_sm2[i])
    expmap = fits.getdata(expmap2[i])
    ClModel = fits.getdata(ClModel2[i])

    imagemodel = ClModel * expmap + backgrd_sm
    image = rawimage - backgrd

    kernel = Gaussian2DKernel(2) #################################################################################### kernel

    nx = header['NAXIS1']
    ny = header['NAXIS2']
    scale = header['CDELT2'] * np.float64(60.0e0)
    r500_pix = r500 / scale

    raw_data = rawimage - backgrd
    sim_data = imagemodel - backgrd_sm
    simu_raw = convolve(raw_data, kernel)
    simu_sim = convolve(sim_data, kernel)
    simu_exp = convolve(expmap, kernel)
    raw = simu_raw / simu_exp ################################################## important!
    sim = simu_sim / simu_exp
    i = (expmap == 0)
    raw[i] = 0
    sim[i] = 0

    center = ndimage.measurements.center_of_mass(sim)
    #ind = unravel_index(corim.argmax(), corim.shape)
    #xcen_peak = np.round(ind[1], 0)
    #ycen_peak = np.round(ind[0], 0)
    xcen_ini0 = np.round(center[1], 0)
    ycen_ini0 = np.round(center[0], 0)
    #print("Peak : ", xcen_peak, ycen_peak)
    #print("Mass Initial : ", xcen_ini, ycen_ini)

    def Center0(img, xcen, ycen, aper) :

        nx = np.size(img[0,])
        ny = np.size(img[:, 0])
        rep = np.repeat(np.float64(1.0e0), ny)
        rep1 = rep.reshape(ny, 1)
        rep2 = np.arange(ny)
        rep3 = rep2.reshape(ny, 1)
        offset = 100

        while offset > 0.01:
            xx = np.multiply(np.arange(nx), rep1) - xcen
            yy = np.multiply(np.repeat(np.float64(1.0e0), nx), rep3) - ycen
            Rad = np.sqrt(xx ** 2 + yy ** 2)
            ap = np.float64(aper)
            ap_ind = np.where(Rad <= ap)
            imared = img[ap_ind]

            centroid1 = float(np.round(np.sum(xx[ap_ind] * imared) / np.sum(imared) + xcen, 1))
            centroid2 = float(np.round(np.sum(yy[ap_ind] * imared) / np.sum(imared) + ycen, 1))

            offset = np.sqrt((xcen - centroid1) ** 2 + (ycen - centroid2) ** 2)
            xcen = centroid1
            ycen = centroid2

        ap_out = np.where(Rad > ap)
        img2 = np.copy(img)
        img2[ap_out] = 0
        ind = unravel_index(img2.argmax(), img2.shape)
        xcen_peak = np.round(ind[1], 0)
        ycen_peak = np.round(ind[0], 0)

        return xcen, ycen, xcen_peak, ycen_peak

    def Concentration0(img, xcen, ycen, aper, frac) :

        nx = np.size(img[0,])
        ny = np.size(img[:,0])
        rep = np.repeat(np.float64(1.0e0), ny)
        rep1 = rep.reshape(ny, 1)
        rep2 = np.arange(ny)
        rep3 = rep2.reshape(ny, 1)
        xx = np.multiply(np.arange(nx), rep1) - xcen
        yy = np.multiply(np.repeat(np.float64(1.0e0), nx), rep3) - ycen
        Rad = np.sqrt(xx**2 + yy**2)
        ap = np.float64(aper)
        ap_ind = np.where(Rad <= ap)
        Flx_large = np.sum(img[ap_ind])
        ap_ind2 = np.where(Rad <= ap*frac)
        Flx_small = np.sum(img[ap_ind2])

        return Flx_small / Flx_large

    def Centroid0(img, xcen, ycen, xpeak, ypeak, aper, num) :

        nx = np.size(img[0,])
        ny = np.size(img[:,0])
        rep = np.repeat(np.float64(1.0e0), ny)
        rep1 = rep.reshape(ny, 1)
        rep2 = np.arange(ny)
        rep3 = rep2.reshape(ny, 1)
        xx = np.multiply(np.arange(nx), rep1) - xcen
        yy = np.multiply(np.repeat(np.float64(1.0e0), nx), rep3) - ycen
        Rad = np.sqrt(xx**2 + yy**2)
        ap = np.float64(aper)
        ap_rad = np.linspace(ap, 0.0, num + 1)
        pdist_app = []
        cen1 = []
        cen2 = []

        for i in ap_rad :
            if i != 0 :
                ap_ind = np.where(Rad <= i)
                imared = img[ap_ind]
                centroid1 = float(np.round(np.sum(xx[ap_ind] * imared) / np.sum(imared) + xcen, 1))
                centroid2 = float(np.round(np.sum(yy[ap_ind] * imared) / np.sum(imared) + ycen, 1))
                cen1.append(centroid1 - xpeak)
                cen2.append(centroid2 - ypeak)

        mean1 = np.mean(cen1)
        mean2 = np.mean(cen2)

        for i in ap_rad :
            if i != 0 :
                ap_ind = np.where(Rad <= i)
                imared = img[ap_ind]
                centroid1 = float(np.round(np.sum(xx[ap_ind] * imared) / np.sum(imared) + xcen))
                centroid2 = float(np.round(np.sum(yy[ap_ind] * imared) / np.sum(imared) + ycen))
                pdist = (((centroid1 - xpeak) - mean1)**2 + ((centroid2 - ypeak) - mean2)**2)
                #pdist = np.sqrt((centroid1 - xcen)**2 + (centroid2 - ycen)**2)
                pdist_app.append(pdist)

        tot = np.sum(pdist_app)
        Centshift = np.sqrt(tot / (num - 1)) / (r500_pix) #np.std(pdist_app) / r500_pix

        return Centshift


    def Gini0(img, xcen, ycen, aper):

        nx = np.size(img[0,])
        ny = np.size(img[:, 0])
        rep = np.repeat(np.float64(1.0e0), ny)
        rep1 = rep.reshape(ny, 1)
        rep2 = np.arange(ny)
        rep3 = rep2.reshape(ny, 1)
        xx = np.multiply(np.arange(nx), rep1) - xcen
        yy = np.multiply(np.repeat(np.float64(1.0e0), nx), rep3) - ycen
        Rad = np.sqrt(xx ** 2 + yy ** 2)
        ap = np.float64(aper)
        ap_ind = np.where(Rad <= ap)

        n = len(img[ap_ind])
        Ksig2 = []
        Fsig2 = []
        Ki2 = []
        Fi2 = []
        Km = []
        Fm = []
        img2 = img[ap_ind]
        Rad2 = Rad[ap_ind] ** 2
        img_new = img2[np.argsort(img2)]
        Rad_new = Rad2[np.argsort(img2)]
        img_rev = img2[np.argsort(-img2)]
        Rad_rev = Rad2[np.argsort(-img2)]

        i = 1

        while i < n + 1:
            Ki = np.abs(img_new[i - 1])  # K_i
            Fi = np.abs(img_new[i - 1]) * Rad_new[i - 1]  # F_i
            Ksig = (2 * i - n - 1) * Ki
            Fsig = (2 * i - n - 1) * Fi

            Ki2.append(Ki)
            Fi2.append(Fi)
            Ksig2.append(Ksig)
            Fsig2.append(Fsig)

            i = i + 1

        Kdot = np.mean(Ki2)
        Fdot = np.mean(Fi2)
        Kisum = np.sum(Ki2)
        Fisum = np.sum(Fi2)

        j = 1

        while j < n + 1:

            K_m = np.abs(img_rev[j - 1])
            F_m = np.abs(img_rev[j - 1]) * Rad_rev[j - 1]
            Km.append(K_m)
            Fm.append(F_m)

            Kmtot = np.sum(Km)
            Fmtot = np.sum(Fm)

            if Kmtot <= 0.2 * Kisum:

                j = j + 1

            else:
                Km.pop()
                Fm.pop()
                break

        sigma_F = np.sum(Fm)
        Ksigmai = np.sum(Ksig2)
        Fsigmai = np.sum(Fsig2)
        Gin = (1 / (Kdot * n * (n - 1))) * Ksigmai
        Gin_m = (1 / (Fdot * n * (n - 1))) * Fsigmai
        M20 = np.log10(sigma_F / Fisum)

        return Gin, Gin_m, M20

    def Asymmetry0(img, xcen, ycen, aper) :
        img = img.astype(np.float64)
        nx = np.size(img[0,])
        ny = np.size(img[:,0])
        rep = np.repeat(np.float64(1.0e0), ny)
        rep1 = rep.reshape(ny, 1)
        rep2 = np.arange(ny)
        rep3 = rep2.reshape(ny, 1)
        xx = np.multiply(np.arange(nx), rep1) - xcen
        yy = np.multiply(np.repeat(np.float64(1.0e0), nx), rep3) - ycen
        Rad = np.sqrt(xx**2 + yy**2)
        ap = np.float64(aper)
        ap_ind = np.where(Rad <= ap)
        center = (xcen, ycen)
        angle90 = 90
        angle180 = 180
        scale = 1.0
        M90 = cv2.getRotationMatrix2D(center, angle90, scale)
        M = cv2.getRotationMatrix2D(center, angle180, scale)
        rotated90 = cv2.warpAffine(img, M90, (nx, ny))
        rotated180 = cv2.warpAffine(img, M, (nx, ny))

        S = np.sum(img[ap_ind])
        T = np.sum(np.abs(img[ap_ind] - rotated180[ap_ind]))
        T90 = np.sum(np.abs(img[ap_ind] - rotated90[ap_ind]))

        #plt.imshow(rotated180, cmap = 'viridis', origin = 'lower', interpolation = 'none', norm = LogNorm())  #viridis or spectral
        #plt.plot(xcen, ycen, '.', color = 'red')
        #plt.colorbar(orientation = 'vertical')
        #plt.xlim(130, 815)
        #plt.ylim(160, 855)
        return T / S, T90 / S

    def Power0(img, xcen, ycen, aper, minorder, maxorder) :

        nx = np.size(img[0,])
        ny = np.size(img[:,0])
        rep = np.repeat(np.float64(1.0e0), ny)
        rep1 = rep.reshape(ny, 1)
        rep2 = np.arange(ny)
        rep3 = rep2.reshape(ny, 1)
        xx = np.multiply(np.arange(nx), rep1) - xcen
        yy = np.multiply(np.repeat(np.float64(1.0e0), nx), rep3) - ycen
        Rad = np.sqrt(xx**2 + yy**2)
        ap = np.float64(aper)
        i = (xx != 0)
        yy[i] /= xx[i]
        phi = np.arctan(yy)
        j = (xx < 0)
        j_n = len(j)
        if j_n > 0 :
            phi[j] += np.pi
        k = (xx > 0) & (yy < 0)
        k_n = len(k)
        if k_n > 0 :
            phi[k] += 2 * np.pi
        l = (xx == 0) & (yy >= 0)
        l_n = len(l)
        if l_n > 0 :
            phi[l] = np.pi / 2
        q = (xx == 0) & (yy < 0)
        q_n = len(q)
        if q_n > 0 :
            phi[q] = (3 * np.pi) / 2

        ap_ind = np.where(Rad <= ap)
        imared = img[ap_ind]
        radred = Rad[ap_ind]
        phired = phi[ap_ind]
        num = np.arange(minorder, maxorder + 1)
        po = []
        a_0 = np.sum(imared)
        P_0 = (a_0 * np.log(ap))**2

        for m in num :

           a_m = np.sum(imared * radred**m * np.cos(m * phired))
           b_m = np.sum(imared * radred**m * np.sin(m * phired))
           P_m = (a_m**2 + b_m**2) / (2 * (m**2) * ap**(2*m))
           po.append(P_m)

        power = po / P_0

        return power

    Cent_raw_1 = Center0(raw, xcen_ini0, ycen_ini0, r500_pix)
    Cent_sim_1 = Center0(sim, xcen_ini0, ycen_ini0, r500_pix)
    #print("Center :", Cent_raw_2[0], Cent_raw_2[1])
    #print("Peak : ", Cent_raw_2[2], Cent_raw_2[3])
    xcen_mass_raw_1 = Cent_raw_1[0]
    ycen_mass_raw_1 = Cent_raw_1[1]
    xcen_peak_raw_1 = Cent_raw_1[2]
    ycen_peak_raw_1 = Cent_raw_1[3]
    xcen_mass_sim_1 = Cent_sim_1[0]
    ycen_mass_sim_1 = Cent_sim_1[1]
    xcen_peak_sim_1 = Cent_sim_1[2]
    ycen_peak_sim_1 = Cent_sim_1[3]

    Cen_raw_1 = Centroid0(raw, xcen_mass_raw_1, ycen_mass_raw_1, xcen_peak_raw_1, ycen_peak_raw_1, r500_pix, 10)
    Cen_sim_1 = Centroid0(sim, xcen_mass_sim_1, ycen_mass_sim_1, xcen_peak_sim_1, ycen_peak_sim_1, r500_pix, 10)
    #print("Centroid Shift : {:.3g}".format(Cen0))

    Con_raw_1 = Concentration0(raw, xcen_mass_raw_1, ycen_mass_raw_1, r500_pix, 0.1)
    Con_sim_1 = Concentration0(sim, xcen_mass_sim_1, ycen_mass_sim_1, r500_pix, 0.1)
    #print("Concentration : {:.3g}".format(Con_raw_1))

    Gin_raw_1 = Gini0(raw, xcen_mass_raw_1, ycen_mass_raw_1, r500_pix)
    Gin_sim_1 = Gini0(sim, xcen_mass_sim_1, ycen_mass_sim_1, r500_pix)
    #print("Gini Coefficient : {:.3g}".format(Gin_raw_1[0]))
    #print("Second Order : {:.3g}".format(Gin_raw_1[1]))
    #print("M20 : {:.4g}".format(Gin_raw_1[2]))

    Asy_raw_1 = Asymmetry0(raw, xcen_mass_raw_1, ycen_mass_raw_1, r500_pix)
    Asy_sim_1 = Asymmetry0(sim, xcen_mass_sim_1, ycen_mass_sim_1, r500_pix)
    #print("Asymmetry Parameter 90 : {:.3g}".format(Asy_raw_1[1]))
    #print("Asymmetry Parameter 180 : {:.3g}".format(Asy_raw_1[0]))

    Pow_raw_1 = Power0(raw, xcen_mass_raw_1, ycen_mass_raw_1, r500_pix, 1, 4)
    Pow_sim_1 = Power0(sim, xcen_mass_sim_1, ycen_mass_sim_1, r500_pix, 1, 4)
    #print("Power Ratio P2/P0 : {:.3g}".format(Pow_raw_1[1]))
    #print("Power Ratio P3/P0 : {:.3g}".format(Pow_raw_1[2]))
    #print("Power Ratio P4/P0 : {:.3g}".format(Pow_raw_1[3]))

    k = 0
    Cen1 = []
    Con1 = []
    Gin1 = []
    M201 = []
    Asy180 = []
    Asy90 = []
    Pow2 = []
    Pow3 = []
    Pow4 = []

    while k < 100 :

        poi_1 = np.random.poisson(np.abs(imagemodel)) ## first
        simu_1 = poi_1 - backgrd_sm
        simu_smo_1 = convolve(simu_1, kernel)
        simu_exp = convolve(expmap, kernel)
        sim_2_1 = simu_smo_1 / simu_exp
        i = (expmap == 0)
        sim_2_1[i] = 0

        center = ndimage.measurements.center_of_mass(sim_2_1)
        #ind = unravel_index(sim_2.argmax(), sim_2.shape)
        #xcen_peak = np.round(ind[1], 0)
        #ycen_peak = np.round(ind[0], 0)
        xcen_ini = np.round(center[1], 0)
        ycen_ini = np.round(center[0], 0)

        def Center(img, xcen, ycen, aper) :

            nx = np.size(img[0,])
            ny = np.size(img[:, 0])
            rep = np.repeat(np.float64(1.0e0), ny)
            rep1 = rep.reshape(ny, 1)
            rep2 = np.arange(ny)
            rep3 = rep2.reshape(ny, 1)
            offset = 100

            while offset > 0.01:
                xx = np.multiply(np.arange(nx), rep1) - xcen
                yy = np.multiply(np.repeat(np.float64(1.0e0), nx), rep3) - ycen
                Rad = np.sqrt(xx ** 2 + yy ** 2)
                ap = np.float64(aper)
                ap_ind = np.where(Rad <= ap)
                imared = img[ap_ind]

                centroid1 = float(np.round(np.sum(xx[ap_ind] * imared) / np.sum(imared) + xcen, 1))
                centroid2 = float(np.round(np.sum(yy[ap_ind] * imared) / np.sum(imared) + ycen, 1))

                offset = np.sqrt((xcen - centroid1) ** 2 + (ycen - centroid2) ** 2)
                xcen = centroid1
                ycen = centroid2
                #print(xcen, ycen)

            ap_out = np.where(Rad > ap)
            img2 = np.copy(img)
            img2[ap_out] = 0
            ind = unravel_index(img2.argmax(), img2.shape)
            xcen_peak = np.round(ind[1], 0)
            ycen_peak = np.round(ind[0], 0)

            return xcen, ycen, xcen_peak, ycen_peak


        def Centroid(img, xcen, ycen, xpeak, ypeak, aper, num):

            nx = np.size(img[0,])
            ny = np.size(img[:, 0])
            rep = np.repeat(np.float64(1.0e0), ny)
            rep1 = rep.reshape(ny, 1)
            rep2 = np.arange(ny)
            rep3 = rep2.reshape(ny, 1)
            xx = np.multiply(np.arange(nx), rep1) - xcen
            yy = np.multiply(np.repeat(np.float64(1.0e0), nx), rep3) - ycen
            Rad = np.sqrt(xx ** 2 + yy ** 2)
            ap = np.float64(aper)
            ap_rad = np.linspace(ap, 0.0, num + 1)
            pdist_app = []
            cen1 = []
            cen2 = []

            for i in ap_rad:
                if i != 0:
                    ap_ind = np.where(Rad <= i)
                    imared = img[ap_ind]
                    centroid1 = float(np.round(np.sum(xx[ap_ind] * imared) / np.sum(imared) + xcen, 1))
                    centroid2 = float(np.round(np.sum(yy[ap_ind] * imared) / np.sum(imared) + ycen, 1))
                    cen1.append(centroid1 - xpeak)
                    cen2.append(centroid2 - ypeak)

            mean1 = np.mean(cen1)
            mean2 = np.mean(cen2)

            for i in ap_rad:
                if i != 0:
                    ap_ind = np.where(Rad <= i)
                    imared = img[ap_ind]
                    centroid1 = float(np.round(np.sum(xx[ap_ind] * imared) / np.sum(imared) + xcen))
                    centroid2 = float(np.round(np.sum(yy[ap_ind] * imared) / np.sum(imared) + ycen))
                    pdist = (((centroid1 - xpeak) - mean1) ** 2 + ((centroid2 - ypeak) - mean2) ** 2)
                    # pdist = np.sqrt((centroid1 - xcen)**2 + (centroid2 - ycen)**2)
                    pdist_app.append(pdist)

            tot = np.sum(pdist_app)
            Centshift = np.sqrt(tot / (num - 1)) / (r500_pix)  # np.std(pdist_app) / r500_pix

            return Centshift


        def Concentration(img, xcen, ycen, aper, frac):

            nx = np.size(img[0,])
            ny = np.size(img[:, 0])
            rep = np.repeat(np.float64(1.0e0), ny)
            rep1 = rep.reshape(ny, 1)
            rep2 = np.arange(ny)
            rep3 = rep2.reshape(ny, 1)
            xx = np.multiply(np.arange(nx), rep1) - xcen
            yy = np.multiply(np.repeat(np.float64(1.0e0), nx), rep3) - ycen
            Rad = np.sqrt(xx ** 2 + yy ** 2)
            ap = np.float64(aper)
            ap_ind = np.where(Rad <= ap)
            Flx_large = np.sum(img[ap_ind])
            ap_ind2 = np.where(Rad <= ap * frac)
            Flx_small = np.sum(img[ap_ind2])
            ap_pnd = np.where(Rad == ap)
            imared = img[ap_pnd]

            # ap_out = np.where(img < np.mean(imared))
            # img[ap_out] = 0

            return Flx_small / Flx_large


        def Gini(img, xcen, ycen, aper):

            nx = np.size(img[0,])
            ny = np.size(img[:, 0])
            rep = np.repeat(np.float64(1.0e0), ny)
            rep1 = rep.reshape(ny, 1)
            rep2 = np.arange(ny)
            rep3 = rep2.reshape(ny, 1)
            xx = np.multiply(np.arange(nx), rep1) - xcen
            yy = np.multiply(np.repeat(np.float64(1.0e0), nx), rep3) - ycen
            Rad = np.sqrt(xx ** 2 + yy ** 2)
            ap = np.float64(aper)
            ap_ind = np.where(Rad <= ap)

            n = len(img[ap_ind])
            Ksig2 = []
            Fsig2 = []
            Ki2 = []
            Fi2 = []
            Km = []
            Fm = []
            img2 = img[ap_ind]
            Rad2 = Rad[ap_ind] ** 2
            img_new = img2[np.argsort(img2)]
            Rad_new = Rad2[np.argsort(img2)]
            img_rev = img2[np.argsort(-img2)]
            Rad_rev = Rad2[np.argsort(-img2)]

            i = 1

            while i < n + 1:
                Ki = np.abs(img_new[i - 1])  # K_i
                Fi = np.abs(img_new[i - 1]) * Rad_new[i - 1]  # F_i
                Ksig = (2 * i - n - 1) * Ki
                Fsig = (2 * i - n - 1) * Fi

                Ki2.append(Ki)
                Fi2.append(Fi)
                Ksig2.append(Ksig)
                Fsig2.append(Fsig)

                i = i + 1

            Kdot = np.mean(Ki2)
            Fdot = np.mean(Fi2)
            Kisum = np.sum(Ki2)
            Fisum = np.sum(Fi2)

            j = 1

            while j < n + 1:

                K_m = np.abs(img_rev[j - 1])
                F_m = np.abs(img_rev[j - 1]) * Rad_rev[j - 1]
                Km.append(K_m)
                Fm.append(F_m)

                Kmtot = np.sum(Km)
                Fmtot = np.sum(Fm)

                if Kmtot <= 0.2 * Kisum:

                    j = j + 1

                else:
                    Km.pop()
                    Fm.pop()
                    break

            sigma_F = np.sum(Fm)
            Ksigmai = np.sum(Ksig2)
            Fsigmai = np.sum(Fsig2)
            Gin = (1 / (Kdot * n * (n - 1))) * Ksigmai
            Gin_m = (1 / (Fdot * n * (n - 1))) * Fsigmai
            M20 = np.log10(sigma_F / Fisum)

            return Gin, Gin_m, M20


        def Asymmetry(img, xcen, ycen, aper):

            img = img.astype(np.float64)
            nx = np.size(img[0,])
            ny = np.size(img[:, 0])
            rep = np.repeat(np.float64(1.0e0), ny)
            rep1 = rep.reshape(ny, 1)
            rep2 = np.arange(ny)
            rep3 = rep2.reshape(ny, 1)
            xx = np.multiply(np.arange(nx), rep1) - xcen
            yy = np.multiply(np.repeat(np.float64(1.0e0), nx), rep3) - ycen
            Rad = np.sqrt(xx ** 2 + yy ** 2)
            ap = np.float64(aper)
            ap_ind = np.where(Rad <= ap)
            center = (xcen, ycen)
            angle90 = 90
            angle180 = 180
            scale = 1.0
            M90 = cv2.getRotationMatrix2D(center, angle90, scale)
            M = cv2.getRotationMatrix2D(center, angle180, scale)
            rotated90 = cv2.warpAffine(img, M90, (nx, ny))
            rotated180 = cv2.warpAffine(img, M, (nx, ny))

            S = np.sum(img[ap_ind])
            T = np.sum(np.abs(img[ap_ind] - rotated180[ap_ind]))
            T90 = np.sum(np.abs(img[ap_ind] - rotated90[ap_ind]))

            # plt.imshow(rotated180, cmap = 'viridis', origin = 'lower', interpolation = 'none', norm = LogNorm())  #viridis or spectral
            # plt.plot(xcen, ycen, '.', color = 'red')
            # plt.colorbar(orientation = 'vertical')
            # plt.xlim(130, 815)
            # plt.ylim(160, 855)
            return T / S, T90 / S


        def Power(img, xcen, ycen, aper, minorder, maxorder):

            nx = np.size(img[0,])
            ny = np.size(img[:, 0])
            rep = np.repeat(np.float64(1.0e0), ny)
            rep1 = rep.reshape(ny, 1)
            rep2 = np.arange(ny)
            rep3 = rep2.reshape(ny, 1)
            xx = np.multiply(np.arange(nx), rep1) - xcen
            yy = np.multiply(np.repeat(np.float64(1.0e0), nx), rep3) - ycen
            Rad = np.sqrt(xx ** 2 + yy ** 2)
            ap = np.float64(aper)
            i = (xx != 0)
            yy[i] /= xx[i]
            phi = np.arctan(yy)
            j = (xx < 0)
            j_n = len(j)

            if j_n > 0:
                phi[j] += np.pi
            k = (xx > 0) & (yy < 0)
            k_n = len(k)

            if k_n > 0:
                phi[k] += 2 * np.pi
            l = (xx == 0) & (yy >= 0)
            l_n = len(l)

            if l_n > 0:
                phi[l] = np.pi / 2
            q = (xx == 0) & (yy < 0)
            q_n = len(q)

            if q_n > 0:
                phi[q] = (3 * np.pi) / 2

            ap_ind = np.where(Rad <= ap)
            imared = img[ap_ind]
            radred = Rad[ap_ind]
            phired = phi[ap_ind]
            num = np.arange(minorder, maxorder + 1)
            po = []
            a_0 = np.sum(imared)
            P_0 = (a_0 * np.log(ap)) ** 2

            for m in num:
                a_m = np.sum(imared * radred ** m * np.cos(m * phired))
                b_m = np.sum(imared * radred ** m * np.sin(m * phired))
                P_m = (a_m ** 2 + b_m ** 2) / (2 * (m ** 2) * ap ** (2 * m))
                po.append(P_m)

            power = po / P_0

            return power

        Cent_1 = Center(sim_2_1, xcen_ini, ycen_ini, r500_pix)
        Cent_2 = Center(sim_2_1, xcen_ini, ycen_ini, r500_pix)
        xcen_mass = Cent_1[0]
        ycen_mass = Cent_1[1]
        xcen_peak = Cent_1[2]
        ycen_peak = Cent_1[3]
        xcen_mass_2 = Cent_2[0]
        ycen_mass_2 = Cent_2[1]
        xcen_peak_2 = Cent_2[2]
        ycen_peak_2 = Cent_2[3]
        Cen_1 = Centroid(sim_2_1, xcen_mass, ycen_mass, xcen_peak, ycen_peak, r500_pix, 10)
        Con_1 = Concentration(sim_2_1, xcen_mass, ycen_mass, r500_pix, 0.1)
        Gin_1 = Gini(sim_2_1, xcen_mass, ycen_mass, r500_pix)
        Asy_1 = Asymmetry(sim_2_1, xcen_mass, ycen_mass, r500_pix)
        Pow_1 = Power(sim_2_1, xcen_mass, ycen_mass, r500_pix, 1, 4)

        Cen1.append(Cen_1)
        Con1.append(Con_1)
        Gin1.append(Gin_1[0])
        M201.append(Gin_1[2])
        Asy180.append(Asy_1[0])
        Asy90.append(Asy_1[1])
        Pow2.append(Pow_1[1])
        Pow3.append(Pow_1[2])
        Pow4.append(Pow_1[3])

        k = k + 1

    Cen_me = stats.trim_mean(Cen_raw_1 * (Cen_sim_1 / Cen1), 0.25)
    Cen_up = np.percentile(Cen_raw_1 * (Cen_sim_1 / Cen1), 50 + 68.3 / 2)  # << upper 1sigma value
    Cen_lo = np.percentile(Cen_raw_1 * (Cen_sim_1 / Cen1), 50 - 68.3 / 2)  # << lower 1sigma value

    Con_me = stats.trim_mean(Con_raw_1 * (Con_sim_1 / Con1), 0.25)
    Con_up = np.percentile(Con_raw_1 * (Con_sim_1 / Con1), 50 + 68.3 / 2)  # << upper 1sigma value
    Con_lo = np.percentile(Con_raw_1 * (Con_sim_1 / Con1), 50 - 68.3 / 2)  # << lower 1sigma value

    Gin_me = stats.trim_mean(Gin_raw_1[0] * (Gin_sim_1[0] / Gin1), 0.25)
    Gin_up = np.percentile(Gin_raw_1[0] * (Gin_sim_1[0] / Gin1), 50 + 68.3 / 2)  # << upper 1sigma value
    Gin_lo = np.percentile(Gin_raw_1[0] * (Gin_sim_1[0] / Gin1), 50 - 68.3 / 2)  # << lower 1sigma value

    M20_me = stats.trim_mean(Gin_raw_1[2] * (Gin_sim_1[2] / M201), 0.25)
    M20_up = np.percentile(Gin_raw_1[2] * (Gin_sim_1[2] / M201), 50 + 68.3 / 2)  # << upper 1sigma value
    M20_lo = np.percentile(Gin_raw_1[2] * (Gin_sim_1[2] / M201), 50 - 68.3 / 2)  # << lower 1sigma value

    Asy180_me = stats.trim_mean(Asy_raw_1[0] * (Asy_sim_1[0] / Asy180), 0.25)
    Asy180_up = np.percentile(Asy_raw_1[0] * (Asy_sim_1[0] / Asy180), 50 + 68.3 / 2)  # << upper 1sigma value
    Asy180_lo = np.percentile(Asy_raw_1[0] * (Asy_sim_1[0] / Asy180), 50 - 68.3 / 2)  # << lower 1sigma value

    Asy90_me = stats.trim_mean(Asy_raw_1[1] * (Asy_sim_1[1] / Asy90), 0.25)
    Asy90_up = np.percentile(Asy_raw_1[1] * (Asy_sim_1[1] / Asy90), 50 + 68.3 / 2)  # << upper 1sigma value
    Asy90_lo = np.percentile(Asy_raw_1[1] * (Asy_sim_1[1] / Asy90), 50 - 68.3 / 2)  # << lower 1sigma value

    Pow2_me = stats.trim_mean(Pow_raw_1[1] * (Pow_sim_1[1] / Pow2), 0.25)
    Pow2_up = np.percentile(Pow_raw_1[1] * (Pow_sim_1[1] / Pow2), 50 + 68.3 / 2)  # << upper 1sigma value
    Pow2_lo = np.percentile(Pow_raw_1[1] * (Pow_sim_1[1] / Pow2), 50 - 68.3 / 2)  # << lower 1sigma value

    Pow3_me = stats.trim_mean(Pow_raw_1[2] * (Pow_sim_1[2] / Pow3), 0.25)
    Pow3_up = np.percentile(Pow_raw_1[2] * (Pow_sim_1[2] / Pow3), 50 + 68.3 / 2)  # << upper 1sigma value
    Pow3_lo = np.percentile(Pow_raw_1[2] * (Pow_sim_1[2] / Pow3), 50 - 68.3 / 2)  # << lower 1sigma value

    Pow4_me = stats.trim_mean(Pow_raw_1[3] * (Pow_sim_1[3] / Pow4), 0.25)
    Pow4_up = np.percentile(Pow_raw_1[3] * (Pow_sim_1[3] / Pow4), 50 + 68.3 / 2)  # << upper 1sigma value
    Pow4_lo = np.percentile(Pow_raw_1[3] * (Pow_sim_1[3] / Pow4), 50 - 68.3 / 2)  # << lower 1sigma value

    fig1 = plt.figure(1)
    sub1 = fig1.add_subplot(3, 1, 1)
    plt.title('Centroid Shift')

    y1, x1, p1 = plt.hist(Cen_raw_1 * (Cen_sim_1 / Cen1), density = True, histtype='stepfilled', alpha = 0.3, bins = 20)
    ax1 = sns.distplot(Cen_raw_1 * (Cen_sim_1 / Cen1), norm_hist=True)
    x10 = ax1.lines[0].get_xdata()
    y10 = ax1.lines[0].get_ydata()
    maxid1 = np.argmax(y10)
    loc1 = x10[maxid1]

    plt.axvline(x = Cen_me, color='red', linestyle='-', linewidth=1.2, label=r'T = {:.3g}'.format(Cen_me))
    plt.axvline(x = Cen_up, color='red', linestyle='--', linewidth=1.2, label=r'+ {:.3g}'.format(float(Cen_up - Cen_me)))
    plt.axvline(x = Cen_lo, color='red', linestyle='--', linewidth=1.2, label=r'$-$ {:.3g}'.format(float(Cen_me - Cen_lo)))

    plt.legend(numpoints=1, frameon=True, loc='upper right', prop={'size': 7.5})
    plt.yticks([])


    sub2 = fig1.add_subplot(3, 1, 2)
    plt.title('Concentration')

    y2, x2, p2 = plt.hist(Con_raw_1 * (Con_sim_1 / Con1), density = True, histtype='stepfilled', alpha = 0.3, bins = 20)
    ax2 = sns.distplot(Con_raw_1 * (Con_sim_1 / Con1), norm_hist=True)
    x20 = ax2.lines[0].get_xdata()
    y20 = ax2.lines[0].get_ydata()
    maxid2 = np.argmax(y20)
    loc2 = x20[maxid2]

    plt.axvline(x = Con_me, color='red', linestyle='-', linewidth=1.2, label=r'T = {:.3g}'.format(Con_me))
    plt.axvline(x = Con_up, color='red', linestyle='--', linewidth=1.2, label=r'+ {:.3g}'.format(float(Con_up - Con_me)))
    plt.axvline(x = Con_lo, color='red', linestyle='--', linewidth=1.2, label=r'$-$ {:.3g}'.format(float(Con_me - Con_lo)))

    plt.legend(numpoints=1, frameon=True, loc='upper right', prop={'size': 7.5})
    plt.yticks([])
    #plt.xlim([Con_me - ((Con_me - Con_lo) * 3), Con_me - ((Con_up - Con_me) * 3)])

    sub3 = fig1.add_subplot(3, 1, 3)
    plt.title('Gini Coefficient')

    y3, x3, p3 = plt.hist(Gin_raw_1[0] * (Gin_sim_1[0] / Gin1), density = True, histtype='stepfilled', alpha = 0.3, bins = 20)
    ax3 = sns.distplot(Gin_raw_1[0] * (Gin_sim_1[0] / Gin1), norm_hist=True)
    x30 = ax3.lines[0].get_xdata()
    y30 = ax3.lines[0].get_ydata()
    maxid3 = np.argmax(y30)
    loc3 = x30[maxid3]

    plt.axvline(x = Gin_me, color='red', linestyle='-', linewidth=1.2, label=r'T = {:.3g}'.format(Gin_me))
    plt.axvline(x = Gin_up, color='red', linestyle='--', linewidth=1.2, label=r'+ {:.3g}'.format(float(Gin_up - Gin_me)))
    plt.axvline(x = Gin_lo, color='red', linestyle='--', linewidth=1.2, label=r'$-$ {:.3g}'.format(float(Gin_me - Gin_lo)))

    plt.legend(numpoints=1, frameon=True, loc='upper right', prop={'size': 7.5})
    plt.yticks([])
    #plt.xlim([Gin_me - ((Gin_me - Gin_lo) * 3), Gin_me - ((Gin_up - Gin_me) * 3)])
    plt.subplots_adjust(hspace=0.9)
    plt.savefig("D:/Master thesis/Final_refill/Fig5(0.2)/{}.png".format(Name + '_1'), dpi = 300)
    plt.close()

    fig2 = plt.figure(2)
    sub4 = fig2.add_subplot(3, 1, 1)
    plt.title('M20')

    y4, x4, p4 = plt.hist(Gin_raw_1[2] * (Gin_sim_1[2] / M201), density = True, histtype='stepfilled', alpha = 0.3, bins = 20)
    ax4 = sns.distplot(Gin_raw_1[2] * (Gin_sim_1[2] / M201), norm_hist=True)
    x40 = ax4.lines[0].get_xdata()
    y40 = ax4.lines[0].get_ydata()
    maxid4 = np.argmax(y40)
    loc4 = x40[maxid4]

    plt.axvline(x = M20_me, color='red', linestyle='-', linewidth=1.2, label=r'T = {:.4g}'.format(M20_me))
    plt.axvline(x = M20_up, color='red', linestyle='--', linewidth=1.2, label=r'+ {:.4g}'.format(float(M20_up - M20_me)))
    plt.axvline(x = M20_lo, color='red', linestyle='--', linewidth=1.2, label=r'$-$ {:.4g}'.format(float(M20_me - M20_lo)))

    plt.legend(numpoints=1, frameon=True, loc='upper right', prop={'size': 7.5})
    plt.yticks([])
    #plt.xlim([M20_me - ((M20_me - M20_lo) * 3), M20_me - ((M20_up - M20_me) * 3)])

    sub5 = fig2.add_subplot(3, 1, 2)
    plt.title('Asymmetry Parameter 180')

    y5, x5, p5 = plt.hist(Asy_raw_1[0] * (Asy_sim_1[0] / Asy180), density = True, histtype='stepfilled', alpha = 0.3, bins = 20)
    ax5 = sns.distplot(Asy_raw_1[0] * (Asy_sim_1[0] / Asy180), norm_hist=True)
    x50 = ax5.lines[0].get_xdata()
    y50 = ax5.lines[0].get_ydata()
    maxid5 = np.argmax(y50)
    loc5 = x50[maxid5]

    plt.axvline(x = Asy180_me, color='red', linestyle='-', linewidth=1.2, label=r'T = {:.3g}'.format(Asy180_me))
    plt.axvline(x = Asy180_up, color='red', linestyle='--', linewidth=1.2, label=r'+ {:.3g}'.format(float(Asy180_up - Asy180_me)))
    plt.axvline(x = Asy180_lo, color='red', linestyle='--', linewidth=1.2, label=r'$-$ {:.3g}'.format(float(Asy180_me - Asy180_lo)))

    plt.legend(numpoints=1, frameon=True, loc='upper right', prop={'size': 7.5})
    plt.yticks([])
    #plt.xlim([Asy180_me - ((Asy180_me - Asy180_lo) * 3), Asy180_me - ((Asy180_up - Asy180_me) * 3)])

    sub6 = fig2.add_subplot(3, 1, 3)
    plt.title('Asymmetry Parameter 90')

    y6, x6, p6 = plt.hist(Asy_raw_1[1] * (Asy_sim_1[1] / Asy90), density = True, histtype='stepfilled', alpha = 0.3, bins = 20)
    ax6 = sns.distplot(Asy_raw_1[1] * (Asy_sim_1[1] / Asy90), norm_hist=True)
    x60 = ax6.lines[0].get_xdata()
    y60 = ax6.lines[0].get_ydata()
    maxid6 = np.argmax(y60)
    loc6 = x60[maxid6]

    plt.axvline(x = Asy90_me, color='red', linestyle='-', linewidth=1.2, label=r'T = {:.3g}'.format(Asy90_me))
    plt.axvline(x = Asy90_up, color='red', linestyle='--', linewidth=1.2, label=r'+ {:.3g}'.format(float(Asy90_up - Asy90_me)))
    plt.axvline(x = Asy90_lo, color='red', linestyle='--', linewidth=1.2, label=r'$-$ {:.3g}'.format(float(Asy90_me - Asy90_lo)))

    plt.legend(numpoints=1, frameon=True, loc='upper right', prop={'size': 7.5})
    plt.yticks([])
    #plt.xlim([Asy90_me - ((Asy90_me - Asy90_lo) * 3), Asy90_me - ((Asy90_up - Asy90_me) * 3)])

    plt.subplots_adjust(hspace=0.9)
    plt.savefig("D:/Master thesis/Final_refill/Fig5(0.2)/{}.png".format(Name + '_2'), dpi = 300)
    plt.close()

    fig3 = plt.figure(3)

    sub7 = fig3.add_subplot(3, 1, 1)
    plt.title('Power Ratio P2/P0')

    y7, x7, p7 = plt.hist(Pow_raw_1[1] * (Pow_sim_1[1] / Pow2), density = True, histtype='stepfilled', alpha = 0.3, bins = 20)
    ax7 = sns.distplot(Pow_raw_1[1] * (Pow_sim_1[1] / Pow2), norm_hist=True)
    x70 = ax7.lines[0].get_xdata()
    y70 = ax7.lines[0].get_ydata()
    maxid7 = np.argmax(y70)
    loc7 = x70[maxid7]

    plt.axvline(x = Pow2_me, color='red', linestyle='-', linewidth=1.2, label=r'T = {:.3g}'.format(Pow2_me))
    plt.axvline(x = Pow2_up, color='red', linestyle='--', linewidth=1.2, label=r'+ {:.3g}'.format(float(Pow2_up - Pow2_me)))
    plt.axvline(x = Pow2_lo, color='red', linestyle='--', linewidth=1.2, label=r'$-$ {:.3g}'.format(float(Pow2_me - Pow2_lo)))

    plt.legend(numpoints=1, frameon=True, loc='upper right', prop={'size': 7.5})
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    plt.yticks([])


    sub8 = fig3.add_subplot(3, 1, 2)
    plt.title('Power Ratio P3/P0')

    y8, x8, p8 = plt.hist(Pow_raw_1[2] * (Pow_sim_1[2] / Pow3), density = True, histtype='stepfilled', alpha = 0.3, bins = 20)
    ax8 = sns.distplot(Pow_raw_1[2] * (Pow_sim_1[2] / Pow3), norm_hist=True)
    x80 = ax8.lines[0].get_xdata()
    y80 = ax8.lines[0].get_ydata()
    maxid8 = np.argmax(y80)
    loc8 = x80[maxid8]

    plt.axvline(x = Pow3_me, color='red', linestyle='-', linewidth=1.2, label=r'T = {:.3g}'.format(Pow3_me))
    plt.axvline(x = Pow3_up, color='red', linestyle='--', linewidth=1.2, label=r'+ {:.3g}'.format(float(Pow3_up - Pow3_me)))
    plt.axvline(x = Pow3_lo, color='red', linestyle='--', linewidth=1.2, label=r'$-$ {:.3g}'.format(float(Pow3_me - Pow3_lo)))

    plt.legend(numpoints=1, frameon=True, loc='upper right', prop={'size': 7.5})
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    plt.yticks([])


    sub9 = fig3.add_subplot(3, 1, 3)
    plt.title('Power Ratio P4/P0')

    y9, x9, p9 = plt.hist(Pow_raw_1[3] * (Pow_sim_1[3] / Pow4), density = True, histtype='stepfilled', alpha = 0.3, bins = 20)
    ax9 = sns.distplot(Pow_raw_1[3] * (Pow_sim_1[3] / Pow4), norm_hist=True)
    x90 = ax9.lines[0].get_xdata()
    y90 = ax9.lines[0].get_ydata()
    maxid9 = np.argmax(y90)
    loc9 = x90[maxid9]

    plt.axvline(x = Pow4_me, color='red', linestyle='-', linewidth=1.2, label=r'T = {:.3g}'.format(Pow4_me))
    plt.axvline(x = Pow4_up, color='red', linestyle='--', linewidth=1.2, label=r'+ {:.3g}'.format(float(Pow4_up - Pow4_me)))
    plt.axvline(x = Pow4_lo, color='red', linestyle='--', linewidth=1.2, label=r'$-$ {:.3g}'.format(float(Pow4_me - Pow4_lo)))

    plt.legend(numpoints=1, frameon=True, loc='upper right', prop={'size': 7.5})
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    plt.yticks([])

    plt.subplots_adjust(hspace=0.9)

    plt.savefig("D:/Master thesis/Final_refill/Fig5(0.2)/{}.png".format(Name + '_3'), dpi = 300)
    plt.close()

    write_ws.append(
            ["{}".format(Name),
             "{:.3g}".format(Cen_me), "{:.3g}".format(float(Cen_up - Cen_me)), "{:.3g}".format(float(Cen_me - Cen_lo)),
             "{:.3g}".format(Con_me), "{:.3g}".format(float(Con_up - Con_me)), "{:.3g}".format(float(Con_me - Con_lo)),
             "{:.3g}".format(Gin_me), "{:.3g}".format(float(Gin_up - Gin_me)), "{:.3g}".format(float(Gin_me - Gin_lo)),
             "{:.4g}".format(M20_me), "{:.4g}".format(float(M20_up - M20_me)), "{:.4g}".format(float(M20_me - M20_lo)),
             "{:.3g}".format(Asy180_me), "{:.3g}".format(float(Asy180_up - Asy180_me)), "{:.3g}".format(float(Asy180_me - Asy180_lo)),
             "{:.3g}".format(Asy90_me), "{:.3g}".format(float(Asy90_up - Asy90_me)), "{:.3g}".format(float(Asy90_me - Asy90_lo)),
             "{:.3g}".format(Pow2_me), "{:.3g}".format(float(Pow2_up - Pow2_me)), "{:.3g}".format(float(Pow2_me - Pow2_lo)),
             "{:.3g}".format(Pow3_me), "{:.3g}".format(float(Pow3_up - Pow3_me)), "{:.3g}".format(float(Pow3_me - Pow3_lo)),
             "{:.3g}".format(Pow4_me), "{:.3g}".format(float(Pow4_up - Pow4_me)), "{:.3g}".format(float(Pow4_me - Pow4_lo))]),
    write_wb.save("D:/Master thesis/Final_refill/(0.2).xlsx")

