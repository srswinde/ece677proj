from astropy.io import fits
from fits_solver import m4k_imclient
import matplotlib.pyplot as plt
import scipy.spatial
import numpy as np
import math
import struct
from pyds9 import DS9;D=DS9()


class quad_set:

    def __init__( self, fname="pointing0064_merged.fits", thresh=5.0 ):
        self.xymat = self.sextract(fname, thresh)
        self.quads = self.make_quads(self.xymat)
        self.rquads = self.allhash(self.quads, self.xymat)


    def sextract(self, fname="pointing0064_merged.fits", thresh=5.0):
        self.fitsfd = fits.open(fname)
        objs = m4k_imclient.getobjects(self.fitsfd[0].data, thresh)
        self.df = m4k_imclient.mkdataframe(objs)

        xymat = self.df[["xpeak", "ypeak"]].as_matrix()

        return xymat


    def make_quads(self, xymat):

        quads = []
        cpxymat = xymat.copy()
        while len(cpxymat) > 4:

            cpquad = scipy.spatial.KDTree(cpxymat).query( cpxymat[0], 4 )[1]
            quad = []
            #get xymat index
            for ii in [0, 1, 2, 3]:
                idx = int(np.where(
                    (xymat == cpxymat[cpquad[ii]]).all(axis=1))[0])
                quad.append( idx )
            quads.append( quad )
            cpxymat = np.delete(cpxymat, cpquad, 0)

        return quads


    def plot_quads(self, quads, xymat):

        color = "rgbcmykw"
        for ii, q in enumerate(quads):
            sub = xymat[q]
            x = sub[:, 0]
            y = sub[:, 1]
            plt.plot(x, y, color[ii % len(color)]+'o-')

    def __getitem__(self, key):
        return self.quadpts(key)

    def quadhash(self, quad_points, return_hash=False):

        dmat = scipy.spatial.distance_matrix(quad_points, quad_points)
        maxd = np.unravel_index(np.argmax(dmat), dmat.shape)

        #origin at A
        quad_points = (quad_points - quad_points[maxd[0]]).astype(float)

        B = quad_points[maxd[1]]

        angle = math.pi/4.0 - math.atan2(B[1], B[0])
        rotmat = np.array(
            [[math.cos(angle), -math.sin(angle)],
             [math.sin(angle), math.cos(angle)]])
        qps = np.dot(rotmat, quad_points.T).T
        qps = qps/qps[maxd[1], 0]

        if return_hash:
            return struct.pack("4f", *qps[1:3].ravel())
        else:
            return qps[1:3].ravel()


    def allhash( self, quads, xymat, return_hash=False ):
        if return_hash:
            bts = b''
            for q in quads:
                bts += self.quadhash(xymat[q], True)
            return bts
        else:
            pts = []
            for q in quads:
                pts.append(self.quadhash(xymat[q]))
            return np.array(pts)

    def binhash(self):
        bts = b''
        for q in self.rquads:
            bts += struct.pack("4f", *q)
        return bts

    def quadpts(self, key):
        return self.xymat[self.quads[key]]

    def display(self, ds=D):
        ds.set_pyfits(self.fitsfd)

    def draw_quad(self, idx, ds=D):
        qps = self.__getitem__(idx)
        for pts in qps:
            D.set("regions command {{circle {} {} 5}}".format(*pts))

    def draw_all(self):
        for pt in self.xymat:
            x, y = pt
            D.set("regions command {{circle {} {} 5}}".format(x, y))


def setup():
    return quad_set(), quad_set("skv625064874090.fits", thresh=20.0)


def compare_quads( rquads1, rquads2 ):
    isclose = []
    import time
    t0 = time.time()
    for i1, q1 in enumerate(rquads1.rquads):
        for i2, q2 in enumerate(rquads2.rquads):
            if np.allclose(q1, q2, atol=0.01, rtol=0) and \
                    not np.allclose(q1, np.array([0., 0., 1., 1.])):

                isclose.append((i1, i2))
        if time.time()-t0 > 4*60:
            break

    return isclose


def check(qs1, qs2, matches):
    m1, m2 = matches[0]
    pts1 = qs1[m1]
    pts2 = qs2[m2]

    minx1 = min(pts1[0, :])-50
    minx2 = min(pts2[0, :])-50

    maxx1 = max(pts1[0, :])+50
    maxx2 = max(pts2[0, :])+50

    miny1 = min(pts1[1, :])-50
    miny2 = min(pts2[1, :])-50

    maxy1 = max(pts1[1, :])+50
    maxy2 = max(pts2[1, :])+50

    fig1=plt.figure()
    plt.imshow(qs1.fitsfd[0].data[minx1:maxx1, miny1:maxy1])
    fig2=plt.figure()
    plt.imshow(qs1.fitsfd[0].data[minx2:maxx2, miny2:maxy2])
