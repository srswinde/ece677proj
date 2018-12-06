from astropy.io import fits
from fits_solver import m4k_imclient
import matplotlib.pyplot as plt
import scipy.spatial
import numpy as np
import math
import struct
from pyds9 import DS9;D=DS9()



class quad_set:
    """
        Class to analyze image and build quads from stellar sources.
    """
    def __init__( self, fname="pointing0064_merged.fits", thresh=5.0 ):
        """
        Load the image and extract stellar sources them using SEP
        any my own repo at:
        https://github.com/srswinde/fits_solver
        """

        #numpy array of pixel coordinates of sources
        self.xymat = self.sextract(fname, thresh)

        self.quads = self.make_quads2(self.xymat)

        # create an array of quads in the geometrically hashed form.
        self.rquads = self.allhash(self.quads, self.xymat)
        self.fname = fname


    def sextract(self, fname="pointing0064_merged.fits", thresh=5.0):
        """
            Do the Source Extraction
        """
        self.fitsfd = fits.open(fname)
        objs = m4k_imclient.getobjects(self.fitsfd[0].data, thresh)
        self.df = m4k_imclient.mkdataframe(objs)

        xymat = self.df[["xpeak", "ypeak"]].as_matrix()

        return xymat


    def make_quads(self, xymat):
        """
        Old make quads. Uses each star only once.
        """
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
            cpxymat = np.delete( cpxymat, cpquad, 0 )

        return quads

    def make_quads2(self, xymat=None):
        """
        creates a 10 x 10 grid of sections of the image
        All quads are taken from within these cells of data.
        """

        if xymat is None: xymat = self.xymat
        quads = []
        cpxymat = xymat.copy()
        ptree = scipy.spatial.KDTree(cpxymat)
        width, height = self.fitsfd[0].data.shape
        radius = math.sqrt((width//10)**2 + (height//10)**2)
        #10x10 grid to make quads in
        xs = scipy.linspace(0, width, 10, dtype=int)
        ys = scipy.linspace(0, height, 10, dtype=int)
        gridpoints = []

        for x in xs.ravel():
            for y in ys.ravel():
                pts = ptree.query_ball_point( (x, y), r=radius)
                for pt1 in pts:
                    for pt2 in pts[:3]:
                        if pt1 == pt2: continue
                        for pt3 in pts[:3]:
                            if len(set((pt1, pt2, pt3))) != 3: continue
                            for pt4 in pts[:3]:
                                if len(set((pt1, pt2, pt3, pt4))) != 4: continue
                                quads.append(np.array([pt1, pt2, pt3, pt4], dtype=int))

        return np.array(quads)


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

    def save(self):
        basename = self.fname.split('.')[0]
        with open(basename+".bin", "wb") as binfd:
            binfd.write(self.binhash())

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
    print(i1, i2)
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


def runtest():
    fname1 = "pointing0064_merged.fits"
    fname2 = "skv625064874090.fits"
    q1 = quad_set(fname1)
    q2 = quad_set(fname2)
    q1.save()
    q2.save()

    #print( compare_quads(q1, q2) )
    return q1, q2
