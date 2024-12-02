import numpy as np
import scipy as sp
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

#class Arrow3D(FancyArrowPatch):
#    def __init__(self, xs, ys, zs, *args, **kwargs):
#        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
#        self._verts3d = xs, ys, zs
#
#    def draw(self, renderer):
#        xs3d, ys3d, zs3d = self._verts3d
#        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
#        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
#        FancyArrowPatch.draw(self, renderer)
class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        super().__init__((0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def do_3d_projection(self, renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))

        return np.min(zs)

def spotBounds(geo,atnums,ecs):
    bounds_x = list()
    bounds_y = list()
    bounds_z = list()
    for i,g1 in enumerate(geo):
        for j,g2 in enumerate(geo):
            if j >= i:
                continue
            sumCovRad = (ecs[atnums[i]][3] + ecs[atnums[j]][3]) * 0.01
            dist = ((g1[0]-g2[0])**2 +
                    (g1[1]-g2[1])**2 +
                    (g1[2]-g2[2])**2)**0.5
            if dist < sumCovRad + 0.1*sumCovRad:
                bounds_x.append((g1[0],g2[0]))
                bounds_y.append((g1[1],g2[1]))
                bounds_z.append((g1[2],g2[2]))
    return bounds_x, bounds_y, bounds_z

def plot_cdft(ofile,cdfts,ecs):
    cdft_fnames = [(r'f$^{+}$','f+'),(r'f$^{-}$','f-'),('\u0394f','delta_f'),
                   (r's$^{+}$','s+'),(r's$^{-}$','s-'),('\u0394s','delta_s'),
                   (u'\u0394\u03C1$_{elec}$','delta_rho_elec'),
                   (u'\u0394\u03C1$_{nuc}$','delta_rho_nuc')]
    scale = 5
    chtypes = list()
    for cht in cdfts[1]:
        chtypes.append(cht)
    print('Enter a number to chose a function:')
    for c,cfn in enumerate(cdft_fnames):
        print(f'-{c+1}: {cfn[1]}')
    n = input() 
    while True:
        try:
            n = int(n)-1
            if n < len(cdft_fnames) and n >= 0:
                break
            else:
                print('Please, enter an integer between 1 and '
                     f'{len(cdft_fnames)}')
                n = input() 
        except:
            print('Please, enter an integer')
            n = input() 
    print('Enter a number to chose the type of charges:')
    for c,cht in enumerate(chtypes):
        print(f'-{c+1}: {cht}')
    cn = input()
    while True:
        try:
            cn = int(cn) - 1
            if cn < len(chtypes) and cn >= 0:
                break
            else:
                print(f'Please, enter an integer between 1 and {len(chtypes)}')
                cn = input() 
        except:
            print('Please, enter an integer')
            cn = input() 
    cht = chtypes[cn]

    xs = list()
    ys = list()
    zs = list()
    size = list()
    labels = list()
    for e,elt in enumerate(ofile.coordinates):
        xs.append(elt[0])
        ys.append(elt[1])
        zs.append(elt[2])
        size.append(float(ecs[ofile.atnums[e]][3])*scale)
        value = round(cdfts[n][cht][e],3)
        labels.append(f"({ofile.atnums[e]}) {e+1}: {value}")

    bxs,bys,bzs = spotBounds(ofile.coordinates,ofile.atnums,ecs)

    cmap = plt.cm.rainbow
    fig = plt.figure()
    title = f"{cdft_fnames[n][0]} with {cht} charges"
    fig.suptitle(title)
    ax = fig.add_subplot(111, projection='3d', facecolor="black")
    # Set axis ratio
    ax.set_box_aspect((np.ptp(xs), np.ptp(ys), np.ptp(zs)))
    scatter = ax.scatter(xs,ys,zs,s=size,c=cdfts[n][cht],cmap=cmap)
    for b,bnd in enumerate(bxs):
        ax.plot(bnd, bys[b], bzs[b], color='black', linewidth=2)
    for l,lab in enumerate(labels):
        ax.text(xs[l],ys[l],zs[l], "{}".format(lab),color='white')
    cbar = plt.colorbar(scatter)
    cbar.ax.set_ylabel(f'{cdft_fnames[n][0]}', labelpad=10, rotation=270)
    plt.show()
    plt.close()
    
def superposeGeos(geo1,geo2,rmsd,ecs,label=True):
    #bounds_x, bounds_y, bounds_z = spotBounds(xs[mol], ys[mol], zs[mol])
    bxs1,bys1,bzs1 = spotBounds(geo1,ecs)
    geos = dict()
    geos['geo1'] = dict()
    geos['geo2'] = dict()
    for l,line1 in enumerate(geo1):
        geos['geo1'][l+1] = [line1.split()[0], 
                             [float(c) for c in line1.split()[1:4]]]
        line2 = geo2[l]
        geos['geo2'][l+1] = [line2.split()[0], 
                             [float(c) for c in line2.split()[1:4]]]
    size = 3
    xs = dict()
    ys = dict()
    zs = dict()
    labels = dict()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d', facecolor="orange")
    clr = ['red','green']
    mol1 = list(geos.keys())[0]
    mol2 = list(geos.keys())[1]
    atomsK = list()
    for at in geos[mol1]:
        atomsK.append(at)
    geoK1 = np.zeros((len(atomsK),3))
    geoK2 = np.zeros((len(atomsK),3))
    for i,at in enumerate(atomsK):
        for j,c in enumerate(geos[mol1][at][1]):
            geoK1[i][j] = c
            geoK2[i][j] = geos[mol2][at][1][j]
        
    rot,rssd = sp.spatial.transform.Rotation.align_vectors(geoK1,geoK2)
    print(rssd)
    rot_mat = rot.as_matrix()

    ageos = dict()
    ageos[mol1] = dict()
    ageos[mol2] = dict()
    for at in geos[mol1]:
        ageos[mol1][at] = geos[mol1][at]
    for at in geos[mol2]:
        ageos[mol2][at] = list()
        ageos[mol2][at].append(geos[mol2][at][0])
        newC = rot_mat.dot(geos[mol2][at][1])
        ageos[mol2][at].append(newC)

    #Sorry
    geo2_bnd = list()
    for at in ageos[mol2]:
        line = f'{ageos[mol2][at][0]} '
        for coord in ageos[mol2][at][1]:
            line += f'{coord} '
        geo2_bnd.append(line)
    bxs2,bys2,bzs2 = spotBounds(geo2_bnd,ecs)
    
    for m,mol in enumerate(ageos):
        xs[mol] = dict()
        ys[mol] = dict()
        zs[mol] = dict()
        labels[mol] = dict()
        #Classify atoms according to their nature
        for at in ageos[mol]:
            g = ageos[mol][at]
            if not g[0] in xs[mol]:
                xs[mol][g[0]] = list()
            xs[mol][g[0]].append(float(g[1][0]))
            if not g[0] in ys[mol]:
                ys[mol][g[0]] = list()
            ys[mol][g[0]].append(float(g[1][1]))
            if not g[0] in zs[mol]:
                zs[mol][g[0]] = list()
            zs[mol][g[0]].append(float(g[1][2]))
            if not g[0] in labels[mol]:
                labels[mol][g[0]] = list()
            labels[mol][g[0]].append(at)
        upLim = 0
        downLim = 0
        for c in [xs[mol],ys[mol],zs[mol]]:
            for t in c:
                cOrd = c[t].copy()
                cOrd.sort()
                if cOrd[0] < downLim:
                    downLim = cOrd[0]
                if cOrd[-1] > upLim:
                    upLim = cOrd[-1]
        ax.set_xlim3d(downLim,upLim)
        ax.set_ylim3d(downLim,upLim)
        ax.set_zlim3d(downLim,upLim)
        for t in xs[mol]:
            for elt in ecs:
                if t == ecs[elt][0]:
                    e = elt
                    break
            ax.scatter(xs[mol][t], ys[mol][t], zs[mol][t],
                       c=clr[m],s=int(ecs[e][3])*size,
                       depthshade=False,alpha=0.5)
            if label:
                for k,l in enumerate(labels[mol][t]):
                    ax.text(xs[mol][t][k],ys[mol][t][k],zs[mol][t][k],
                            "{}".format(l),color='white')
        ax.annotate(f'rsmd = {rmsd}', xy=(0.05, 0.95), 
                    xycoords='axes fraction')
        for b,bnd in enumerate(bxs1):
            ax.plot(bnd, bys1[b], bzs1[b], color='black', linewidth=2)
        for b,bnd in enumerate(bxs2):
            ax.plot(bnd, bys2[b], bzs2[b], color='black', linewidth=2)
    plt.show()
    plt.close()

def plot_nv(ofile,ecs):
    scale = 5
    scale_vec = 2
    xs = list()
    ys = list()
    zs = list()
    size = list()
    colors = list()
    nvs = list()
    labels = list()
    for e,elt in enumerate(ofile.coordinates.rstrip().split('\n')):
        nvx = (float(elt.split()[1]),
               float(elt.split()[1])+ofile.normvec[e][1][0]*scale_vec)
        nvy = (float(elt.split()[2]),
               float(elt.split()[2])+ofile.normvec[e][1][1]*scale_vec)
        nvz = (float(elt.split()[3]),
               float(elt.split()[3])+ofile.normvec[e][1][2]*scale_vec)
        nvs.append((nvx,nvy,nvz)) 
        xs.append(float(elt.split()[1]))
        ys.append(float(elt.split()[2]))
        zs.append(float(elt.split()[3]))
        for f in ecs:
            if elt.split()[0] == ecs[f][0]:
                ne = f
                break
        size.append(float(ecs[ne][3])*scale)
        colors.append(ecs[ne][2])
        labels.append(f"{e+1}")

    bxs,bys,bzs = spotBounds(ofile.coordinates.rstrip().split('\n'),ecs)

    fig = plt.figure()
    #title = f"{cdft_fnames[n][0]} with {cht} charges"
    #fig.suptitle(title)
    ax = fig.add_subplot(111, projection='3d', facecolor="black")
    # Set axis ratio
    ax.set_box_aspect((np.ptp(xs), np.ptp(ys), np.ptp(zs)))
    scatter = ax.scatter(xs,ys,zs,s=size,c=colors,alpha=1)
    for b,bnd in enumerate(bxs):
        ax.plot(bnd, bys[b], bzs[b], color='black', linewidth=2)
    for nv in nvs:
        #ax.plot(nv[0],nv[1],nv[2],color='red',lw=scale_vec)
        a = Arrow3D(nv[0],nv[1],nv[2],color='r',lw=scale_vec, arrowstyle="-|>",
                    mutation_scale=10)
        ax.add_artist(a)
    for l,lab in enumerate(labels):
        ax.text(xs[l],ys[l],zs[l], "{}".format(lab),color='white')
    plt.show()
    plt.close()
