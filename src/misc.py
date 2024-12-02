import math
import numpy as np

def format_numbers_as_ranges(numbers):
    """Format numbers for constraints.inp file"""
    ranges = []
    current_range = [numbers[0], numbers[0]]

    for i in range(1, len(numbers)):
        if numbers[i] == current_range[1] + 1:
            current_range[1] = numbers[i]
        else:
            ranges.append(current_range)
            current_range = [numbers[i], numbers[i]]

    ranges.append(current_range)

    formatted_ranges = []
    for start, end in ranges:
        if start == end:
            formatted_ranges.append(str(start))
        else:
            formatted_ranges.append(f"{start}-{end}")

    return ', '.join(formatted_ranges)

def calculate_angle(point1, point2, point3):

    vector1 = [p1 - p2 for p1, p2 in zip(point1, point2)]
    vector2 = [p3 - p2 for p3, p2 in zip(point3, point2)]

    dot_product = sum(v1 * v2 for v1, v2 in zip(vector1, vector2))

    magnitude1 = math.sqrt(sum(v ** 2 for v in vector1))
    magnitude2 = math.sqrt(sum(v ** 2 for v in vector2))

    cosine_angle = dot_product / (magnitude1 * magnitude2)

    angle_radians = math.acos(cosine_angle)
    angle_degrees = math.degrees(angle_radians)

    return angle_degrees

# Calculate dihedral angle between four points
# https://gist.github.com/hypnopump/30d6bfdb8358d3a57d010c9a501fda56
def calculate_dihedral_angle(point1, point2, point3, point4):

    point1, point2, point3, point4 = np.array(point1), np.array(point2),\
                                     np.array(point3), np.array(point4)

    a1 = point2 - point1
    a2 = point3 - point2
    a3 = point4 - point3

    v1 = np.cross(a1, a2)
    v1 = v1 / (v1 * v1).sum(-1)**0.5
    v2 = np.cross(a2, a3)
    v2 = v2 / (v2 * v2).sum(-1)**0.5
    porm = np.sign((v1 * a3).sum(-1))
    angle_rad = np.arccos((v1*v2).sum(-1) / ((v1**2).sum(-1)\
                * (v2**2).sum(-1))**0.5)
    if not porm == 0:
        angle_rad = angle_rad * porm

    angle_deg = np.degrees(angle_rad)

    return angle_deg

def calc_cdft_gds(ofiles,r=6):
    cdft_gds = []
    energies = []
    for of in ofiles:
        energies.append((of.charge,of.last_SCF))
    energies.sort(key=lambda k: k[0]) 
    ionpot = round(energies[2][1]-energies[1][1],r)
    cdft_gds.append(ionpot)
    elecaf = round(energies[1][1]-energies[0][1],r)
    cdft_gds.append(elecaf)
    mu = round(-0.5*(ionpot+elecaf),r)
    cdft_gds.append(mu)
    mup = -elecaf
    cdft_gds.append(mup)
    mum = -ionpot
    cdft_gds.append(mum)
    eta = round(ionpot - elecaf,r)
    cdft_gds.append(eta)
    omega = round(mu**2/(2*eta),r)
    cdft_gds.append(omega)
    return cdft_gds
    

def calc_cdft_funcs(ofiles,cdft_gds,r=6):
    cdft_fnames = [(r'f$^{+}$','f+'),(r'f$^{-}$','f-'),('\u0394f','delta_f'),
                   (r's$^{+}$','s+'),(r's$^{-}$','s-'),('\u0394s','delta_s'),
                   (u'\u0394\u03C1$_{elec}$','delta_rho_elec'),
                   (u'\u0394\u03C1$_{nuc}$','delta_rho_nuc')]

    charges = []
    for of in ofiles:
        charges.append((of.charge,of.charges))
    charges.sort(key=lambda k: k[0]) 

    chtypes = list()
    for cht in ofiles[1].charges:
        chtypes.append(cht)
    cdfts = list()
    for f,func in enumerate(cdft_fnames):
        cdfts.append({})
        for cht in chtypes:
            cdfts[f][cht] = list()
    for cht in chtypes:
        for n,at in enumerate(charges[1][1][cht]):
            fpres = round(at[1] - charges[0][1][cht][n][1],r)
            fmres = round(charges[2][1][cht][n][1] - at[1],r)
            fdres = round(fpres - fmres,6)
            mup = cdft_gds[3]
            mum = cdft_gds[4]
            eta = cdft_gds[5]
            dreres = round(-(mup/eta)*fpres + 0.5*(mup/eta)**2*fdres,r)
            drnres = round((mum/eta)*fmres + 0.5*(mum/eta)**2*fdres,r)
            for f,func in enumerate(cdft_fnames):
                #match f:
                #    case 0:
                #        cdfts[f][cht].append(fpres)
                #    case 1:
                #        cdfts[f][cht].append(fmres)
                #    case 2:
                #        cdfts[f][cht].append(fdres)
                #    case 3:
                #        cdfts[f][cht].append(round(fpres/eta,r))
                #    case 4:
                #        cdfts[f][cht].append(round(fmres/eta,r))
                #    case 5:
                #        cdfts[f][cht].append(round(fdres/eta**2,r))
                #    case 6:
                #        cdfts[f][cht].append(dreres)
                #    case 7:
                #        cdfts[f][cht].append(drnres)
                if f == 0:
                    cdfts[f][cht].append(fpres)
                elif f == 1:
                    cdfts[f][cht].append(fmres)
                elif f == 2:
                    cdfts[f][cht].append(fdres)
                elif f == 3:
                    cdfts[f][cht].append(round(fpres/eta,r))
                elif f == 4:
                    cdfts[f][cht].append(round(fmres/eta,r))
                elif f == 5:
                    cdfts[f][cht].append(round(fdres/eta**2,r))
                elif f == 6:
                    cdfts[f][cht].append(dreres)
                elif f == 7:
                    cdfts[f][cht].append(drnres)
    return cdfts

def calcelectrophil(NPA):
    NPA_C = {}
    NPA_E = []
    Local_E = []
    for i in range(1,4):
        NPA_C[i] = [item.split()[2] for item in NPA[i-1][1]]
        NPA_C[i] = [float(item) for item in NPA_C[i]]

        NPA_E.append(NPA[i-1][0])       
    
    NPA_E = [float(item) for item in NPA_E]

    for i in range(0,len(NPA_C[1])):
        if NPA[0][1][i].split()[0] != 'H':
            Elec_local = ((NPA_E[1]-NPA_E[0])*(NPA_C[2][i]-NPA_C[1][i]))\
                         /(NPA_E[2]-2*NPA_E[1]+NPA_E[0])\
                         +math.pow((NPA_E[1]-NPA_E[0])/(NPA_E[2]-2*NPA_E[1]
                                   +NPA_E[0]),2)\
                         *((2*NPA_C[2][i]-NPA_C[1][i]-NPA_C[3][i])/2)
            Atom = NPA[0][1][i].split()[0]
            Number = NPA[0][1][i].split()[1]
            Elec_local = Atom, Number, Elec_local
            Local_E.append(Elec_local)
    
    Local_E.sort(reverse=True, key=lambda item: item[2])
    return Local_E

def compnv(refnv,compnv):
    mode = 'scal_prod'
    if mode == 'norm_thr':
        res = True
        thr1 = 0.6
        rnorms = list()
        cnorms = list()
        for a,at in enumerate(refnv):
            rnorms.append(np.linalg.norm(at[1]))
            cnorms.append(np.linalg.norm(compnv[a][1]))
        print(rnorms,cnorms)
        rmax = sorted(rnorms)[-1]
        cmax = sorted(cnorms)[-1]
        rats = list()
        cats = list()
        for n,norm in enumerate(rnorms):
            if norm > thr1*rmax:
                rats.append(n)
            if cnorms[n] > thr1*cmax:
                cats.append(n)
        if cats != rats:
            res = False
        print(rats,cats)
        print(res)
    elif mode == 'norm_top':
        res = True
        thr1 = 0.5
        thr2 = 0.7
        rnorms = dict()
        cnorms = dict()
        rnormtot = 0
        cnormtot = 0
        for a,at in enumerate(refnv):
            rnorms[a] = np.linalg.norm(at[1])
            cnorms[a] = np.linalg.norm(compnv[a][1])
            rnormtot += np.linalg.norm(at[1])
            cnormtot += np.linalg.norm(compnv[a][1])
        srnorms = {k: v/rnormtot for k, v in sorted(rnorms.items(), 
                   key=lambda item: item[1], reverse = True)}
        scnorms = {k: v/cnormtot for k, v in sorted(cnorms.items(), 
                   key=lambda item: item[1], reverse = True)}
        prop = 0
        rmainat = list()
        for at in srnorms:
            rmainat.append(at)
            prop += srnorms[at]
            if prop > thr1:
                break
        prop = 0
        cmainat = list()
        for at in scnorms:
            cmainat.append(at)
            prop += scnorms[at]
            if prop > thr2:
                break
    elif mode == "scal_prod":
        thr = 0.5
        res = True
        reftnv = list()
        comptnv = list()
        for a,at in enumerate(refnv):
            for v in at[1]:
                reftnv.append(v)
            for v in compnv[a][1]:
                comptnv.append(v)
        sp = np.dot(reftnv,comptnv)
        print(sp)
        if sp < thr:
            res = False
    return res
        
def checknv(compnv,activ_ats):
    res = True
    thr2 = 0.5
    lim = 3
    cnorms = dict()
    cnormtot = 0
    for a,at in enumerate(compnv):
        cnorms[a+1] = np.linalg.norm(at[1])
        cnormtot += np.linalg.norm(at[1])
    scnorms = {k: v/cnormtot for k, v in sorted(cnorms.items(), 
               key=lambda item: item[1], reverse = True)}
    #prop = 0
    #cmainat = list()
    #for at in scnorms:
    #    cmainat.append(at)
    #    prop += scnorms[at]
    #    if prop > thr2:
    #        break
    #for actat in activ_ats:
    #    if actat not in cmainat:
    #        res = False
    for actat in activ_ats:
        if actat not in list(scnorms.keys())[:lim]:
            res = False
            break
    return res
    
