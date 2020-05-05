import numpy as np
import matplotlib.pyplot as plt


d2y = 365.2425 #days in a year
r2au = 0.0046491 #solar radii to AU
au = 1.49597e13 #AU in cm
TWOPI = 2*np.pi

def user_rc(lw=1.5, fontsize=16, usetex=True):
    """Set plotting RC parameters"""
    plt.rc('lines', linewidth=lw)
    plt.rc('axes', lw=lw, labelsize=fontsize, titlesize=22)
    if usetex:
        plt.rc('text', usetex=True)
        plt.rc('font', **{'family': 'serif'})#, 'serif': ['Computer Modern Roman', 'Times']})
    else:
        plt.rc('text', usetex=False)
        plt.rc('font', **{'family': 'serif'})#, 'serif': ['Ubuntu'], 'monospace': ['Ubuntu Mono']})#['Computer Modern']})
    plt.rc('xtick', labelsize=max(18, int(fontsize*0.8)))
    plt.rc('ytick', labelsize=max(18, int(fontsize*0.8)))
    plt.rc('figure', titlesize=22, figsize=(10,8))
    return

user_rc()
  
#class Tatooine(object):
#    def __init__(self, obs_dur=1.0, m1=1.0, m2=1.0, 
#                 r1=1.0, r2=1.0):
#        self.obs_dur = obs_dur
#        self.m1 = m1
#        self.m2 = m2
#        self.r1 = r1
#        self.r2 = r2
        

def get_a(P, msum):
    """Compute semi-major axis given P (d) and sum of masses (Msun)"""
    return ((P/d2y)**2 * (msum))**(1/3.)

def get_P(a, msum):
    """Compute period given semi-major axis (AU) and sum of masses(Msun)"""
    return (a**3/msum)**(0.5) * d2y

def get_Pcrit(Pb, mA, mB, e):
    """Compute critical stability period given binary period, masses, and eccentricity"""
    a_b = ((Pb/d2y)**2 * (mA+mB))**(1./3.) #in AU
    mu = mB/(mA+mB)
    a_c = (1.6 + (5.1*e)+(-2.22*e**2) + 4.12*mu - 4.27*e*mu - 5.09*mu**2 + 4.61*e**2*mu**2)*a_b
    P_c = np.sqrt(a_c**3 / (mA+mB)) * d2y
    return P_c

def get_Tprec(Pp, m1, m2, ap, ab, mutual_inc):
    """Compute precession orbital period (d) in the circular binary limit 
    (Eqn 14)"""
    return abs(Pp*4./(3*np.cos(mutual_inc)) * (ap/ab)**2 * (m1+m2)**2 / (m1*m2))

def get_Tprec2(Pp, Pb, m1, m2, mutual_inc):
    """Compute precession orbital period (d) in the circular binary limit, 
    using only Ps and Ms."""
    return abs(4/3. * (Pp**7/Pb**4)**(1/3.) * (m1+m2)**2/(m1*m2) / np.cos(mutual_inc))
    
def get_Omegadot(Tprec):
    """Compute dOmega/dt (radians / day) of CBP. Omega is the 
    longitude of ascending node of the planet."""
    return TWOPI/Tprec

def get_f1(r1, ap, mutual_inc, ib):
    """Equation 11 of Li et al 2016; NOTE ib is |90-inclination of binary|; 
    ap in AU, r1 in Rsun, angles in rad"""
    cond = (-r1*r2au/ap + np.sin(abs(ib)))/np.sin(mutual_inc)
    if abs(cond)<1.:
        return np.arcsin(cond)
    elif cond<-1:
        return -np.pi/2.
    else:
        return np.pi/2.
        
def get_f1_vec(r1, ap, mutual_inc, ib):
    """Vectorized version of 
    Equation 11 of Li et al 2016; NOTE ib is |90-inclination of binary|; 
    ap in AU, r1 in Rsun, angles in rad"""
    cond = (-r1*r2au/ap + np.sin(abs(ib)))/np.sin(mutual_inc)
    results = np.zeros(cond.shape) + np.pi/2.
    state1 = (abs(cond)<1.)
    results[state1] = np.arcsin(cond[state1])
    state2 = (cond<-1)
    results[state2] = -np.pi/2.
    return results

def get_f2(r1, ap, mutual_inc, ib):
    """Equation 12 of Li et al 2016; NOTE ib is |90-inc of binary|;
    ap in AU, r1 in Rsun, angles in rad"""
    cond = (r1*r2au/ap + np.sin(abs(ib)))/np.sin(mutual_inc)
    if abs(cond)<1.:
        return np.arcsin(cond)
    elif cond<-1:
        return -np.pi/2.
    else:
        return np.pi/2.

def get_f2_vec(r1, ap, mutual_inc, ib):
    """Vectorized version of 
    Equation 12 of Li et al 2016; NOTE ib is |90-inc of binary|;
    ap in AU, r1 in Rsun, angles in rad"""
    cond = (r1*r2au/ap + np.sin(abs(ib)))/np.sin(mutual_inc)
    results = np.zeros(cond.shape) + np.pi/2.
    state1 = (abs(cond)<1.)
    results[state1] = np.arcsin(cond[state1])
    state2 = (cond<-1)
    results[state2] = -np.pi/2.
    return results

def get_Delta_Omega1(ab1, ap, f1, f2):
    """Range of Omega that allows planet-star orbit crossings 
    (Equation 10 of Li et al 2016)"""
    asin = np.arcsin(ab1/ap)
#     if type(ab1) == np.ndarray or type(ab1) == list:
#         other_min = np.zeros(len(ab1))+TWOPI
#     else:
#         other_min = TWOPI
#     if (np.pi/2.-f2 > asin) and (np.pi/2.+f1 > asin):
#         return np.min(np.vstack((2.*(f2-f1) + 4*asin, 
#                                  other_min)), axis=0)
#     elif f1+np.pi/2. > asin:
#         return np.nanmin(np.vstack((2.*(f2-f1) + 2*(np.pi/2.-f2) + 2*asin, 
#                                  other_min)), axis=0)
#     else:
#         return np.nanmin(np.vstack((2.*(f2-f1) + 2*(np.pi/2.-f2) + 2*(f1+np.pi/2.), 
#                                  other_min)), axis=0)
    if ((np.pi/2. - f2) > asin) and ((f1 + np.pi/2.) > asin):
        return min(2*(f2-f1) + 4*asin, TWOPI)
    elif (f1 + np.pi/2.) > asin:
        return min(2*(f2-f1) + 2*(np.pi/2. - f2) + 2*asin, TWOPI)
    else:
        return min(2*(f2-f1) + 2*(np.pi/2. - f2) + 2*(f1+np.pi/2.), TWOPI)

def get_Delta_Omega1_vec(ab1, ap, f1, f2):
    """Vectorized version of 
    Range of Omega that allows planet-star orbit crossings 
    (Equation 10 of Li et al 2016)"""
    asin = np.arcsin(ab1/ap)
    state1 = (np.pi/2.-f2 > asin) * (np.pi/2.+f1 > asin)
    result = np.zeros(state1.shape) + TWOPI
    other_min = np.zeros(state1.shape) + TWOPI
    result[state1] = (2.*(f2-f1) + 4*asin)[state1]
    state2 = (~state1) * (f1+np.pi/2. > asin)
    result[state2] = (2.*(f2-f1) + 2*(np.pi/2.-f2) + 2*asin)[state2]
    state3 = (~state1) * (~state2)
    result[state3] = (2.*(f2-f1) + 2*(np.pi/2.-f2) + 2*(f1+np.pi/2.))[state3]
    return np.nanmin(np.stack((result, other_min), axis=0), axis=0)


def get_del_Omega_prec(Omegadot, Tobs):
    """Total change in CBP Omega (long. of ascending node) due to precession 
    during obs. time period, where Omegadot is in rad/d, Tobs in d."""
    return Omegadot * Tobs

def get_del_Omega1(Delta_Omega1, del_Omega_prec, ab1, ap, f1, f2):
    """Total range in long of asc. node during obs. time period Tobs
    (Eqn 15 of Li et al 2016)"""
    if f2>f1:
        if (np.pi/2. - f2) > del_Omega_prec/2.:
            return Delta_Omega1 + 2.*del_Omega_prec
        elif (np.pi/2. - f2) > np.arcsin(ab1/ap):
            return Delta_Omega1 + del_Omega_prec + 2*(np.pi/2. - f2)
        else:
            return Delta_Omega1 + del_Omega_prec
    else:
        return 0.

def get_del_Omega1_vec(Delta_Omega1, del_Omega_prec, ab1, ap, f1, f2):
    """Vectorized version of 
    Total range in long of asc. node during obs. time period Tobs
    (Eqn 15 of Li et al 2016)"""
    good = (f2>f1)
    state1 = ((np.pi/2. - f2) > del_Omega_prec/2.)
    state2 = (~state1) * ((np.pi/2. - f2) > np.arcsin(ab1/ap))
    state3 = (~state1) * (~state2)
    result = np.zeros(state1.shape)
    result[good*state1] = (Delta_Omega1 + 2.*del_Omega_prec)[good*state1]
    result[good*state2] = (Delta_Omega1 + del_Omega_prec + 2*(np.pi/2. - f2))[good*state2]
    result[good*state3] = (Delta_Omega1 + del_Omega_prec)[good*state3]
    return result

def get_Pcross1(del_Omega1):
    """Probability to cross stellar orbit of single component in binary
    (Eqn 16 in Li et al 2016)"""
#     if type(del_Omega1) == np.ndarray or type(del_Omega1) == list:
#         other_min = np.zeros(len(del_Omega1)) + 2*np.pi
#     else:
#         other_min = 2*np.pi
#     return np.min(np.vstack((del_Omega1, other_min)), axis=0)/(2*np.pi)
    return min(del_Omega1, TWOPI)/TWOPI

def get_Pcross1_vec(del_Omega1):
    """Vectorized version of 
    Probability to cross stellar orbit of single component in binary
    (Eqn 16 in Li et al 2016)"""
    other_min = np.zeros(del_Omega1.shape)+ TWOPI
    return np.nanmin(np.stack((del_Omega1, other_min), axis=0), axis=0)/(TWOPI)

def get_crossing_duration(r1, vp, mutual_inc):
    """Time (d) it takes for planet to cross stellar orbit, vp in au/d, 
    r1 in rsun"""
    return np.pi/2. * r1 * r2au / (vp * np.sin(mutual_inc))

def get_vel(P, a):
    """Compute velocity of planet, P in d, a in AU"""
    return TWOPI*a/P

def get_dl_same(r1, vp, mutual_inc, v1):
    """Relative displacement (AU) of planet and star, when planet and star 
    moving in same direction, as planet and star both 
    are located toward the observer relative to system's COM. r1 in Rsun,
    vp in AU/d, v1 in AU/d. NOTE 
    this has nothing to do w/ primary or secondary component, but 
    rather location of a stellar component relative to planet and COM."""
    xing_dur = get_crossing_duration(r1, vp, mutual_inc)
    return xing_dur * abs(vp*np.cos(mutual_inc) - 2*v1/np.pi) + r1*r2au

def get_dl_opp(r1, vp, mutual_inc, v1):
    """Relative displacement of planet and star, when planet and star
    are moving in opposite directions, as star is on the other side of
    the system's COM. r1 in Rsun, vp & v1 in AU/d. NOTE 
    this has nothing to do w/ primary or secondary component, but 
    rather location of a stellar component relative to planet and COM."""
    xing_dur = get_crossing_duration(r1, vp, mutual_inc)
    ## XX TYPO in paper, should be +r1, not +2r1 (symmetric w/dl_same)
    return xing_dur * abs(vp*np.cos(mutual_inc) + 2*v1/np.pi) + r1 * r2au

def get_Pstar1(dl_same, dl_opp, ab1):
    """Probability that planet transits star 1, assuming that planet crosses
    stellar orbit w/projected size 2*ab1. All inputs in AU."""
    if (dl_same+dl_opp)/2. > 2*ab1:
        return 1.
    else:
        return (dl_same+dl_opp)/(4.*ab1)

def get_Pstar1_vec(dl_same, dl_opp, ab1):
    """Vectorized version of 
    Probability that planet transits star 1, assuming that planet crosses
    stellar orbit w/projected size 2*ab1. All inputs in AU."""
    result = (dl_same+dl_opp)/(4.*ab1)
    state1 = (dl_same+dl_opp)/2. > 2*ab1
    result[state1] = 1.0
    return result

def get_n1(Tobs, Pp, Delta_Omega1, del_Omega_prec, Omegadot, ab1, ap, f2):
    """Number of planet-stellar orbit crossings (Eqn 18 Li et al 2016)"""
    asin = np.arcsin(ab1/ap)
    numerator = Delta_Omega1/2. + del_Omega_prec
    #ratio of change in omega (relative to change due to precession)
    dOmega_ratio = (Delta_Omega1 / (Omegadot * Pp))
#     if type(Delta_Omega1) == np.ndarray or type(Delta_Omega1) == list:
#         other_min = np.zeros(len(del_Omega1)) + Tobs/Pp
#     else:
#         other_min = Tobs/Pp
#     if (np.pi/2. - f2 > asin) and (numerator > np.pi):
#         return np.min(np.vstack((other_min, 
#                                  numerator/np.pi * dOmega_ratio/2.)), axis=0)
#     elif np.pi/2. - f2 > asin:
#         return np.min(np.vstack((other_min, 
#                                  dOmega_ratio/2.)), axis=0)
#     elif Delta_Omega1 + del_Omega_prec > 2*np.pi:
#         return np.min(np.vstack((other_min,
#                                 (Delta_Omega1+del_Omega_prec)/(2*np.pi) * dOmega_ratio)), axis=0)
#     else:
#         return np.min(np.vstack((other_min, dOmega_ratio)), axis=0)
    if ((np.pi/2. - f2) > asin) and (numerator > np.pi):
        return min(Tobs/Pp, numerator/np.pi * dOmega_ratio/2.)
    elif (np.pi/2. - f2) > asin:
        return min(Tobs/Pp, dOmega_ratio/2.)
    elif (Delta_Omega1+del_Omega_prec) > TWOPI:
        return min(Tobs/Pp, (Delta_Omega1+del_Omega_prec)/TWOPI * dOmega_ratio)
    else:
        return min(Tobs/Pp, dOmega_ratio)
    
def get_n1_vec(Tobs, Pp, Delta_Omega1, del_Omega_prec, Omegadot, ab1, 
               ap, f2):
    """Number of planet-stellar orbit crossings (Eqn 18 Li et al 2016)"""
    asin = np.arcsin(ab1/ap)
    numerator = Delta_Omega1/2. + del_Omega_prec
    #ratio of change in omega (relative to change due to precession)
    dOmega_ratio = (Delta_Omega1 / (Omegadot * Pp))
    state1 = (np.pi/2. - f2 > asin) * (numerator > np.pi)
    
    other_min = np.zeros(state1.shape) + Tobs/Pp
    result = np.zeros(state1.shape) + TWOPI
    
    result[state1] = (numerator/np.pi * dOmega_ratio/2.)[state1]
    state2 = (~state1) * (np.pi/2. - f2 > asin)
    result[state2] = (dOmega_ratio/2.)[state2]
    state3 = (~state1) * (~state2) * (Delta_Omega1 + del_Omega_prec > 2*np.pi)
    result[state3] = ((Delta_Omega1+del_Omega_prec)/(2*np.pi) * dOmega_ratio)[state3]
    state4 = (~state1) * (~state2) * (~state3)
    result[state4] = dOmega_ratio[state4]
    return np.nanmin(np.stack((other_min, result), axis=0), axis=0)

def binomial(M, k, p0):
    """Computes the binomial function for M choose k events with 
    p0 probability (M=avg orbit crossings, k=N of transits, 
    p0=prob of transit"""
    assert k<=M, "M choose k means k<=M!"
    dummy = np.arange(k-1)+1
    if k == 0:
        Mchoosek=1.
    elif k == 1:
        Mchoosek=M*1.0
    else:
        numerator = np.prod(M-dummy)*M
        denom = np.prod(np.arange(1,k+1))
        Mchoosek = numerator/denom
    return Mchoosek * p0**k * (1.-p0)**(M-k)
    
def binomial_vec(M, k, p0):
    """Vectorized version of 
    Computes the binomial function for M choose k events with 
    p0 probability (M=avg orbit crossings, k=N of transits, 
    p0=prob of transit. M is array."""
    dummy = np.arange(k-1)+1
    if k == 0:
        Mchoosek=np.ones(M.shape)
    elif k == 1:
        Mchoosek=np.ones(M.shape)*M
    else:
        # product over dummy (the factorial) axis
        numerator = np.prod(M[:, np.newaxis]-dummy[np.newaxis,:], axis=-1)*M
        denom = np.prod(np.arange(1,k+1))
        Mchoosek = numerator/denom
    return Mchoosek * p0**k * (1.-p0)**(M-k)


def get_Pk(M, k, p0):
    """Computes total probability of i=0...k events (Prob of at least k 
    events = 1 - get_Pk())."""
    return np.sum([binomial(M, i, p0) for i in range(0, k)])
    
def get_Pk_vec(M, k, p0):
    """Vectorized version of 
    Computes total probability of i=0...k events (Prob of at least k 
    events = 1 - get_Pk())."""
    return np.sum([binomial_vec(M, i, p0) for i in range(0, k)], axis=0)    
    
def get_Ptran_star1(Pcross1, Pstar1, n1, fduty=1., k=3):
    """Computes the transit probability for *at least k transits* of planet 
    across one stellar component of the binary, given prob of stellar 
    orbit crossing (Pcross1), prob of stellar disk crossing (Pstar1), 
    and number of orbit crossings (n1), duty cycle of data (fduty)"""
    vectorized, arr_size = check_input_size([Pcross1, Pstar1, n1])
    if vectorized:
        P_k = get_Pk_vec(n1, k, fduty*Pstar1)
    else:
        P_k = get_Pk(n1, k, fduty*Pstar1)
    return Pcross1 * (1.-P_k)
    
def check_input_size(inputs):
    """Checks non-scalar elements in inputs have same length/size. Returns 
    False, 0 if all inputs scalar, True, N if input contains lists/arrays of size N."""
    N = len(inputs)
    whichone = [((type(ip) == list) or (type(ip) == np.ndarray)) for ip in inputs]
    if np.any(whichone):
        input_sizes = [len(inputs[wo]) for wo in np.arange(N)[whichone]]
        assert np.all(np.diff(input_sizes)==0), "Non-scalar inputs must have same size {}".format(input_sizes)
        return True, input_sizes[0]
    return False, 0
    
def get_Ptran_binary(Pcross1, Pcross2, Pstar1, Pstar2, n1, n2, fduty=1., k=3):
    """Computes transit probability for *at least k transits* of planet across
    both stellar components of binary, given prob of stellar orbit crossings 
    (Pcross1, Pcross2), prob of stellar disk crossings (Pstar1, Pstar2), and 
    number of orbit crossings (n1, n2), duty cycle of data (fduty)."""
    # if arrays present
    vectorized, arr_size = check_input_size([Pcross1, Pcross2, 
                                             Pstar1, Pstar2, n1, n2])
    if vectorized:
        switch = Pcross1>Pcross2
        if switch.sum()>0:
            Pstar1[switch], Pstar2[switch] = Pstar2[switch], Pstar1[switch]
            n1[switch], n2[switch] = n2[switch], n1[switch]
            Pcross1[switch], Pcross2[switch] = Pcross2[switch], Pcross1[switch]
        pk_star2 = get_Pk_vec(n2, k, fduty*Pstar2)
        pk_star1 = get_Pk_vec(n1, k, fduty*Pstar1) * Pcross1/Pcross2 + \
                                            (Pcross2-Pcross1)/Pcross2
        zerox = (Pcross1 == 0)*(Pcross2 == 0)
        result = Pcross2 * (1. - pk_star2 * pk_star1)
        result[zerox] = 0.
    # if only scalars
    else:
        if Pcross1==0 and Pcross2==0:
            return 0
            
        # if star1 sweeps out larger orbit area (b/c Q>1)
        if Pcross1>Pcross2:
#            print("boop_P")
            Pcross1, Pcross2 = Pcross2, Pcross1
            Pstar1, Pstar2 = Pstar2, Pstar1
            n1, n2 = n2, n1
            
        pk_star2 = get_Pk(n2, k, fduty*Pstar2)
        pk_star1 = get_Pk(n1, k, fduty*Pstar1)*Pcross1/Pcross2 + \
                                (Pcross2 - Pcross1)/Pcross2
                                
        result = Pcross2 * (1. - pk_star2 * pk_star1)
    return result    
 
def get_Qtran(m1, m2, r1, r2, Pb, Pp, i_b, di, Tobs, fduty, k, 
              components='both'):
    """Computes the total probability for at least k CBP transits during Tobs, 
    with observing duty fduty.
    
    Parameters:
    ----------
    m1, m2 : masses of primary and secondary in Msun
    r1, r2 : radii of primary and secondary in Rsun
    Pb, Pp : periods of binary and planet in days
    i_b : binary inclination (relative to observer) in rad
    di : mutual inclination between planet + binary in rad
    Tobs : total duration of observation in days
    fduty : duty cycle of observations (1 = no gaps, 0 = all missing cadences)
    k : minimum number of transits 
    components : compute for 'both' stellar components (default) or star '1' 
                    or star '2'
    
    Returns:
    -------
    Qt : total CBP transit probability
    Nt : average number of transits
    """
    
    # do calculations that don't need to differentiate btw. scalar/vector first
    
    # get system orbital parameters
    ab = get_a(Pb, m1+m2)
    ab1 = m2/(m1+m2) * ab
    ab2 = m1/(m1+m2) * ab
    ap = get_a(Pp, m1+m2)
    d_ib = np.pi/2. - i_b

    # precession timescale
    Tprec = get_Tprec(Pp, m1, m2, ap, ab, di)
    Omegadot = get_Omegadot(Tprec)
    del_Omega_prec = get_del_Omega_prec(Omegadot, Tobs)
    
    # orbit properties to compute prob to cross stellar disk
    vp = get_vel(Pp, ap)
    v1 = get_vel(Pb, ab1)
    v2 = get_vel(Pb, ab2)

    # relative displacement of star + planet traveling in same/opp directions 
    # during crossing
    dl_same1 = get_dl_same(r1, vp, di, v1)
    dl_opp1 = get_dl_opp(r1, vp, di, v1)
    dl_same2 = get_dl_same(r2, vp, di, v2)
    dl_opp2 = get_dl_opp(r2, vp, di, v2)
        
    # now differentiate...
    vectorized, arr_size = check_input_size([m1, m2, r1, r2, Pb, Pp, 
                                             i_b, di, Tobs, fduty, k])

    if vectorized:
        # vectorized
    
        # angles of precession wrt star 1
        fone1 = get_f1_vec(r1, ap, di, d_ib)
        ftwo1 = get_f2_vec(r1, ap, di, d_ib)
        Delta_Omega1 = get_Delta_Omega1_vec(ab1, ap, fone1, ftwo1)
    
        # angles of precession wrt star 2
        fone2 = get_f1_vec(r2, ap, di, d_ib)
        ftwo2 = get_f2_vec(r2, ap, di, d_ib)
        Delta_Omega2 = get_Delta_Omega1_vec(ab2, ap, fone2, ftwo2)
    
        # precession evolution relative to star 1 and 2
        del_Omega1 = get_del_Omega1_vec(Delta_Omega1, del_Omega_prec, ab1, 
                                    ap, fone1, ftwo1)
        del_Omega2 = get_del_Omega1_vec(Delta_Omega2, del_Omega_prec, ab2, 
                                    ap, fone2, ftwo2)
        
        
        # number of stellar orbit crossings
        n1 = get_n1_vec(Tobs, Pp, Delta_Omega1, del_Omega_prec, Omegadot, 
                    ab1, ap, ftwo1)
        n2 = get_n1_vec(Tobs, Pp, Delta_Omega2, del_Omega_prec, Omegadot, 
                    ab2, ap, ftwo2)
                    
        # prob to cross stellar disks
        Pstar1 = get_Pstar1_vec(dl_same1, dl_opp1, ab1)
        Pstar2 = get_Pstar1_vec(dl_same2, dl_opp2, ab2)

        # prob to cross stellar orbits
        Pcross1 = get_Pcross1_vec(del_Omega1)
        Pcross2 = get_Pcross1_vec(del_Omega2)

    else:

        # angles of precession wrt star 1
        fone1 = get_f1(r1, ap, di, d_ib)
        ftwo1 = get_f2(r1, ap, di, d_ib)
        Delta_Omega1 = get_Delta_Omega1(ab1, ap, fone1, ftwo1)
    
        # angles of precession wrt star 2
        fone2 = get_f1(r2, ap, di, d_ib)
        ftwo2 = get_f2(r2, ap, di, d_ib)
        Delta_Omega2 = get_Delta_Omega1(ab2, ap, fone2, ftwo2)
    
        # precession evolution relative to star 1 and 2
        del_Omega1 = get_del_Omega1(Delta_Omega1, del_Omega_prec, ab1, 
                                    ap, fone1, ftwo1)
        del_Omega2 = get_del_Omega1(Delta_Omega2, del_Omega_prec, ab2, 
                                    ap, fone2, ftwo2)
        
        
        # number of stellar orbit crossings
        n1 = get_n1(Tobs, Pp, Delta_Omega1, del_Omega_prec, Omegadot, 
                    ab1, ap, ftwo1)
        n2 = get_n1(Tobs, Pp, Delta_Omega2, del_Omega_prec, Omegadot, 
                    ab2, ap, ftwo2)
    
        # prob to cross stellar disks
        Pstar1 = get_Pstar1(dl_same1, dl_opp1, ab1)
        Pstar2 = get_Pstar1(dl_same2, dl_opp2, ab2)
    
        # prob to cross stellar orbits
        Pcross1 = get_Pcross1(del_Omega1)
        Pcross2 = get_Pcross1(del_Omega2)
    
#    print("Pc1, Pc2", Pcross1, Pcross2)
    components = components.lower()
    if components == 'both':
        Qt = get_Ptran_binary(Pcross1, Pcross2, Pstar1, Pstar2, n1, n2, 
                          fduty=fduty, k=k)
        Nt = get_Ntransit_avg(Pcross1, Pcross2, Pstar1, Pstar2, n1, n2)
    elif components == '1' or components == 'a':
        Qt = get_Ptran_star1(Pcross1, Pstar1, n1, fduty=fduty, k=k)
        Nt = n1*Pstar1
    elif components == '2' or components == 'b':
        Qt = get_Ptran_star1(Pcross2, Pstar2, n2, fduty=fduty, k=k)
        Nt = n2*Pstar2
    else:
        print("Options are 'both', '1', or '2'... defaulting to 'both'.")
        Qt = get_Ptran_binary(Pcross1, Pcross2, Pstar1, Pstar2, n1, n2, 
                          fduty=fduty, k=k)
        Nt = get_Ntransit_avg(Pcross1, Pcross2, Pstar1, Pstar2, n1, n2)
    return Qt, Nt

def get_Ptran_star1_at_least_once(Pcross1, Pstar1, n1):
    """Probability that planet will transit star m1 at least once, assuming 
    orbit crossings are independent from each other (Eqn 19 Li et al 2016)"""
    return Pcross1 * (1 - (1-Pstar1)**n1)

def get_Ptran_binary_at_least_once(Pcross1, Pcross2, Pstar1, Pstar2, n1, n2):
    """Probability that planet will transit two stars at least once, assuming
    that m1>m2"""
    if Pcross1==0 and Pcross2==0:
        return 0
    if Pcross1>Pcross2:
#        print("boop_P")
        Pcross1, Pcross2 = Pcross2, Pcross1
        Pstar1, Pstar2 = Pstar2, Pstar1
        n1, n2 = n2, n1
    return Pcross2 * (1-(1-Pstar2)**n2 * \
                      (Pcross1/Pcross2 * (1-Pstar1)**n1 + (Pcross2-Pcross1)/Pcross2))
    
def get_Ntransit_avg(Pcross1, Pcross2, Pstar1, Pstar2, n1, n2):
    """Average number of transits given that a system transits, assuming 
    that m1>m2."""
    vectorized, arr_size = check_input_size([Pcross1, Pcross2, Pstar1, 
                                             Pstar2, n1, n2])
    if vectorized:
        switch = Pcross1>Pcross2
        if switch.sum()>0:
            Pstar1[switch], Pstar2[switch] = Pstar2[switch], Pstar1[switch]
            n1[switch], n2[switch] = n2[switch], n1[switch]
            Pcross1[switch], Pcross2[switch] = Pcross2[switch], Pcross1[switch]
        zerox = (Pcross1==0.)*(Pcross2==0.)
        result = n2*Pstar2 + n1*Pstar1*Pcross1/Pcross2
        result[zerox] = 0.
    else:
        if Pcross1==0 and Pcross2==0:
            return 0
        if Pcross1>Pcross2:
#            print("boop_N")
            Pcross1, Pcross2 = Pcross2, Pcross1
            Pstar1, Pstar2 = Pstar2, Pstar1
            n1, n2 = n2, n1
        result = n2*Pstar2 + n1*Pstar1*Pcross1/Pcross2
        
    return result

def get_ipc1(ib, ab1, r1, ap):
    """Critical planet inclination above which transits cannot occur b/c planet
    does not cross the stellar orbits. ab1 & ap in AU, r1 in Rsun"""
    return abs(ib) - np.arcsin((ab1*np.sin(ib) + r1*r2au)/ap)

def get_Ptran_winc(Ptran, mutual_inc, ipc1, ipc2):
    """Probability that planet will transit two stars NOT directly edge-on 
    at least once."""
    vectorized, arr_size = check_input_size([Ptran, mutual_inc, ipc1, ipc2])
    if vectorized:
        bad = (mutual_inc <= np.nanmin([ipc1, ipc2], axis=0))
        result = Ptran*1.0
        result[bad] = 0.
    else:
        if mutual_inc <= min(ipc1, ipc2):
            return 0.
    return Ptran
    
    
def unit_test_get_all(m1, m2, r1, r2, Pb, Pp, i_b, di, Tobs, fduty, k):
    # do calculations that don't need to differentiate btw. scalar/vector first
    
    # get system orbital parameters
    ab = get_a(Pb, m1+m2)
    ab1 = m2/(m1+m2) * ab
    ab2 = m1/(m1+m2) * ab
    ap = get_a(Pp, m1+m2)
    d_ib = np.pi/2. - i_b

    # precession timescale
    Tprec = get_Tprec(Pp, m1, m2, ap, ab, di)
    Omegadot = get_Omegadot(Tprec)
    del_Omega_prec = get_del_Omega_prec(Omegadot, Tobs)
    
    # orbit properties to compute prob to cross stellar disk
    vp = get_vel(Pp, ap)
    v1 = get_vel(Pb, ab1)
    v2 = get_vel(Pb, ab2)

    # relative displacement of star + planet traveling in same/opp directions 
    # during crossing
    dl_same1 = get_dl_same(r1, vp, di, v1)
    dl_opp1 = get_dl_opp(r1, vp, di, v1)
    dl_same2 = get_dl_same(r2, vp, di, v2)
    dl_opp2 = get_dl_opp(r2, vp, di, v2)
        
    # now differentiate...
    vectorized, arr_size = check_input_size([m1, m2, r1, r2, Pb, Pp, 
                                             i_b, di, Tobs, fduty, k])

    if vectorized:
        # vectorized
    
        # angles of precession wrt star 1
        fone1 = get_f1_vec(r1, ap, di, d_ib)
        ftwo1 = get_f2_vec(r1, ap, di, d_ib)
        Delta_Omega1 = get_Delta_Omega1_vec(ab1, ap, fone1, ftwo1)
    
        # angles of precession wrt star 2
        fone2 = get_f1_vec(r2, ap, di, d_ib)
        ftwo2 = get_f2_vec(r2, ap, di, d_ib)
        Delta_Omega2 = get_Delta_Omega1_vec(ab2, ap, fone2, ftwo2)
    
        # precession evolution relative to star 1 and 2
        del_Omega1 = get_del_Omega1_vec(Delta_Omega1, del_Omega_prec, ab1, 
                                    ap, fone1, ftwo1)
        del_Omega2 = get_del_Omega1_vec(Delta_Omega2, del_Omega_prec, ab2, 
                                    ap, fone2, ftwo2)
        
        
        # number of stellar orbit crossings
        n1 = get_n1_vec(Tobs, Pp, Delta_Omega1, del_Omega_prec, Omegadot, 
                    ab1, ap, ftwo1)
        n2 = get_n1_vec(Tobs, Pp, Delta_Omega2, del_Omega_prec, Omegadot, 
                    ab2, ap, ftwo2)
                    
        # prob to cross stellar disks
        Pstar1 = get_Pstar1_vec(dl_same1, dl_opp1, ab1)
        Pstar2 = get_Pstar1_vec(dl_same2, dl_opp2, ab2)

        # prob to cross stellar orbits
        Pcross1 = get_Pcross1_vec(del_Omega1)
        Pcross2 = get_Pcross1_vec(del_Omega2)
    else:
        # angles of precession wrt star 1
        fone1 = get_f1(r1, ap, di, d_ib)
        ftwo1 = get_f2(r1, ap, di, d_ib)
        Delta_Omega1 = get_Delta_Omega1(ab1, ap, fone1, ftwo1)
    
        # angles of precession wrt star 2
        fone2 = get_f1(r2, ap, di, d_ib)
        ftwo2 = get_f2(r2, ap, di, d_ib)
        Delta_Omega2 = get_Delta_Omega1(ab2, ap, fone2, ftwo2)
    
        # precession evolution relative to star 1 and 2
        del_Omega1 = get_del_Omega1(Delta_Omega1, del_Omega_prec, ab1, 
                                    ap, fone1, ftwo1)
        del_Omega2 = get_del_Omega1(Delta_Omega2, del_Omega_prec, ab2, 
                                    ap, fone2, ftwo2)
        
        
        # number of stellar orbit crossings
        n1 = get_n1(Tobs, Pp, Delta_Omega1, del_Omega_prec, Omegadot, 
                    ab1, ap, ftwo1)
        n2 = get_n1(Tobs, Pp, Delta_Omega2, del_Omega_prec, Omegadot, 
                    ab2, ap, ftwo2)
    
        # prob to cross stellar disks
        Pstar1 = get_Pstar1(dl_same1, dl_opp1, ab1)
        Pstar2 = get_Pstar1(dl_same2, dl_opp2, ab2)
    
        # prob to cross stellar orbits
        Pcross1 = get_Pcross1(del_Omega1)
        Pcross2 = get_Pcross1(del_Omega2)
        
    return fone1, ftwo1, Delta_Omega1, fone2, ftwo2, Delta_Omega2, \
            del_Omega1, del_Omega2, n1, n2, Pstar1, Pstar2, \
            Pcross1, Pcross2
            
def unit_test_check_all(m1, m2, r1, r2, Pb, Pp, i_b, di, Tobs, fduty, k, 
                   components='both', tol=1e-6):
    ball1 = unit_test_get_all(m1, m2, r1, r2, Pb, Pp, i_b, di, Tobs, fduty, k)
    fone1, ftwo1, Delta_Omega1, fone2, ftwo2, Delta_Omega2, \
            del_Omega1, del_Omega2, n1, n2, Pstar1, Pstar2, \
            Pcross1, Pcross2 = ball1
    _fone1, _ftwo1, _Delta_Omega1, _fone2, _ftwo2, _Delta_Omega2, \
        _del_Omega1, _del_Omega2, _n1, _n2, _Pstar1, _Pstar2, \
        _Pcross1, _Pcross2 = fone1*0., ftwo1*0., Delta_Omega1*0., \
                            fone2*0., ftwo2*0., Delta_Omega2*0., \
                            del_Omega1*0., del_Omega2*0., n1*0., n2*0., \
                            Pstar1*0., Pstar2*0., Pcross1*0., Pcross2*0.
    for ii in range(len(fone1)):
        _fone1[ii], _ftwo1[ii], _Delta_Omega1[ii], \
            _fone2[ii], _ftwo2[ii], _Delta_Omega2[ii], \
            _del_Omega1[ii], _del_Omega2[ii], _n1[ii], _n2[ii], \
            _Pstar1[ii], _Pstar2[ii], \
            _Pcross1[ii], _Pcross2[ii] = unit_test_get_all(m1, m2, r1, r2, 
                                                            Pb, Pp, i_b, di[ii], 
                                                            Tobs, fduty, k)
    passtest=True
    if (abs(_fone1-fone1)>tol).sum()>0:
        print("get_f1; star1 failed")
        passtest*=False
    if (abs(_fone2-fone2)>tol).sum()>0:
        print("get_f1; star2 failed")
        passtest*=False
    if (abs(_ftwo1-ftwo1)>tol).sum()>0:
        print("get_f2; star1 failed")
        passtest*=False
    if (abs(_ftwo2-ftwo2)>tol).sum()>0:
        print("get_f2; star2 failed")
        passtest*=False
    if (abs(_Delta_Omega1-Delta_Omega1)>tol).sum()>0:
        print("get_Delta_Omega1 failed")
        passtest*=False
    if (abs(_Delta_Omega2-Delta_Omega2)>tol).sum()>0:
        print("get_Delta_Omega2 failed")
        passtest*=False
    if (abs(_del_Omega1-del_Omega1)>tol).sum()>0:
        print("get_del_Omega1; star1 failed")
        passtest*=False
    if (abs(_del_Omega2-del_Omega2)>tol).sum()>0:
        print("get_del_Omega1; star2 failed")
        passtest*=False
    if (abs(_n1-n1)>tol).sum()>0:
        print("get_n1 failed")    
        passtest*=False
    if (abs(_n2-n2)>tol).sum()>0:
        print("get_n2 failed")   
        passtest*=False
    if (abs(_Pstar1-Pstar1)>tol).sum()>0:
        print("get_Pstar1 failed")    
        passtest*=False
    if (abs(_Pstar2-Pstar2)>tol).sum()>0:
        print("get_Pstar2 failed")  
        passtest*=False
    if (abs(_Pcross1-Pcross1)>tol).sum()>0:
        print("get_Pcross1 failed")    
        passtest*=False
    if (abs(_Pcross2-Pcross2)>tol).sum()>0:
        print("get_Pcross2 failed") 
        passtest*=False
    return passtest