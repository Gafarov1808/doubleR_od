from sqlalchemy import select, desc
from datetime import datetime
import numpy as np

from kiamdb.meas import TrackletMeta, Tracklet, Station, SessionMeas
import pyorbs
import doubleR_module

r1, r2 = 25, 50

def get_data():
    
    t_list = []
    ra_list = []
    de_list = []
    rotation_mats = []
    r_xyz = []
    n = 4

    sq1 = select(TrackletMeta.unit_id, TrackletMeta.ip_id).where(
                    TrackletMeta.t_start > datetime(2026, 3, 1), 
                    TrackletMeta.t_end < datetime(2026, 3, 12)).order_by(
                    desc(TrackletMeta.t_end - TrackletMeta.t_start)).limit(n)
    with SessionMeas() as session:
        res1 = session.execute(sq1).all() 
    unit_id, ip_id = res1[n-1].unit_id, res1[n-1].ip_id

    sq2 = select(Tracklet.time, Tracklet.ra_rad, Tracklet.de_rad).where(Tracklet.unit_id == unit_id)
    with SessionMeas() as session:
        results = session.execute(sq2).all()
        for time, ra_rad, de_rad in results:
            t_list.append(time)
            ra_list.append(ra_rad)
            de_list.append(de_rad)

    sq3 = select(Station.xyz_m).where(Station.ip_id == ip_id)
    with SessionMeas() as session:
        res = session.execute(sq3).all()
    wgs = res[0][0]

    for t in t_list[0]:
        rotation_mats.append(pyorbs.auxilary.eme2itrf(t))

    for rot_mat in rotation_mats:
        r_xyz.append(rot_mat.T @ wgs * 1e-6)

    return t_list[0], ra_list[0], de_list[0], r_xyz 

def get_L_mat(ra, dec):
    L = []
    L.append(np.cos(dec)*np.cos(ra))
    L.append(np.cos(dec)*np.sin(ra))
    L.append(np.sin(dec))
    return L

def main():

    times, ra, dec, xyz = get_data()
    second_pos = len(times) // 2

    t1, t2, t3 = times[0], times[second_pos], times[-1]
    ra1, ra2, ra3 = ra[0], ra[second_pos], ra[-1]
    dec1, dec2, dec3 = dec[0], dec[second_pos], dec[-1]
    L1 = get_L_mat(ra1, dec1)
    L2 = get_L_mat(ra2, dec2)
    L3 = get_L_mat(ra3, dec3)

    r_st1, r_st2, r_st3 = xyz[0], xyz[second_pos], xyz[-1]
    tau1 = (t1-t2).total_seconds() * 1e-3
    tau3 = (t3-t2).total_seconds() * 1e-3

    determ = doubleR_module.DoubleRIteration(L1, L2, L3, tau1, tau3, r_st1, r_st2, r_st3)
    determ.solver(r1, r2)
    state_v = determ.get_state()
    elements = determ.get_elements()
    print(f'state = {state_v}')
    print(f'elements = {elements}')

if __name__ == "__main__":
    main()
