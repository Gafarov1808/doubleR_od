from sqlalchemy import select, desc
from datetime import datetime
import numpy as np

from kiamdb.meas import TrackletMeta, Tracklet, Station, SessionMeas
import pyorbs
import doubleR_module

r1, r2 = 39, 42


def get_data():
    
    t_list = []
    ra_list = []
    de_list = []
    rotation_mats = []
    r_xyz = []

    sq1 = select(TrackletMeta.unit_id, TrackletMeta.ip_id).where(
                    TrackletMeta.t_start > datetime(2026, 3, 1), 
                    TrackletMeta.t_end < datetime(2026, 3, 12)).order_by(
                    desc(TrackletMeta.t_end - TrackletMeta.t_start)).limit(1)
    
    with SessionMeas() as session:
        res = session.execute(sq1).all() 
    unit_id, ip_id = res[0].unit_id, res[0].ip_id

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
        #print(rot_mat.T @ wgs* 1e-6)
        r_xyz.append(rot_mat.T @ wgs * 1e-6)

    return t_list[0], ra_list[0], de_list[0], r_xyz 

def get_JD(t1, t2, t3):
    date_t1 = float(t1.strftime('%Y%m%d'))
    time_t1 = float(t1.strftime('%H%M%S.%f'))
    JD1, _ = pyorbs.bal.datetime_jdt(date_t1, time_t1)
    date_t2 = float(t2.strftime('%Y%m%d'))
    time_t2 = float(t2.strftime('%H%M%S.%f'))
    JD2, _ = pyorbs.bal.datetime_jdt(date_t2, time_t2)
    date_t3 = float(t3.strftime('%Y%m%d'))
    time_t3 = float(t3.strftime('%H%M%S.%f'))
    JD3, _ = pyorbs.bal.datetime_jdt(date_t3, time_t3)

    return JD1, JD2, JD3

def get_L_mat(ra, dec):
    L = []
    L.append(np.cos(dec)*np.cos(ra))
    L.append(np.cos(dec)*np.sin(ra))
    L.append(np.sin(dec))
    return L

def main():

    times, ra, dec, xyz = get_data()

    t1, t2, t3 = times[0], times[len(times)//2], times[-1]
    JD1, JD2, JD3 = get_JD(t1, t2, t3)

    ra1, ra2, ra3 = ra[0], ra[len(times)//2], ra[-1]
    dec1, dec2, dec3 = dec[0], dec[len(times)//2], dec[-1]
    L1 = get_L_mat(ra1, dec1)
    L2 = get_L_mat(ra2, dec2)
    L3 = get_L_mat(ra3, dec3)

    r_st1, r_st2, r_st3 = xyz[0], xyz[len(times)//2], xyz[-1]

    determ = doubleR_module.DoubleRIteration(L1, L2, L3, JD1, JD2, JD3, r_st1, r_st2, r_st3)
    determ.solver(r1, r2)
    state_v = determ.get_state()
    print(state_v)

if __name__ == "__main__":
    main()
