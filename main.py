#!/usr/bin/env python3
from sqlalchemy import select, desc
from datetime import datetime

from kiamdb.meas import TrackletMeta, Tracklet, Station, SessionMeas
import pyorbs
import ODdoubleR

def get_data(n):
    
    t_list = []
    ra_list = []
    de_list = []
    rotation_mats = []
    r_xyz = []

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

def main():

    for count in range(1, 200):
        times, ra, dec, xyz = get_data(count)
        central = len(times) - 15

        t1, t2, t3 = times[0], times[central], times[-1]
        Ra, Dec = [ra[0], ra[central], ra[-1]], [dec[0], dec[central], dec[-1]]
        r_st1, r_st2, r_st3 = xyz[0], xyz[central], xyz[-1]
        tau1, tau3 = (t1-t2).total_seconds() * 1e-3, (t3-t2).total_seconds() * 1e-3

        RIt = ODdoubleR.DoubleRIteration(Ra, Dec, r_st1, r_st2, r_st3, tau1, tau3)
        RIt.solver(2 + 6.371, 50 + 6.371)
        RIt.print_elements()

if __name__ == "__main__":
    main()
