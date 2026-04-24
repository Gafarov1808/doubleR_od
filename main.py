#!/usr/bin/env python3
from sqlalchemy import select, desc
from datetime import datetime, timedelta
import numpy as np

from kiamdb.meas import TrackletMeta, Tracklet, Station, SessionMeas
from kiamdb.orbits import OrbitSolution, SessionOrbits
import pyorbs
import ODdoubleR

def get_data(n):
    
    t_list = []
    ra_list, de_list = [], []
    rotation_mats = []
    r_xyz = []

    sq1 = select(TrackletMeta.unit_id, TrackletMeta.ip_id, TrackletMeta.obj_id,
                 TrackletMeta.t_start, TrackletMeta.t_end).where(
                    TrackletMeta.t_start > datetime(2026, 2, 28), 
                    TrackletMeta.t_end < datetime(2026, 4, 24)).order_by(
                    desc(TrackletMeta.t_end - TrackletMeta.t_start)).limit(n)
    with SessionMeas() as session:
        res1 = session.execute(sq1).all() 
    unit_id, ip_id = res1[n-1].unit_id, res1[n-1].ip_id
    obj, start, end = res1[n-1].obj_id, res1[n-1].t_start, res1[n-1].t_end

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

    sq4 = select(OrbitSolution.state, OrbitSolution.epoch).where(
        OrbitSolution.obj_id == obj, OrbitSolution.time_obtained > start - timedelta(days=5), 
        OrbitSolution.time_obtained < end + timedelta(days=3))
    
    with SessionOrbits() as session:
        res_orb = session.execute(sq4).all()

    if res_orb == []:
        raise ValueError("Не удалось найти истинную орбиту в таблице solutions.")
    
    state, epoch = np.array(res_orb[0].state), res_orb[0].epoch
    return t_list[0], ra_list[0], de_list[0], r_xyz, state, epoch

def main():
    amount = 100
    k = 0
    for count in range(1, amount):
        try:
            times, ra, dec, xyz, v, t = get_data(count)
        except ValueError as e:
            print(f"{count}) {e}")
            continue
        central = len(times) - 5

        t1, t2, t3 = times[0], times[central], times[-1]
        Ra, Dec = [ra[0], ra[central], ra[-1]], [dec[0], dec[central], dec[-1]]
        r_st1, r_st2, r_st3 = xyz[0], xyz[central], xyz[-1]
        tau1, tau3 = (t1-t2).total_seconds() * 1e-3, (t3-t2).total_seconds() * 1e-3

        RIt = ODdoubleR.DoubleRIteration(Ra, Dec, r_st1, r_st2, r_st3, tau1, tau3)
        RIt.solver(2 + 6.371, 50 + 6.371)
        state = RIt.get_state()

        orb = pyorbs.pyorbs.orbit()
        orb.state_v, orb.time = v, pyorbs.pyorbs.ephem_time(t)
        orb.setup_parameters()
        try:
            orb.move(t2)
        except pyorbs.exceptions.ReentryException as e:
            print(f'{count}) {e}')
            continue

        if np.linalg.norm(state) < 1.1 * np.linalg.norm(orb.state_v):
            k += 1
    print(k/(amount-1))

if __name__ == "__main__":
    main()
