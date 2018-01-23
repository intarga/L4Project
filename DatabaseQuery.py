import eagleSqlTools as sql
import numpy as np
import matplotlib.pyplot as plt

mySims = np.array([('RefL0100N1504', 100.)])
con = sql.connect('aabraham', password='LM277HBz')

for sim_name, sim_size in mySims:
    print sim_name

    myQuery = 'SELECT \
                    AP.ApertureSize, \
                    AP.Mass_Star, \
                    AP.SFR, \
                    SH.GalaxyID \
                FROM \
                    RefL0100N1504_SubHalo as SH, \
                    RefL0100N1504_Aperture as AP \
                WHERE \
                    SH.SnapNum = 28 \
                    and SH.GalaxyID = AP.GalaxyID \
                    and AP.GalaxyID = 1'

    myData = sql.execute_query(con, myQuery)

    print myData