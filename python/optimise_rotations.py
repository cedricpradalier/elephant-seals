#!/usr/bin/python

import math
import os
import sqlite3
import numpy

conn = sqlite3.connect('elephant-seals.db')
c = conn.cursor()

# Create table (no primary key, rowid is sufficient)
# c.execute('''CREATE TABLE IF NOT EXISTS ddt
#              (individual INTEGER, date REAL, 
#                 ax REAL, ay REAL, az REAL, 
#                 mx REAL, my REAL, mz REAL,
#                 depth REAL, vel REAL, has_vel INTEGER)''')
# 
# c.execute('''CREATE TABLE IF NOT EXISTS raw_estimate
#         (id INTEGER, dive_status INTEGER, 
#         x REAL, y REAL, z REAL, 
#         roll REAL, pitch REAL, yaw REAL,
#         velocity REAL,
#         FOREIGN KEY(id) REFERENCES ddt(rowid))''')
# 
# c.execute('''CREATE TABLE IF NOT EXISTS dives
#         (start INTEGER, end INTEGER, 
#         FOREIGN KEY(start) REFERENCES ddt(rowid),
#         FOREIGN KEY(end) REFERENCES ddt(rowid))''')


# file generated from octave
c.execute('''DROP TABLE orientations;''')
c.execute('''DROP TABLE propulsion;''')
c.execute('''CREATE TABLE IF NOT EXISTS orientations
        (id INTEGER NOT NULL, use_quaternions INTEGER, 
        x REAL, y REAL, z REAL, w REAL, 
        FOREIGN KEY(id) REFERENCES ddt(rowid))''')
c.execute('''CREATE TABLE IF NOT EXISTS propulsion
        (id INTEGER NOT NULL, float P, 
        FOREIGN KEY(id) REFERENCES ddt(rowid))''')

c.execute("select rowid,start,end from dives")
dives = c.fetchall()
for (dive_id,start,end) in dives:
    query="""
        select
        ddt.rowid,date,ax,ay,az,mx,my,mz,depth,vel,has_vel,x,y,z,roll,pitch,yaw,velocity,dive_status 
        from ddt left join raw_estimate on raw_estimate.id = ddt.rowid
        where (id not null) and (ddt.rowid>=?) and (ddt.rowid<=?)
        """
    c.execute(query,(start,end))
    l=[list(x) for x in c.fetchall()]
    m=[[math.floor(x[1]),x[1]-math.floor(x[1])]+x[2:] for x in l] + [int(dive_id)];
    m=[" ".join([str(y) for y in x]) for x in m]
    f=open("rot_input.txt","w");
    f.write("\n".join(m))
    f.close()
    print "Launching optimisation for dive [%d , %d] (%d rows)" % (start,end,end-start+1)
    res = os.system("../ceres/cerise/bin/optimise_rotations -input rot_input.txt -output rot_output.txt -num_threads 4 -robustify -use_quaternions -use_local_parameterization")
    print "Optimisation completed: %d " % res
    if res == 0:
        f = open("rot_output.txt","r")
        rot = [l.strip() for l in f.readlines()]
        quat = [[float(x) for x in l.split(" ")] for l in rot if l[0] != '#']
        header= [l[1:].split(" ") for l in rot if l[0]=='#']
        dive_prop = {}
        for h in header: 
            if h[0] == 'UseQuaternion':
                dive_prop['UseQuaternion'] = bool(h[1])
            elif h[0] == 'Bscale':
                dive_prop['Bscale'] = [float(x) for x in h[1:]]
            elif h[0] == 'Kdepth':
                dive_prop['Kdepth'] = float(h[1])
        c.execute("insert into dive_properties VALUES (?,?,?,?,?)",
                (r,dive_prop['Bscale'][0],dive_prop['Bscale'][1],dive_prop['Bscale'][2],dive_prop['Kdepth']))

        if dive_prop['UseQuaternion']:
            inp = [(x[0],q[2],q[3],q[4],q[1]) for (x,q) in zip(l,quat)]
            c.executemany("insert into orientations VALUES (?,1,?,?,?,?)",inp)
            inp = [(x[0],q[5]) for (x,q) in zip(l,quat)]
            c.executemany("insert into propulsion VALUES (?,?)",inp)
        else:
            inp = [(x[0],q[1],q[2],q[3]) for (x,q) in zip(l,quat)]
            c.executemany("insert into orientations VALUES (?,1,?,?,?,0)",inp)
            inp = [(x[0],q[4]) for (x,q) in zip(l,quat)]
            c.executemany("insert into propulsion VALUES (?,?)",inp)
        conn.commit()

# We can also close the connection if we are done with it.
# Just be sure any changes have been committed or they will be lost.
conn.close()


