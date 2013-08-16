#!/usr/bin/python

import sqlite3

conn = sqlite3.connect('elephant-seals.db')
c = conn.cursor()

# Create table (no primary key, rowid is sufficient)
c.execute('''CREATE TABLE IF NOT EXISTS ddt
             (individual INTEGER NOT NULL, date REAL, 
                ax REAL, ay REAL, az REAL, 
                mx REAL, my REAL, mz REAL,
                depth REAL, vel REAL, has_vel INTEGER)''')

c.execute('''CREATE TABLE IF NOT EXISTS raw_estimate
        (id INTEGER NOT NULL, dive_status INTEGER, 
        x REAL, y REAL, z REAL, 
        roll REAL, pitch REAL, yaw REAL,
        velocity REAL,
        FOREIGN KEY(id) REFERENCES ddt(rowid))''')

c.execute('''CREATE TABLE IF NOT EXISTS dives
        (start INTEGER NOT NULL, end INTEGER NOT NULL, 
        FOREIGN KEY(start) REFERENCES ddt(rowid),
        FOREIGN KEY(end) REFERENCES ddt(rowid))''')

c.execute('''CREATE TABLE IF NOT EXISTS dive_properties
        (id INTEGER NOT NULL,
        bscale_x REAL, bscale_y REAL, bscale_z REAL,
        k_depth REAL,
        FOREIGN KEY(id) REFERENCES dives(rowid))''')

c.execute('''CREATE TABLE IF NOT EXISTS orientations
        (id INTEGER NOT NULL, use_quaternions INTEGER, 
        x REAL, y REAL, z REAL, w REAL, 
        FOREIGN KEY(id) REFERENCES ddt(rowid))''')

c.execute('''CREATE TABLE IF NOT EXISTS gps
        (individual INTEGER NOT NULL, date REAL, latitude REAL, longitude REAL,
        dive_start INTEGER,
        easting REAL, northing REAL, zone INTEGER, 
        bx REAL, by REAL, bz REAL,
        FOREIGN KEY(dive_start) REFERENCES dives(rowid))''')

# file generated from octave

f = open("data/preload_mat_all.txt","r")
l = f.readlines()
f.close()
l = [[float(y) for y in x.strip().split()] for x in l]
l_ddt = [tuple([0,x[0]+x[1]] + x[2:10] + [int(x[10])]) for x in l]
l_raw = [tuple([i,int(x[18])] + x[11:18]) for i,x in enumerate(l)]


c.executemany('INSERT INTO ddt VALUES ('+'?,'*(len(l_ddt[0])-1)+'?)', l_ddt)
c.executemany('INSERT INTO raw_estimate VALUES ('+'?,'*(len(l_raw[0])-1)+'?)', l_raw)

# Save (commit) the changes
conn.commit()

f = open("data/gps_mat_all.txt","r")
l = f.readlines()
f.close()
l = [[float(y) for y in x.strip().split()] for x in l]
# [ ts tsf lat lon ind err east north zone bx by bz]
l_gps = [tuple([0,x[0]+x[1]] + x[2:4] + [None] + x[6:]) for x in l]
c.executemany('INSERT INTO gps VALUES ('+'?,'*(len(l_gps[0])-1)+'?)', l_gps)

# Save (commit) the changes
conn.commit()
# query = """
#     select gps.rowid,ddt.rowid,min(abs(ddt.date-gps.date)) 
#     from ddt, gps 
#     where abs(ddt.date-gps.date)<30./(24*60) 
#     group by gps.rowid;
# """
# for g_rowid,d_rowid,_ in c.execute(query):
#     c.execute("UPDATE gps SET 
    

# We can also close the connection if we are done with it.
# Just be sure any changes have been committed or they will be lost.
conn.close()


