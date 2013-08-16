#!/usr/bin/python

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

query = """
select id,individual,date,dive_status 
from ddt left join raw_estimate on raw_estimate.id = ddt.rowid
where id not null
order by individual,date;
"""
c.execute(query)
l = c.fetchall();
rows = [r for (r,i,t,d) in l]
indiv = numpy.array([i for (r,i,t,d) in l]);
dindiv = [0] + list(numpy.diff(indiv));
# List of change of individual in the sequence
l_new_record = [(i,r) for ((r,_,_,_),(i,delta)) in zip (l,enumerate(dindiv)) if delta!=0]
ds = numpy.array([d for (r,i,t,d) in l]);
dds = [0] + list(numpy.diff(ds));
# List of departure from the surface
d_start = [(i,r) for ((r,_,_,d),(i,delta)) in zip(l,enumerate(dds)) if (d==1) and (delta==1)]
# List of arrival at the surface
d_end = [(i,r) for ((r,_,_,d),(i,delta)) in zip(l,enumerate(dds)) if (d==0) and (delta==-3)]

# symbols, for reference
events = {0:'start', 1:'new individual', 2:'dive starts', 3: 'dive ends', 4:'end'}

# Create a list of symbols
evlist = [(0,0), (len(l)-1,4)] + [(i,1) for (i,r) in l_new_record] \
        + [(i,2) for (i,r) in d_start] + [(i,3) for (i,r) in d_end]
evlist.sort()
# Removes duplicates
dev = [1] + list(numpy.diff(numpy.array([e for (i,e) in evlist])))
rep = set([i for (i,d) in enumerate(dev) if d == 0]);
evlist = [e for (i,e) in enumerate(evlist) if i not in rep]

# Now we can extract the dive limits

assert(len(evlist)>=2) # Otherwise, it does not make sense to extract dives

#  evlist[0] is a start event
(i_first,t_first) = evlist[1];
if t_first == 2:
    # We first see a dive, which means we were on the surface initially
    # Replace the event with the beginning of the sequence
    evlist[1] = (0,2);
if t_first == 3:
    # We first see a surface event, which means we were diving before.
    evlist = [evlist[0], (0,2)] + evlist[1:]

(i_last, t_last) = evlist[-2];
n = len(l)-1;
if t_last == 2:
    # The last thing we see is a dive starting, so we add an end at the end
    evlist = evlist[:-1] + [(n,3), (n,4)]
if t_last == 3:
    # The last thing we see is the end of a dive, so we finish on the surface
    evlist[-2] = (n,3)
        
el=[]
for (i,t) in evlist:
    if t==1:
        el += [(i,3),(i,2)]
    else:
        el.append((i,t))
evlist = el


# Now the list is [(0,0) (x,2), (x,3), .....(x,2), (x,3), (n,4)]
assert(len(evlist)%2 == 0)
r_end = [r for (r,t) in evlist if t == 3]
r_start = [r for (r,t) in evlist if t == 2]
surf_time = zip(r_end[:-1],r_start[1:])
dive_limits = [(s+e)/2 for (s,e) in surf_time]
dive_starts = [rows[r_start[0]]] + [rows[dl] for dl in dive_limits]
dive_ends = [rows[dl-1] for dl in dive_limits] + [rows[r_end[-1]]]

c.executemany('INSERT INTO dives VALUES (?,?)', zip(dive_starts,dive_ends))

# query = """
#     select gps.rowid,ddt.rowid,min(abs(ddt.date-gps.date)) 
#     from ddt, gps 
#     where abs(ddt.date-gps.date)<30./(24*60) 
#     group by gps.rowid;
# """
# for g_rowid,d_rowid,_ in c.execute(query):
#     c.execute("UPDATE gps SET 
    
conn.commit()

# We can also close the connection if we are done with it.
# Just be sure any changes have been committed or they will be lost.
conn.close()


