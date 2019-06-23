#/usr/local/bin/python3
# -*- coding: utf-8 -*-

# This code transforms a time series (S1(t),S2(t),S3(t)) in a quaternion time series (a(t),b(t),c(t),d(t))
# This code is provided for sake of reproducibility, to complement a paper submitted to a conference on Machine Learning


import numpy as np

def axisangle_to_q(v, theta):
    ''' Convert axis-angle parameters to a unit-norm quaternion '''
    n=np.linalg.norm(v)
    v = v/max(n,1e-16) # tolerance of 1e-16 to avoid division by 0
    x, y, z = v
    w = np.cos(theta/2)
    x = x * np.sin(theta/2)
    y = y * np.sin(theta/2)
    z = z * np.sin(theta/2)
    return w, x, y, z

def quat_rotate_2p(v_init, v_end):
    ''' compute unit-norm quaternion to rotate from v_init to v_end '''
    x2,y2,z2=v_end
    x1,y1,z1=v_init
    # dot product (cos(angle))
    dp=(x1*x2)+(y1*y2)+(z1*z2)
    dp=np.clip(dp,-1, 1)
    # cross product (non-normalized axis)
    x=(y1*z2)-(z1*y2)
    y=(z1*x2)-(x1*z2)
    z=(x1*y2)-(y1*x2)
    # axis-angle to quaternion
    (a,b,c,d) = axisangle_to_q(np.array((x,y,z)),np.arccos(dp))
    return a,b,c,d

def rotate_q(q,v):
    ''' rotate 3D vector using a unit-norm quaternion q '''
    a,b,c,d = q
    v1,v2,v3 = v
    t2 =   a*b
    t3 =   a*c
    t4 =   a*d
    t5 =  -b*b
    t6 =   b*c
    t7 =   b*d
    t8 =  -c*c
    t9 =   c*d
    t10 = -d*d
    v1new = 2*( (t8 + t10)*v1 + (t6 -  t4)*v2 + (t3 + t7)*s3 ) + v1
    v2new = 2*( (t4 +  t6)*v1 + (t5 + t10)*v2 + (t9 - t2)*s3 ) + v2
    v3new = 2*( (t7 -  t3)*v1 + (t2 +  t9)*v2 + (t5 + t8)*s3 ) + v3
    return v1new, v2new, v3new

# open IO files
file = open("input.txt","r")
newfile = open("output.txt","w")

# file format (txt)
# col 1 : id
# col 2 : timestamps
# col 3 : S1 (X)
# col 4 : S2 (Y)
# col 5 : S3 (Z)

k=0        # line number in input file
cnt=0      # number of point in the current time series
notsp=0    # counter to verify how often short path inversion is used
id_old=-2. # initial id

for line in file:
    k=k+1

    # read data in current line
    # id, t, S1, S2, S3
    temp=line.rstrip('\n\r')
    if (k==1):
        newfile.write("id" + '\t' + "t" + '\t' +  "a" + '\t' + "X" + '\t' + "Y" + '\t' + "Z" + '\n')
    if (k>1):
        data = temp.rstrip('\n\r').split("\t")
        id=float(data[0])
        t =float(data[1])
        s1=float(data[2])
        s2=float(data[3])
        s3=float(data[4])

        # make sure that (s1,s2,s3) is on the unit sphere
        n=np.linalg.norm((s1,s2,s3))
        s1=s1/n
        s2=s2/n
        s3=s3/n

    # when a new trajectory starts (with different id) reset to North pole
    # compute quaternion to rotate (s1,s2,s3) to North pole (0,0,1)
    if ((id!=id_old) and (k>1)):
        print("new id=", id)
        cnt=0
        q0,q1,q2,q3 = quat_rotate_2p((s1,s2,s3),(0,0,1))
        s1,s2,s3=rotate_q((q0,q1,q2,q3), (s1,s2,s3))

    # convert current trajectory to quaternion series
    if ((id==id_old) and (k>1)):
        cnt=cnt+1

        # apply the rotation that reset the trajectory to North pole
        s1,s2,s3=rotate_q((q0,q1,q2,q3), (s1,s2,s3))

        # compute the quaternion representing the rotation between the previous
        # point and the current point
        p0,p1,p2,p3 = quat_rotate_2p((s1_old,s2_old,s3_old), (s1,s2,s3))

        # select canonical form
		# other option: condition on sign at position arg max(|p0|,|p1|,|p2|,|p3|)
        if (p0<0):
            p0, p1, p2, p3=-p0, -p1, -p2, -p3

        # shortest path
        if cnt>1:
            dp=np.dot((p0,p1,p2,p3),quat_old)  # = cos(angle)
            if (dp<0):
                notsp=notsp+1
                print(" sp= ", cnt, " ", "notsp = ", notsp)
                p0, p1, p2, p3=-p0, -p1, -p2, -p3

        # write id, t, quaternion
        if (id<10): newfile.write('0000' + str(int(id)) + '\t')
        if ((id<100) and (id>=10)): newfile.write('000' + str(int(id)) + '\t')
        if ((id<1000) and (id>=100)): newfile.write('00' + str(int(id)) + '\t')
        if ((id<10000) and (id>=1000)): newfile.write('0' + str(int(id)) + '\t')
        if ((id<100000) and (id>=10000)): newfile.write(  str(int(id)) + '\t')
        newfile.write(str(t) + '\t')
        newfile.write(str(p0) + '\t')
        newfile.write(str(p1) + '\t')
        newfile.write(str(p2) + '\t')
        newfile.write(str(p3) + '\n')

        # save quaternion
        quat_old=(p0, p1, p2, p3)

    #save data
    if (k>1):
        id_old, s1_old, s2_old, s3_old =id, s1, s2, s3
