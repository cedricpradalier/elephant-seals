function R = rpy(roll,pitch,yaw)
    R = rotz(yaw)*roty(pitch)*rotx(roll);
