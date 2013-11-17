#include <stdlib.h>
#include <ros/ros.h>

#include <tf/LinearMath/Vector3.h>
#include <tf/LinearMath/Matrix3x3.h>
#include <tf/LinearMath/Transform.h>

double U() {
    return (-1.0 + (2.0*random())/RAND_MAX);
}

int main(int argc,char *argv[]) {
    ros::init(argc,argv,"es_sim");
    ros::Time::init();
    ros::Time t0 = ros::Time::now();
    ros::Rate sim_rate(100);
    tf::Vector3 G(0.0,0.0,-10.0);
    tf::Vector3 B(cos(M_PI/3),sin(M_PI/3),+5.0);

    tf::Vector3 r(10.0,10.0,5.0);

    printf("#t X0 X1 X2  V0 V1 V2  roll pitch yaw  a0 a1 a2  b0 b1 b2\n");
    while (ros::ok()) {
        double t = (ros::Time::now() - t0).toSec();
        tf::Vector3 X(r[0]*cos(t*M_PI/5.3),
                r[1]*sin(t*M_PI/4.1),
                r[2]*sin(t*M_PI/6.2));
        // printf("X %.2f %.2f %.2f\n",X[0],X[1],X[2]);
        // printf("V %.2f %.2f %.2f\n",V[0],V[1],V[2]);
        tf::Vector3 V = tf::Vector3(-r[0]*M_PI/5.3*sin(t*M_PI/5.3),
                r[1]*M_PI/4.1*cos(t*M_PI/4.1),
                r[2]*M_PI/6.2*cos(t*M_PI/6.2));
        tf::Vector3 vx = V.normalized();
        tf::Vector3 vy = (-X).normalized();
        tf::Vector3 vz = (vx.cross(vy)).normalized();
        vy = vz.cross(vx).normalized();
        tf::Matrix3x3 R(vx[0],vy[0],vz[0],
                vx[1],vy[1],vz[1],
                vx[2],vy[2],vz[2]);
        // printf("R %.2f %.2f %.2f\n",R[0][0],R[0][1],R[0][2]);
        // printf("  %.2f %.2f %.2f\n",R[1][0],R[1][1],R[1][2]);
        // printf("  %.2f %.2f %.2f\n",R[2][0],R[2][1],R[2][2]);
        double roll,pitch,yaw;
        R.getRPY(roll,pitch,yaw);
        tf::Vector3 a = R.transpose() * G;
        tf::Vector3 b = R.transpose() * B;
        // assert(fabs(a.norm() - 10) < 2.0);

        printf("%f  %e %e %e  %e %e %e  %e %e %e  %e %e %e  %e %e %e\n",
                t,X[0],
                X[1],X[2], V[0],V[1],V[2], 
                roll,pitch,yaw, 
                a[0]+0.1*U(),a[1]+0.1*U(),a[2]+0.1*U(), 
                b[0]+0.1*U(),b[1]+0.1*U(),b[2]+0.1*U());
        sim_rate.sleep();
    }


    return 0;
}

