/* C functions for the test particle module */
/* Incase you forget how to compile the shared lib file
 * gcc -o functions.so -shared -fPIC functions.c
 */



#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define E(c,i,j)  E[c + (i)*3 + (j)*nx*3]
#define B(c,i,j)  B[c + (i)*3 + (j)*nx*3]

#define Ep(c,p,it) Ep[c + 3*(p) + (it)*3*npart]
#define Bp(c,p,it) Bp[c + 3*(p) + (it)*3*npart]

#define pos(c,p,it) pos[c + 3*(p) + (it)*3*npart]
#define vel(c,p,it) vel[c + 3*(p) + (it)*3*npart]



void fieldinterp(float *E,
                 float *B,
                 unsigned int nx, unsigned int ny,
                 double x, double y,
                 double dx, double dy,
                 double xmin, double ymin,
                 double Ei[3], double Bi[3])
{
    double id, jd, di, dj;
    unsigned int i, j;
    double w1, w2, w3, w4;

    id = (x-xmin)/dx;
    jd = (y-ymin)/dy;

    i   = (unsigned int) id;   /* index of the cell */
    j   = (unsigned int) jd;

    di = id - i;              /* position between x and dx*/
    dj = jd - j;              /* position betwee y and dy */


    /* interpolation coefs */
    w1 = (1.0-di)*(1.0-dj);
    w2 = (1.0-di)*(    dj);
    w3 = (    di)*(    dj);
    w4 = (    di)*(1.0-dj);
    
    //printf("Starting Finterp\n");
    /* printf("xy(%lf %lf), Bzi(%lf)\n",x,y,Bi[2]); 
    printf("xy(%lf %lf), Bz(%lf)\n",x,y,B(2,i,j)); */
    //printf("id= %i, jd= %i \n",i,j);
 

    Ei[0] = w1*E(0,i,j) + w2*E(0,i,j+1) + w3*E(0,i+1,j+1) + w4*E(0,i+1,j);
    Ei[1] = w1*E(1,i,j) + w2*E(1,i,j+1) + w3*E(1,i+1,j+1) + w4*E(1,i+1,j);
    Ei[2] = w1*E(2,i,j) + w2*E(2,i,j+1) + w3*E(2,i+1,j+1) + w4*E(2,i+1,j);

    Bi[0] = w1*B(0,i,j) + w2*B(0,i,j+1) + w3*B(0,i+1,j+1) + w4*B(0,i+1,j);
    Bi[1] = w1*B(1,i,j) + w2*B(1,i,j+1) + w3*B(1,i+1,j+1) + w4*B(1,i+1,j);
    Bi[2] = w1*B(2,i,j) + w2*B(2,i,j+1) + w3*B(2,i+1,j+1) + w4*B(2,i+1,j);

    //printf("Done Finterp\n");
    /* printf("xy(%lf %lf), Bzi(%lf)\n",x,y,Bi[2]); 
    printf("xy(%lf %lf), Bz(%lf)\n",x,y,B(2,i,j)); */

}







void accelerate(double r[3], double v[3],
                float *E ,
                float *B,
                double dx, double dy, double xmin,double ymin,
                unsigned int nx, unsigned int ny,
                double dt, double charge, double mass,
                double vnew[3],int direction)
{
    double uplus[3];
    double uminus[3];
    double uprime[3];
    double a, t[3],s[3],t2;
    double Ei[3], Bi[3];
    double c2 = 225.0;
    double gamma;

    fieldinterp(E,
                B,
                nx,ny,
                r[0], r[1],
                dx, dy,
                xmin, ymin,
                Ei, Bi);  /* interpolate the fields at the particle's position*/

    a = 0.5*charge/mass * dt;

    /*B is multiplied by 'direction' to go backward or
     * forward and keep VxB correct.*/

    /*Now we are trying make this relativistic */

    Bi[0] *= direction;
    Bi[1] *= direction;
    Bi[2] *= direction;

    gamma = sqrt(1. - (v[0]*v[0] + v[1]*v[1] + v[2]*v[2])/c2);

    uminus[0] = v[0]/gamma + a*Ei[0];
    uminus[1] = v[1]/gamma + a*Ei[1];
    uminus[2] = v[2]/gamma + a*Ei[2];

    gamma = sqrt(1. + (uminus[0]*uminus[0] + 
                       uminus[1]*uminus[1] +
                       uminus[2]*uminus[2])/c2);


    t[0] = a*Bi[0]/gamma;
    t[1] = a*Bi[1]/gamma;
    t[2] = a*Bi[2]/gamma;

    t2 = t[0]*t[0] + t[1]*t[1] + t[2]*t[2];

    s[0] = 2*t[0]/(1. + t2);
    s[1] = 2*t[1]/(1. + t2);
    s[2] = 2*t[2]/(1. + t2);

    uprime[0] = uminus[0] + (uminus[1]*t[2] - uminus[2]*t[1]);
    uprime[1] = uminus[1] + (uminus[2]*t[0] - uminus[0]*t[2]);
    uprime[2] = uminus[2] + (uminus[0]*t[1] - uminus[1]*t[0]);

    uplus[0] = uminus[0] + (uprime[1]*s[2] - uprime[2]*s[1]);
    uplus[1] = uminus[1] + (uprime[2]*s[0] - uprime[0]*s[2]);
    uplus[2] = uminus[2] + (uprime[0]*s[1] - uprime[1]*s[0]);

    vnew[0] = uplus[0] + a*Ei[0];
    vnew[1] = uplus[1] + a*Ei[1];
    vnew[2] = uplus[2] + a*Ei[2];

    gamma = sqrt(1. + (vnew[0]*vnew[0] + 
                       vnew[1]*vnew[1] + 
                       vnew[2]*vnew[2])/c2);

    vnew[0] = vnew[0]/gamma;
    vnew[1] = vnew[1]/gamma;
    vnew[2] = vnew[2]/gamma;

}





void move(double r[3], double v[3], double dt, double rnew[3])
{
    unsigned int c;

    for (c=0; c<3; c++)
    {
        rnew[c] = r[c] + v[c]*dt;
    }
}






void moveall(double *pos,
             double *vel,
             float *E,
             float *B,
             unsigned int nx, unsigned int ny,
             double dx, double dy,
             double xmin, double ymin,
             double dt,
             double charge, double mass,
             unsigned int nt, unsigned int it0,
             unsigned int npart)
{

    int it, p;
    double r[3], v[3], rnew[3], vnew[3];
    double xmax, ymax;
    xmax = xmin+dx*(nx-2);
    ymax = ymin+dy*(ny-2);

/*
    printf("movallC : npart(%d), nt(%d),ito(%d),charge(%lf),mass(%lf),dt(%lf)xymin(%lf,%lf),nxy(%d,%d)\n",
            npart,nt,it0,charge,mass,dt,xmin,ymin,nx,ny);

    printf("selection velocity and position : %lf %lf %lf %lf %lf %lf\n",
                    pos(0,0,it0),pos(1,0,it0),pos(2,0,it0),
                    vel(0,0,it0),vel(1,0,it0),vel(2,0,it0));

    printf("Size of the B array is %lf\n",sizeof(B[0]));
    printf("We just want to see a Bz of 1 any where! %lf, %lf, %lf\n",B[0][0][0],B[1][0][0],B[2][0][0]);
*/
    /* reverse the velocity at t=it0 to go backward */
    for (p=0; p<npart; p++)
    {
        vel(0,p,it0) *= -1;
        vel(1,p,it0) *= -1;
        vel(2,p,it0) *= -1;
    }


    for (it = it0-1; it>=0; it--) /* backward integration */
    {
        for (p=0; p < npart; p++)
        {
            /* position and velocity used to 
             * accelerate the particle
             * we are calculating positions and velocities
             * at time 'it' and backward, so 'previous' positions and
             * velocities are those at it+1
             *
             * store the new velocity backward in the array
            */

            /* positions and velocity at t = it+1 */
            r[0] = pos(0,p,it+1);
            r[1] = pos(1,p,it+1);
            r[2] = pos(2,p,it+1);

            v[0] = vel(0,p,it+1);
            v[1] = vel(1,p,it+1);
            v[2] = vel(2,p,it+1);

            /* get v(it) and store it */
            accelerate(r,v,E,B,dx,dy,xmin,ymin,nx,ny,dt,charge,mass,vnew,-1);

            vel(0,p,it) = vnew[0];
            vel(1,p,it) = vnew[1];
            vel(2,p,it) = vnew[2];

            /* get r(it) based on v(it) */
            move(r, vnew, dt, rnew);

            /* double periodic wrapping */
            if (rnew[0] > xmax){ rnew[0] = xmin;}
            if (rnew[0] < xmin){ rnew[0] = xmax;}
            if (rnew[1] > ymax){ rnew[1] = ymin;}
            if (rnew[1] < ymin){ rnew[1] = ymax;}
            /* open periodict wrapping */
            //if (rnew[1] > ymax){ rnew[1] = pos(1,p,it-1); vel(1,p,it) = -1.*vel(1,p,it);}
            //if (rnew[1] > ymin){ rnew[1] = pos(1,p,it-1); vel(1,p,it) = -1.*vel(1,p,it);}

            pos(0,p,it) = rnew[0];
            pos(1,p,it) = rnew[1];
            pos(2,p,it) = rnew[2];
        }
    }/* end of backward loop */

    /*
    printf("selection velocity and position : %lf %lf %lf %lf %lf %lf\n",
                    pos(0,0,it0),pos(1,0,it0),pos(2,0,it0),
                    vel(0,0,it0),vel(1,0,it0),vel(2,0,it0));
*/


    /* reverse the velocities from t=0 to t=it0 included 
     * to go forward from it0 and to have a good sign for previous v's values*/
    for (it = 0; it<= it0; it++)
    {
        for (p=0; p<npart; p++)
        {
            vel(0,p,it) *= -1;
            vel(1,p,it) *= -1;
            vel(2,p,it) *= -1;
        }
    }

    /*printf("now go forward...\n");*/

    for (it = it0+1; it < nt; it++) /* forward integration */
    {
        for (p=0; p < npart; p++)
        {
            /* so we now know the positions and velocities at it-1
             * and want to calculate those at it*/
            /* positions and velocity at t = it+1 */
            r[0] = pos(0,p,it-1);
            r[1] = pos(1,p,it-1);
            r[2] = pos(2,p,it-1);

            v[0] = vel(0,p,it-1);
            v[1] = vel(1,p,it-1);
            v[2] = vel(2,p,it-1);

            /* get v(it) and store it */
            accelerate(r,v,E,B,dx,dy,xmin,ymin,nx,ny,dt,charge,mass,vnew,1);

            vel(0,p,it) = vnew[0];
            vel(1,p,it) = vnew[1];
            vel(2,p,it) = vnew[2];

            /* get r(it) based on v(it) */
            move(r, vnew, dt, rnew);

            /* double periodic wrapping */
            if (rnew[0] > xmax){ rnew[0] = xmin;}
            if (rnew[0] < xmin){ rnew[0] = xmax;}
            if (rnew[1] > ymax){ rnew[1] = ymin;}
            if (rnew[1] < ymin){ rnew[1] = ymax;}
            /* open periodict wrapping */
            //if (rnew[1] > ymax){ rnew[1] = pos(1,p,it-1); vel(1,p,it) = -1.*vel(1,p,it);}
            //if (rnew[1] > ymin){ rnew[1] = pos(1,p,it-1); vel(1,p,it) = -1.*vel(1,p,it);}

            pos(0,p,it) = rnew[0];
            pos(1,p,it) = rnew[1];
            pos(2,p,it) = rnew[2];
        }
    }/* end of backward loop */

    /*
    printf("selection velocity and position : %lf %lf %lf %lf %lf %lf\n",
                    pos(0,0,it0),pos(1,0,it0),pos(2,0,it0),
                    vel(0,0,it0),vel(1,0,it0),vel(2,0,it0));
*/

}






void pfields(double *pos,
             float *E,
             float *B,
             unsigned int nx, unsigned int ny,
             unsigned int npart,
             unsigned int nt,
             double xmin, double ymin,
             double dx, double dy,
             double *Ep,
             double *Bp)
{
    unsigned int p, it;
    double Ei[3], Bi[3];

    for (p=0;  p < npart; p++)
    {
        for (it=0; it < nt; it++)
        {
            fieldinterp(E,
                        B,
                        nx, ny,
                        pos(0,p,it), pos(1,p,it),
                        dx, dy,
                        xmin, ymin,
                        Ei, Bi);  /* interpolate the fields at the particle's position*/

            /* put that in the arrays */
            Ep(0,p,it) = Ei[0];
            Ep(1,p,it) = Ei[1];
            Ep(2,p,it) = Ei[2];
            Bp(0,p,it) = Bi[0];
            Bp(1,p,it) = Bi[1];
            Bp(2,p,it) = Bi[2];

        }
    }
}




















