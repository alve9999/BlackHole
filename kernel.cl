#define A2 (a*a)
#define R (1.0 + sqrt(1.0-A2)+0.00001)
#define DT 0.01
#define Z1  (1.0 + pow(1.0-A2,(1.0/3.0))*(pow(1.0+a,(1.0/3.0))+pow(1.0-a,(1.0/3.0))))
#define Z2  pow(3*A2+Z1*Z1,1/2)
#define SR  (3.0 + Z2 - sqrt((3.0-Z1)*(3.0+Z1+2.0*Z2)))
#define ER  (20.0)
#define HEIGHT 1360
#define WIDTH 2048


typedef struct {
    double rad;
    double theta;
    double phi;
    double pr;
    double ptheta;
}state;

static state func(state y,double K, double L,double a) {
    state dxdy;
    double sintheta = sin(y.theta);
    double costheta = cos(y.theta);
    if(fabs(sintheta) < 1e-8)
    {
        sintheta = sign(sintheta)*1e-8;
    }
    if(fabs(costheta) < 1e-8)
    {
        costheta = sign(costheta)*1e-8;
    }
    double r2 = y.rad * y.rad;

    double delta = r2 - 2.0*y.rad + A2;
    double sigma = r2 - A2*costheta*costheta;

    dxdy.rad = -y.pr*delta/sigma;
    dxdy.theta = -y.ptheta/sigma;
    dxdy.phi = -((2.0*a*y.rad + (sigma-2.0*y.rad)*L/(sintheta*sintheta))/(sigma*delta));
    dxdy.pr = -(((y.rad-1.0)*(-K) + 2.0*y.rad*(r2+A2) - 2.0*a*L)/(sigma*delta) - 2.0*y.pr*y.pr*(y.rad-1.0)/sigma);
    dxdy.ptheta = -(sintheta*costheta*(L*L/(pow(sintheta,4))-A2)/sigma);
    return dxdy;
}

static state increment(state y,state dxdy,double dt,double K,double L){
    y.rad += dt*dxdy.rad;
    y.theta += dt*dxdy.theta;
    y.phi += dt*dxdy.phi;
    y.pr += dt*dxdy.pr;
    y.ptheta += dt*dxdy.ptheta;
    return y;
}

typedef struct{
    state out;
    state err;
}CKRes;

static CKRes CashKarp(state y,double h,double K,double L,double a){
    const double firstOrder[] = {1.0/5.0,3.0/40.0,3.0/10.0,-11.0/54.0,1631.0/55296.0,37.0/378.0};
    const double secondOrder[] = {9.0/40.0,-9.0/10.0,5.0/2.0,175.0/512.0};
    const double thirdOrder[] = {6.0/5.0,-70.0/27.0,575.0/13824.0,250.0/621.0};
    const double fourthOrder[] = {35.0/27.0,44275.0/110592.0,125.0/594.0};
    const double fifthOrder[] = {253.0/4096.0,512.0/1771.0};
    const double err[] = {((37.0/378.0)-(2825.0/27648.0)),((250.0/621.0)-(18575.0/48384.0)),((125.0/594.0)-(13525.0/55296.0)),-(277.0/14336.0),((512.0/1771.0)-0.25)};
    state k0,k1,k2,k3,k4;
    state dxdy = func(y,K,L,a);
    CKRes res;

    //first order
    {
        double hdxdyrad = h * dxdy.rad;
        double rad = y.rad;
        k0.rad = rad + firstOrder[0]*hdxdyrad;
        k1.rad = rad + firstOrder[1]*hdxdyrad;
        k2.rad = rad + firstOrder[2]*hdxdyrad;
        k3.rad = rad + firstOrder[3]*hdxdyrad;
        k4.rad = rad + firstOrder[4]*hdxdyrad;
        res.out.rad= rad + firstOrder[5]*hdxdyrad;
        res.err.rad = err[0]*hdxdyrad;
    }
    {
        double hdxdytheta = h * dxdy.theta;
        double theta = y.theta;
        k0.theta = theta + firstOrder[0]*hdxdytheta;
        k1.theta = theta + firstOrder[1]*hdxdytheta;
        k2.theta = theta + firstOrder[2]*hdxdytheta;
        k3.theta = theta + firstOrder[3]*hdxdytheta;
        k4.theta = theta + firstOrder[4]*hdxdytheta;
        res.out.theta= theta + firstOrder[5]*hdxdytheta;
        res.err.theta = err[0]*hdxdytheta;
    }
    {
        double hdxdyphi = h * dxdy.phi;
        double phi = y.phi;
        k0.phi = phi + firstOrder[0]*hdxdyphi;
        k1.phi = phi + firstOrder[1]*hdxdyphi;
        k2.phi = phi + firstOrder[2]*hdxdyphi;
        k3.phi = phi + firstOrder[3]*hdxdyphi;
        k4.phi = phi + firstOrder[4]*hdxdyphi;
        res.out.phi= phi + firstOrder[5]*hdxdyphi;
        res.err.phi = err[0]*hdxdyphi;
    }
    {
        double hdxdypr = h * dxdy.pr;
        double pr = y.pr;
        k0.pr = pr + firstOrder[0]*hdxdypr;
        k1.pr = pr + firstOrder[1]*hdxdypr;
        k2.pr = pr + firstOrder[2]*hdxdypr;
        k3.pr = pr + firstOrder[3]*hdxdypr;
        k4.pr = pr + firstOrder[4]*hdxdypr;
        res.out.pr= pr + firstOrder[5]*hdxdypr;
        res.err.pr = err[0]*hdxdypr;
    }
    {
        double hdxdyptheta = h * dxdy.ptheta;
        double ptheta = y.ptheta;
        k0.ptheta = ptheta + firstOrder[0]*hdxdyptheta;
        k1.ptheta = ptheta + firstOrder[1]*hdxdyptheta;
        k2.ptheta = ptheta + firstOrder[2]*hdxdyptheta;
        k3.ptheta = ptheta + firstOrder[3]*hdxdyptheta;
        k4.ptheta = ptheta + firstOrder[4]*hdxdyptheta;
        res.out.ptheta= ptheta + firstOrder[5]*hdxdyptheta;
        res.err.ptheta = err[0]*hdxdyptheta;
    }

    //second order
    dxdy = func(k0,K,L,a);
    {
        double hdxdyrad = h * dxdy.rad;
        k1.rad += secondOrder[0]*hdxdyrad;
        k2.rad += secondOrder[1]*hdxdyrad;
        k3.rad += secondOrder[2]*hdxdyrad;
        k4.rad += secondOrder[3]*hdxdyrad;
    }
    {
        double hdxdytheta = h * dxdy.theta;
        k1.theta += secondOrder[0]*hdxdytheta;
        k2.theta += secondOrder[1]*hdxdytheta;
        k3.theta += secondOrder[2]*hdxdytheta;
        k4.theta += secondOrder[3]*hdxdytheta;
    }
    {
        double hdxdyphi = h * dxdy.phi;
        k1.phi += secondOrder[0]*hdxdyphi;
        k2.phi += secondOrder[1]*hdxdyphi;
        k3.phi += secondOrder[2]*hdxdyphi;
        k4.phi += secondOrder[3]*hdxdyphi;
    }
    {
        double hdxdypr = h * dxdy.pr;
        k1.pr += secondOrder[0]*hdxdypr;
        k2.pr += secondOrder[1]*hdxdypr;
        k3.pr += secondOrder[2]*hdxdypr;
        k4.pr += secondOrder[3]*hdxdypr;
    }
    {
        double hdxdyptheta = h * dxdy.ptheta;
        k1.ptheta += secondOrder[0]*hdxdyptheta;
        k2.ptheta += secondOrder[1]*hdxdyptheta;
        k3.ptheta += secondOrder[2]*hdxdyptheta;
        k4.ptheta += secondOrder[3]*hdxdyptheta;
    }

    //third order
    dxdy = func(k1,K,L,a);
    {
        double hdxdyrad = h * dxdy.rad;
        k2.rad += thirdOrder[0]*hdxdyrad;
        k3.rad += thirdOrder[1]*hdxdyrad;
        k4.rad += thirdOrder[2]*hdxdyrad;
        res.out.rad += thirdOrder[3]*hdxdyrad;
        res.err.rad += err[1]*hdxdyrad;
    }
    {
        double hdxdytheta = h * dxdy.theta;
        k2.theta += thirdOrder[0]*hdxdytheta;
        k3.theta += thirdOrder[1]*hdxdytheta;
        k4.theta += thirdOrder[2]*hdxdytheta;
        res.out.theta += thirdOrder[3]*hdxdytheta;
        res.err.theta += err[1]*hdxdytheta;
    }
    {
        double hdxdyphi = h * dxdy.phi;
        k2.phi += thirdOrder[0]*hdxdyphi;
        k3.phi += thirdOrder[1]*hdxdyphi;
        k4.phi += thirdOrder[2]*hdxdyphi;
        res.out.phi += thirdOrder[3]*hdxdyphi;
        res.err.phi += err[1]*hdxdyphi;
    }
    {
        double hdxdypr = h * dxdy.pr;
        k2.pr += thirdOrder[0]*hdxdypr;
        k3.pr += thirdOrder[1]*hdxdypr;
        k4.pr += thirdOrder[2]*hdxdypr;
        res.out.pr += thirdOrder[3]*hdxdypr;
        res.err.pr += err[1]*hdxdypr;
    }
    {
        double hdxdyptheta = h * dxdy.ptheta;
        k2.ptheta += thirdOrder[0]*hdxdyptheta;
        k3.ptheta += thirdOrder[1]*hdxdyptheta;
        k4.ptheta += thirdOrder[2]*hdxdyptheta;
        res.out.ptheta += thirdOrder[3]*hdxdyptheta;
        res.err.ptheta += err[1]*hdxdyptheta;
    }

    //fourth order
    dxdy = func(k2,K,L,a);
    {
        double hdxdyrad = h * dxdy.rad;
        k3.rad += fourthOrder[0]*hdxdyrad;
        k4.rad += fourthOrder[1]*hdxdyrad;
        res.out.rad += fourthOrder[2]*hdxdyrad;
        res.err.rad += err[2]*hdxdyrad;
    }
    {
        double hdxdytheta = h * dxdy.theta;
        k3.theta += fourthOrder[0]*hdxdytheta;
        k4.theta += fourthOrder[1]*hdxdytheta;
        res.out.theta += fourthOrder[2]*hdxdytheta;
        res.err.theta += err[2]*hdxdytheta;
    }
    {
        double hdxdyphi = h * dxdy.phi;
        k3.phi += fourthOrder[0]*hdxdyphi;
        k4.phi += fourthOrder[1]*hdxdyphi;
        res.out.phi += fourthOrder[2]*hdxdyphi;
        res.err.phi += err[2]*hdxdyphi;
    }
    {
        double hdxdypr = h * dxdy.pr;
        k3.pr += fourthOrder[0]*hdxdypr;
        k4.pr += fourthOrder[1]*hdxdypr;
        res.out.pr += fourthOrder[2]*hdxdypr;
        res.err.pr += err[2]*hdxdypr;
    }
    {
        double hdxdyptheta = h * dxdy.ptheta;
        k3.ptheta += fourthOrder[0]*hdxdyptheta;
        k4.ptheta += fourthOrder[1]*hdxdyptheta;
        res.out.ptheta += fourthOrder[2]*hdxdyptheta;
        res.err.ptheta += err[2]*hdxdyptheta;
    }

    //fifth order
    dxdy = func(k3,K,L,a);
    {
        double hdxdyrad = h * dxdy.rad;
        k4.rad += fifthOrder[0]*hdxdyrad;
        res.err.rad += err[3]*hdxdyrad;
    }
    {
        double hdxdytheta = h * dxdy.theta;
        k4.theta += fifthOrder[0]*hdxdytheta;
        res.err.theta += err[3]*hdxdytheta;
    }
    {
        double hdxdyphi = h * dxdy.phi;
        k4.phi += fifthOrder[0]*hdxdyphi;
        res.err.phi += err[3]*hdxdyphi;
    }
    {
        double hdxdypr = h * dxdy.pr;
        k4.pr += fifthOrder[0]*hdxdypr;
        res.err.pr += err[3]*hdxdypr;
    }
    {
        double hdxdyptheta = h * dxdy.ptheta;
        k4.ptheta += fifthOrder[0]*hdxdyptheta;
        res.err.ptheta += err[3]*hdxdyptheta;
    }

    //final out
    dxdy = func(k4,K,L,a);
    {
        double hdxdyrad = h * dxdy.rad;
        res.out.rad += fifthOrder[1]*hdxdyrad;
        res.err.rad += err[4]*hdxdyrad;
    }
    {
        double hdxdytheta = h * dxdy.theta;
        res.out.theta += fifthOrder[1]*hdxdytheta;
        res.err.theta += err[4]*hdxdytheta;
    }
    {
        double hdxdyphi = h * dxdy.phi;
        res.out.phi += fifthOrder[1]*hdxdyphi;
        res.err.phi += err[4]*hdxdyphi;
    }
    {
        double hdxdypr = h * dxdy.pr;
        res.out.pr += fifthOrder[1]*hdxdypr;
        res.err.pr += err[4]*hdxdypr;
    }
    {
        double hdxdyptheta = h * dxdy.ptheta;
        res.out.ptheta += fifthOrder[1]*hdxdyptheta;
        res.err.ptheta += err[4]*hdxdyptheta;
    }
    return res;
}
static double calculate_next_stepsize(state err,double h){
    double avg_err = (err.rad+err.theta+err.pr+err.ptheta+err.phi)/5;
    return h*avg_err;
}

static state RungeKutta4(state y,double K,double L,double a){
    state k1,k2,k3,k4;
    state y_temp;

    k1 = func(y,K,L,a);

    y_temp = increment(y,k1,DT/2.0,K,L);

    k2 = func(y_temp,K,L,a);

    y_temp = increment(y,k2,DT/2.0,K,L);

    k3 = func(y_temp,K,L,a);

    y_temp = increment(y,k3,DT,K,L);

    k4 = func(y_temp,K,L,a);

    y.rad = y.rad + DT / 6.0 * (k1.rad + 2.0 * k2.rad + 2.0 * k3.rad + k4.rad);
    y.theta = y.theta + DT / 6.0 * (k1.theta + 2.0 * k2.theta + 2.0 * k3.theta + k4.theta);
    y.phi = y.phi + DT / 6.0 * (k1.phi + 2.0 * k2.phi + 2.0 * k3.phi + k4.phi);
    y.pr = y.pr + DT / 6.0 * (k1.pr + 2.0 * k2.pr + 2.0 * k3.pr + k4.pr);
    y.ptheta = y.ptheta + DT / 6.0 * (k1.ptheta + 2.0 * k2.ptheta + 2.0 * k3.ptheta + k4.ptheta);
    return y;
}

float lcg(ulong x) {
    const ulong a = 166452ul;
    const ulong c = 1013904ul;
    return (float)((a * x + c)%100000);
}

float generate_uniform_random(ulong seed) {
    ulong x = get_global_id(0);
    ulong y = get_global_id(1);

    ulong state = seed + x ^ (y + 0x9e3779b9 + (x << 6) + (x >> 2));

    return (lcg(state)/100000.0)*2.0-1.0;
}

__kernel void compute_ray(__write_only image2d_t result,double r0,double theta0,double a, int seed){
    const float4 temp[] = {(float4)(1.0, 0.2196078431372549, 0.0,1.0),(float4)(1.0, 0.3254901960784314, 0.0,1.0),(float4)(1.0, 0.396078431372549, 0.0,1.0),(float4)(1.0, 0.45098039215686275, 0.0,1.0),(float4)(1.0, 0.49411764705882355, 0.0,1.0),(float4)(1.0, 0.5372549019607843, 0.07058823529411765,1.0),(float4)(1.0, 0.5764705882352941, 0.17254901960784313,1.0),(float4)(1.0, 0.615686274509804, 0.24705882352941178,1.0),(float4)(1.0, 0.6470588235294118, 0.30980392156862746,1.0),(float4)(1.0, 0.6784313725490196, 0.3686274509803922,1.0),(float4)(1.0, 0.7058823529411765, 0.4196078431372549,1.0),(float4)(1.0, 0.7333333333333333, 0.47058823529411764,1.0),(float4)(1.0, 0.7568627450980392, 0.5176470588235295,1.0),(float4)(1.0, 0.7803921568627451, 0.5607843137254902,1.0),(float4)(1.0, 0.8, 0.6,1.0),(float4)(1.0, 0.8196078431372549, 0.6392156862745098,1.0),(float4)(1.0, 0.8352941176470589, 0.6784313725490196,1.0),(float4)(1.0, 0.8509803921568627, 0.7137254901960784,1.0),(float4)(1.0, 0.8666666666666667, 0.7450980392156863,1.0),(float4)(1.0, 0.8823529411764706, 0.7764705882352941,1.0),(float4)(1.0, 0.8941176470588236, 0.807843137254902,1.0),(float4)(1.0, 0.9098039215686274, 0.8352941176470589,1.0),(float4)(1.0, 0.9215686274509803, 0.8627450980392157,1.0),(float4)(1.0, 0.9333333333333333, 0.8901960784313725,1.0),(float4)(1.0, 0.9411764705882353, 0.9137254901960784,1.0),(float4)(1.0, 0.9529411764705882, 0.9372549019607843,1.0),(float4)(1.0, 0.9607843137254902, 0.9607843137254902,1.0),(float4)(1.0, 0.9725490196078431, 0.984313725490196,1.0),(float4)(0.996078431372549, 0.9764705882352941, 1.0,1.0),(float4)(0.9764705882352941, 0.9647058823529412, 1.0,1.0),(float4)(0.9607843137254902, 0.9529411764705882, 1.0,1.0),(float4)(0.9411764705882353, 0.9450980392156862, 1.0,1.0),(float4)(0.9294117647058824, 0.9372549019607843, 1.0,1.0),(float4)(0.9137254901960784, 0.9294117647058824, 1.0,1.0),(float4)(0.9019607843137255, 0.9215686274509803, 1.0,1.0),(float4)(0.8901960784313725, 0.9137254901960784, 1.0,1.0),(float4)(0.8784313725490196, 0.9058823529411765, 1.0,1.0),(float4)(0.8666666666666667, 0.9019607843137255, 1.0,1.0),(float4)(0.8549019607843137, 0.8941176470588236, 1.0,1.0),(float4)(0.8470588235294118, 0.8901960784313725, 1.0,1.0),(float4)(0.8392156862745098, 0.8823529411764706, 1.0,1.0),(float4)(0.8274509803921568, 0.8784313725490196, 1.0,1.0),(float4)(0.8196078431372549, 0.8745098039215686, 1.0,1.0),(float4)(0.8117647058823529, 0.8666666666666667, 1.0,1.0),(float4)(0.807843137254902, 0.8627450980392157, 1.0,1.0),(float4)(0.8, 0.8588235294117647, 1.0,1.0),(float4)(0.792156862745098, 0.8549019607843137, 1.0,1.0),(float4)(0.788235294117647, 0.8509803921568627, 1.0,1.0),(float4)(0.7803921568627451, 0.8470588235294118, 1.0,1.0),(float4)(0.7764705882352941, 0.8470588235294118, 1.0,1.0),(float4)(0.7686274509803922, 0.8431372549019608, 1.0,1.0),(float4)(0.7647058823529411, 0.8392156862745098, 1.0,1.0),(float4)(0.7607843137254902, 0.8352941176470589, 1.0,1.0),(float4)(0.7568627450980392, 0.8313725490196079, 1.0,1.0),(float4)(0.7529411764705882, 0.8313725490196079, 1.0,1.0),(float4)(0.7490196078431373, 0.8274509803921568, 1.0,1.0),(float4)(0.7450980392156863, 0.8235294117647058, 1.0,1.0),(float4)(0.7411764705882353, 0.8235294117647058, 1.0,1.0),(float4)(0.7372549019607844, 0.8196078431372549, 1.0,1.0),(float4)(0.7333333333333333, 0.8196078431372549, 1.0,1.0),(float4)(0.7294117647058823, 0.8156862745098039, 1.0,1.0),(float4)(0.7254901960784313, 0.8156862745098039, 1.0,1.0),(float4)(0.7215686274509804, 0.8117647058823529, 1.0,1.0),(float4)(0.7176470588235294, 0.8117647058823529, 1.0,1.0),(float4)(0.7176470588235294, 0.807843137254902, 1.0,1.0),(float4)(0.7137254901960784, 0.807843137254902, 1.0,1.0),(float4)(0.7098039215686275, 0.803921568627451, 1.0,1.0),(float4)(0.7098039215686275, 0.803921568627451, 1.0,1.0),(float4)(0.7058823529411765, 0.8, 1.0,1.0),(float4)(0.7019607843137254, 0.8, 1.0,1.0),(float4)(0.7019607843137254, 0.8, 1.0,1.0),(float4)(0.6980392156862745, 0.796078431372549, 1.0,1.0),(float4)(0.6980392156862745, 0.796078431372549, 1.0,1.0),(float4)(0.6941176470588235, 0.792156862745098, 1.0,1.0),(float4)(0.6941176470588235, 0.792156862745098, 1.0,1.0),(float4)(0.6901960784313725, 0.792156862745098, 1.0,1.0),(float4)(0.6862745098039216, 0.788235294117647, 1.0,1.0),(float4)(0.6862745098039216, 0.788235294117647, 1.0,1.0),(float4)(0.6862745098039216, 0.788235294117647, 1.0,1.0),(float4)(0.6823529411764706, 0.788235294117647, 1.0,1.0),(float4)(0.6823529411764706, 0.7843137254901961, 1.0,1.0),(float4)(0.6784313725490196, 0.7843137254901961, 1.0,1.0),(float4)(0.6784313725490196, 0.7843137254901961, 1.0,1.0),(float4)(0.6745098039215687, 0.7803921568627451, 1.0,1.0),(float4)(0.6745098039215687, 0.7803921568627451, 1.0,1.0),(float4)(0.6745098039215687, 0.7803921568627451, 1.0,1.0),(float4)(0.6705882352941176, 0.7803921568627451, 1.0,1.0),(float4)(0.6705882352941176, 0.7764705882352941, 1.0,1.0),(float4)(0.6666666666666666, 0.7764705882352941, 1.0,1.0),(float4)(0.6666666666666666, 0.7764705882352941, 1.0,1.0),(float4)(0.6666666666666666, 0.7764705882352941, 1.0,1.0),(float4)(0.6627450980392157, 0.7764705882352941, 1.0,1.0),(float4)(0.6627450980392157, 0.7725490196078432, 1.0,1.0),(float4)(0.6627450980392157, 0.7725490196078432, 1.0,1.0),(float4)(0.6627450980392157, 0.7725490196078432, 1.0,1.0),(float4)(0.6588235294117647, 0.7725490196078432, 1.0,1.0),(float4)(0.6588235294117647, 0.7725490196078432, 1.0,1.0),(float4)(0.6588235294117647, 0.7686274509803922, 1.0,1.0),(float4)(0.6549019607843137, 0.7686274509803922, 1.0,1.0),(float4)(0.6549019607843137, 0.7686274509803922, 1.0,1.0),(float4)(0.6549019607843137, 0.7686274509803922, 1.0,1.0),(float4)(0.6549019607843137, 0.7686274509803922, 1.0,1.0),(float4)(0.6509803921568628, 0.7647058823529411, 1.0,1.0),(float4)(0.6509803921568628, 0.7647058823529411, 1.0,1.0),(float4)(0.6509803921568628, 0.7647058823529411, 1.0,1.0),(float4)(0.6509803921568628, 0.7647058823529411, 1.0,1.0),(float4)(0.6470588235294118, 0.7647058823529411, 1.0,1.0),(float4)(0.6470588235294118, 0.7647058823529411, 1.0,1.0),(float4)(0.6470588235294118, 0.7647058823529411, 1.0,1.0),(float4)(0.6470588235294118, 0.7607843137254902, 1.0,1.0),(float4)(0.6431372549019608, 0.7607843137254902, 1.0,1.0),(float4)(0.6431372549019608, 0.7607843137254902, 1.0,1.0),(float4)(0.6431372549019608, 0.7607843137254902, 1.0,1.0),(float4)(0.6431372549019608, 0.7607843137254902, 1.0,1.0),(float4)(0.6431372549019608, 0.7607843137254902, 1.0,1.0),(float4)(0.6392156862745098, 0.7607843137254902, 1.0,1.0),(float4)(0.6392156862745098, 0.7568627450980392, 1.0,1.0),(float4)(0.6392156862745098, 0.7568627450980392, 1.0,1.0),(float4)(0.6392156862745098, 0.7568627450980392, 1.0,1.0),(float4)(0.6392156862745098, 0.7568627450980392, 1.0,1.0),(float4)(0.6392156862745098, 0.7568627450980392, 1.0,1.0),(float4)(0.6352941176470588, 0.7568627450980392, 1.0,1.0),(float4)(0.6352941176470588, 0.7568627450980392, 1.0,1.0),(float4)(0.6352941176470588, 0.7568627450980392, 1.0,1.0),(float4)(0.6352941176470588, 0.7568627450980392, 1.0,1.0),(float4)(0.6352941176470588, 0.7529411764705882, 1.0,1.0),(float4)(0.6352941176470588, 0.7529411764705882, 1.0,1.0),(float4)(0.6313725490196078, 0.7529411764705882, 1.0,1.0),(float4)(0.6313725490196078, 0.7529411764705882, 1.0,1.0),(float4)(0.6313725490196078, 0.7529411764705882, 1.0,1.0),(float4)(0.6313725490196078, 0.7529411764705882, 1.0,1.0),(float4)(0.6313725490196078, 0.7529411764705882, 1.0,1.0),(float4)(0.6313725490196078, 0.7529411764705882, 1.0,1.0),(float4)(0.6313725490196078, 0.7529411764705882, 1.0,1.0),(float4)(0.6274509803921569, 0.7529411764705882, 1.0,1.0),(float4)(0.6274509803921569, 0.7490196078431373, 1.0,1.0),(float4)(0.6274509803921569, 0.7490196078431373, 1.0,1.0),(float4)(0.6274509803921569, 0.7490196078431373, 1.0,1.0),(float4)(0.6274509803921569, 0.7490196078431373, 1.0,1.0),(float4)(0.6274509803921569, 0.7490196078431373, 1.0,1.0),(float4)(0.6274509803921569, 0.7490196078431373, 1.0,1.0),(float4)(0.6274509803921569, 0.7490196078431373, 1.0,1.0),(float4)(0.6235294117647059, 0.7490196078431373, 1.0,1.0),(float4)(0.6235294117647059, 0.7490196078431373, 1.0,1.0),(float4)(0.6235294117647059, 0.7490196078431373, 1.0,1.0)};

    double K = 0.0;
    double L = 0.0;

    int2 coords = (int2)(get_global_id(0),get_global_id(1));

    double nx = ((double)get_global_id(0)*2.0 + 1*generate_uniform_random(seed) - ((double)WIDTH)) / (double)HEIGHT;
    double ny = ((double)get_global_id(1)*2.0 + 1*generate_uniform_random(2*seed) - ((double)HEIGHT)) / (double)HEIGHT;
    state y;
    y.rad = r0;
    y.theta= theta0;
    y.phi=0;
    double t = 0.0;

    double sintheta = sin(y.theta);
    double costheta = cos(y.theta);
    double cos2 = costheta * costheta;
    double sin2 = sintheta * sintheta;
    double dr0 = cos(ny) * cos(ny);
    double dtheta0 = sin(ny)/r0;
    double r2 = y.rad * y.rad;
    double sigma = r2 + A2 * cos2;
    double delta = r2 - 2.0 * y.rad + A2;
    double s1 = sigma - 2.0 * y.rad;
    double dphi0 = cos(ny) * sin(nx) / (r0*sintheta);
    double epsilon2 = s1 * (dr0 * dr0 / delta + dtheta0 * dtheta0) + delta * sin2 * dphi0 * dphi0;
    double energy = sqrt(epsilon2);

    y.pr = dr0 * sigma / (delta*energy);
    y.ptheta = dtheta0 * sigma/energy;
    L = ((sigma * delta * dphi0 - 2.0 * a * y.rad * energy) * sin2 / s1) / energy;
    K = y.ptheta * y.ptheta + A2 * sin2 + L*L/sin2;
    double h = 0.01;
    while (1) {
        CKRes res = CashKarp(y,h,K,L,a);
        y=res.out;
        h = calculate_next_stepsize(res.err,h);
        if((y.rad>r0) || (y.rad<R) || t>15000.0){
            write_imagef(result, coords, (float4)(0.0, 0.0, 0.0, 0.0));
            return;
        }
        if((y.rad <= ER) && (y.rad >= SR) && (0.01>fabs(y.theta-M_PI/2.0))||(0.01>fabs(y.theta+M_PI/2.0))) {

            double heat = pow((y.rad*10-20),-3.0/4.0);
            double normalised_temp=clamp(heat/pow(SR,-3.0/4.0),0.0,1.0);

            write_imagef(result, coords, temp[(int)(normalised_temp*145.0)]);
            return;
        }

    t += DT;
    }
}
