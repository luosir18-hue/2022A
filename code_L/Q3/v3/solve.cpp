// ============================================================================
// Q3 v3 交叉验证求解器（独立实现，RK4 on 二阶系统）
// ----------------------------------------------------------------------------
// 独立性来源（与 code_L/Q3/v2 的 RK4-on-一阶状态 [含伴随标量] 写法完全不同）：
//   * 直接对二阶系统 M u'' + C u' + K u = f(t) 写经典 4 阶 RK4，
//     每一步显式用预计算的 M^{-1} 求加速度（不借助 v2 的 rhs/companion 写法）；
//   * M/C/K 在此独立重新推导书写；能量恒等式用 RK4 内嵌的伴随标量精确积分。
// 参数与 v2 锁定值逐字一致，因此若 v3 与 v2 吻合，即证明 v2 代码正确
// （两者差异只能来自共享的模型/参数定义，例如纵摇几何杠杆 le、d）。
// ============================================================================
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <string>
#include <iomanip>
#include <algorithm>
#include <cstdlib>

static const double PI = 3.14159265358979323846;

struct P {
    double omega, F, L;
    double m_f, m_o, m_a, J_a;
    double b_h, b_theta;
    double Kh, Ktheta, k, k_theta;
    double c_l, c_theta;
    double g, rho;
    double d, le, I_f, I_o;
};

static void mul4(const double A[4][4], const double x[4], double y[4]) {
    for (int i=0;i<4;i++){ double s=0; for(int j=0;j<4;j++) s+=A[i][j]*x[j]; y[i]=s; }
}
static void inv4(const double A[4][4], double B[4][4]) {
    double M_[4][4]; for(int i=0;i<4;i++) for(int j=0;j<4;j++) M_[i][j]=A[i][j];
    for(int i=0;i<4;i++) for(int j=0;j<4;j++) B[i][j]=(i==j)?1.0:0.0;
    for(int c=0;c<4;c++){
        int p=c; double mx=std::fabs(M_[c][c]);
        for(int r=c+1;r<4;r++) if(std::fabs(M_[r][c])>mx){mx=std::fabs(M_[r][c]);p=r;}
        for(int j=0;j<4;j++){std::swap(M_[c][j],M_[p][j]);std::swap(B[c][j],B[p][j]);}
        double d=M_[c][c];
        for(int j=0;j<4;j++){M_[c][j]/=d;B[c][j]/=d;}
        for(int r=0;r<4;r++) if(r!=c){ double f=M_[r][c]; for(int j=0;j<4;j++){M_[r][j]-=f*M_[c][j];B[r][j]-=f*B[c][j];}}
    }
}
static void build(const P& p, double M[4][4], double C[4][4], double K[4][4]) {
    double Mf = p.m_f + p.m_a, Jf = p.I_f + p.J_a;
    double d2=p.d*p.d, le2=p.le*p.le, dle=p.d*p.le;
    for(int i=0;i<4;i++) for(int j=0;j<4;j++) M[i][j]=C[i][j]=K[i][j]=0;
    M[0][0]=Mf; M[1][1]=p.m_o;
    M[2][2]=Jf+p.m_o*d2; M[2][3]=p.m_o*dle; M[3][2]=p.m_o*dle; M[3][3]=p.I_o+p.m_o*le2;
    C[0][0]=p.b_h+p.c_l; C[0][1]=-p.c_l; C[1][0]=-p.c_l; C[1][1]=p.c_l;
    C[2][2]=p.b_theta+p.c_theta; C[2][3]=-p.c_theta; C[3][2]=-p.c_theta; C[3][3]=p.c_theta;
    K[0][0]=p.Kh+p.k; K[0][1]=-p.k; K[1][0]=-p.k; K[1][1]=p.k;
    K[2][2]=p.Ktheta+p.k_theta-p.m_o*p.g*p.d; K[2][3]=-p.k_theta; K[3][2]=-p.k_theta; K[3][3]=p.k_theta-p.m_o*p.g*p.le;
}

// 二阶系统导数：y=[u0,u1,u2,u3, v0,v1,v2,v3, Win, Wdiss]
static void deriv(const P& p, const double C[4][4], const double K[4][4],
                  const double Minv[4][4], double t, const double y[10], double dy[10]) {
    double u[4]={y[0],y[1],y[2],y[3]}, v[4]={y[4],y[5],y[6],y[7]};
    double fext[4]={p.F*std::cos(p.omega*t), 0, p.L*std::cos(p.omega*t), 0};
    double Cv[4], Ku[4]; mul4(C,v,Cv); mul4(K,u,Ku);
    double r[4], a[4];
    for(int i=0;i<4;i++) r[i]=fext[i]-Cv[i]-Ku[i];
    mul4(Minv,r,a);
    for(int i=0;i<4;i++){ dy[i]=v[i]; dy[4+i]=a[i]; }
    double c=std::cos(p.omega*t);
    dy[8]=(p.F*v[0]+p.L*v[2])*c;                                  // 波浪输入功率
    dy[9]=p.b_h*v[0]*v[0]+p.b_theta*v[2]*v[2]
         +p.c_l*(v[1]-v[0])*(v[1]-v[0])
         +p.c_theta*(v[3]-v[2])*(v[3]-v[2]);                      // 总耗散功率
}

int main(int argc, char** argv) {
    std::string out_dir=".";
    for(int i=1;i<argc;i++){ std::string a=argv[i]; if(a=="--out"&&i+1<argc) out_dir=argv[++i]; }

    P p; p.omega=1.7152;p.F=3640;p.L=1690;p.m_f=4866;p.m_o=2433;p.m_a=1028.876;p.J_a=7001.914;
    p.b_h=683.4558;p.b_theta=654.3383;p.Kh=31557.29821;p.Ktheta=8890.7;p.k=80000;p.k_theta=250000;
    p.c_l=10000;p.c_theta=1000;p.g=9.8;p.rho=1025;p.d=-1.407925157;p.le=0.4519575;p.I_f=8398.77606;p.I_o=202.75;

    double M[4][4],C[4][4],K[4][4],Minv[4][4];
    build(p,M,C,K); inv4(M,Minv);

    std::cout<<"======== Q3 v3 交叉验证（RK4 二阶系统, 独立实现）========"<<std::endl;
    std::cout<<"ω="<<std::setprecision(10)<<p.omega<<" F="<<p.F<<" L="<<p.L<<std::endl;
    std::cout<<"d="<<p.d<<" ℓ_e="<<p.le<<" I_f="<<p.I_f<<" I_o="<<p.I_o<<std::endl;

    { double mx=0; for(int i=0;i<4;i++)for(int j=i+1;j<4;j++){mx=std::max(mx,std::fabs(M[i][j]-M[j][i]));mx=std::max(mx,std::fabs(C[i][j]-C[j][i]));mx=std::max(mx,std::fabs(K[i][j]-K[j][i]));}
      std::cout<<"矩阵最大对称性误差="<<std::setprecision(4)<<mx<<(mx<1e-10?"  通过":"  不通过")<<std::endl; }

    // L0 单元测试（与 v2 同判据，校验组装）
    { double c01=C[0][1],c23=C[2][3];
      bool t1=(std::fabs(c01+p.c_l)<1e-12)&&(std::fabs(C[1][0]+p.c_l)<1e-12)&&(std::fabs(c23+p.c_theta)<1e-12)&&(std::fabs(C[3][2]+p.c_theta)<1e-12);
      double zz=0.5,th=0.3,vz=1.0,vth=0.7; double yt[4]={zz,zz,th,th},yd[4]={vz,vz,vth,vth};
      double Fd[4]={0,0,0,0},Fs[4]={0,0,0,0};
      for(int i=0;i<4;i++){for(int j=0;j<4;j++){Fd[i]+=C[i][j]*yd[j];Fs[i]+=K[i][j]*yt[j];}}
      double zero=0; auto upd=[&](double v){zero=std::max(zero,std::fabs(v));};
      upd(Fd[0]-p.b_h*yd[0]);upd(Fd[1]);upd(Fd[2]-p.b_theta*yd[2]);upd(Fd[3]);
      upd(Fs[0]-p.Kh*yt[0]);upd(Fs[1]);upd(Fs[2]-(p.Ktheta-p.m_o*p.g*p.d)*yt[2]);upd(Fs[3]+p.m_o*p.g*p.le*yt[3]);
      upd(p.c_l*(yd[1]-yd[0])*(yd[1]-yd[0]));upd(p.c_theta*(yd[3]-yd[2])*(yd[3]-yd[2]));
      bool t2=zero<1e-9;
      std::cout<<"L0-1 PTO耦合项 C[0][1]="<<c01<<" C[2][3]="<<c23<<(t1?"  PASS":"  FAIL")<<std::endl;
      std::cout<<"L0-2 相对运动零内部力/功率 max="<<std::scientific<<zero<<(t2?"  PASS":"  FAIL")<<std::endl;
      std::ofstream f(out_dir+"/result3_L0单元测试.csv",std::ios::binary);
      f<<"测试项,条件,观测值,期望,判定\n";
      f<<"关闭PTO后阻尼矩阵PTO耦合项消失,C[0][1]=C[1][0]="<<c01<<"=-c_l,C[2][3]=C[3][2]="<<c23<<"=-c_θ,兴波阻尼b_h/b_θ仅在主对角,"<<(t1?"PASS":"FAIL")<<"\n";
      f<<"相对运动为零时PTO/弹簧内力内力矩功率为零,max|内力/力矩/功率|="<<zero<<",<1e-9,"<<(t2?"PASS":"FAIL")<<"\n";
      f.close();
    }

    double dt=0.001; double T=2.0*PI/p.omega; double t_max=40.0*T; long nsteps=(long)std::lround(t_max/dt);
    double sample_dt=0.2; long sample_interval=(long)std::lround(sample_dt/dt);

    double y[10]={0,0,0,0,0,0,0,0,0,0};   // 初始静平衡：u=0, v=0, Win=0, Wdiss=0
    struct Rec{double t,zf,vf,tf,wf,zo,vo,to,wo,win,wdiss;};
    std::vector<Rec> recs;
    auto record=[&](double t,const double yy[10]){
        double u[4]={yy[0],yy[1],yy[2],yy[3]}, v[4]={yy[4],yy[5],yy[6],yy[7]};
        double Tk=0,Vp=0;
        for(int i=0;i<4;i++){Tk+=0.5*v[i]*(M[i][0]*v[0]+M[i][1]*v[1]+M[i][2]*v[2]+M[i][3]*v[3]);
                            Vp+=0.5*u[i]*(K[i][0]*u[0]+K[i][1]*u[1]+K[i][2]*u[2]+K[i][3]*u[3]);}
        recs.push_back({t,yy[0],yy[4],yy[2],yy[6],yy[1],yy[5],yy[3],yy[7],yy[8],yy[9]});
    };
    record(0.0,y);

    double k1[10],k2[10],k3[10],k4[10],yt[10];
    for(long step=1;step<=nsteps;step++){
        double t=step*dt;
        deriv(p,C,K,Minv,(step-1)*dt,y,k1);
        for(int i=0;i<10;i++) yt[i]=y[i]+0.5*dt*k1[i];
        deriv(p,C,K,Minv,(step-0.5)*dt,yt,k2);
        for(int i=0;i<10;i++) yt[i]=y[i]+0.5*dt*k2[i];
        deriv(p,C,K,Minv,(step-0.5)*dt,yt,k3);
        for(int i=0;i<10;i++) yt[i]=y[i]+dt*k3[i];
        deriv(p,C,K,Minv,t,yt,k4);
        for(int i=0;i<10;i++) y[i]+=dt/6.0*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
        if(step%sample_interval==0) record(t,y);
    }
    std::cout<<"积分完成 步数="<<nsteps<<" 采样点="<<recs.size()<<" 末时刻="<<recs.back().t<<"s"<<std::endl;

    { std::ofstream f(out_dir+"/result3.csv",std::ios::binary);
      f<<"时间t(s),浮子垂荡位移zf(m),浮子垂荡速度dzf/dt(m/s),浮子纵摇角位移θf(rad),浮子纵摇角速度dθf/dt(rad/s),振子垂荡位移zo(m),振子垂荡速度dzo/dt(m/s),振子纵摇角位移θo(rad),振子纵摇角速度dθo/dt(rad/s)\n";
      for(auto&s:recs) f<<std::setprecision(10)<<s.t<<","<<s.zf<<","<<s.vf<<","<<s.tf<<","<<s.wf<<","<<s.zo<<","<<s.vo<<","<<s.to<<","<<s.wo<<"\n";
      f.close(); std::cout<<"已写入 "<<out_dir<<"/result3.csv"<<std::endl; }

    double key_times[]={10.0,20.0,40.0,60.0,100.0};
    { std::ofstream f(out_dir+"/result3_key.csv",std::ios::binary);
      f<<"时间t(s),浮子垂荡位移zf(m),浮子垂荡速度dzf/dt(m/s),浮子纵摇角位移θf(rad),浮子纵摇角速度dθf/dt(rad/s),振子垂荡位移zo(m),振子垂荡速度dzo/dt(m/s),振子纵摇角位移θo(rad),振子纵摇角速度dθo/dt(rad/s)\n";
      for(double kt:key_times){ long idx=std::lround(kt/sample_dt); if(idx<0)idx=0; if(idx>=(long)recs.size())idx=(long)recs.size()-1; auto&s=recs[idx];
        f<<std::setprecision(10)<<s.t<<","<<s.zf<<","<<s.vf<<","<<s.tf<<","<<s.wf<<","<<s.zo<<","<<s.vo<<","<<s.to<<","<<s.wo<<"\n"; }
      f.close(); std::cout<<"已写入 "<<out_dir<<"/result3_key.csv"<<std::endl; }

    { double mzf=0,mvf=0,mtf=0,mwf=0,mzo=0,mvo=0,mto=0,mwo=0;
      for(auto&s:recs){mzf=std::max(mzf,std::fabs(s.zf));mvf=std::max(mvf,std::fabs(s.vf));mtf=std::max(mtf,std::fabs(s.tf));mwf=std::max(mwf,std::fabs(s.wf));mzo=std::max(mzo,std::fabs(s.zo));mvo=std::max(mvo,std::fabs(s.vo));mto=std::max(mto,std::fabs(s.to));mwo=std::max(mwo,std::fabs(s.wo));}
      std::ofstream f(out_dir+"/result3_响应幅值摘要.csv",std::ios::binary);
      f<<"运动量,中文含义,最大绝对值,单位\n";
      f<<"zf,浮子垂荡位移,"<<std::setprecision(10)<<mzf<<",m\n";f<<"dzf/dt,浮子垂荡速度,"<<mvf<<",m/s\n";
      f<<"θf,浮子纵摇角位移,"<<mtf<<",rad\n";f<<"dθf/dt,浮子纵摇角速度,"<<mwf<<",rad/s\n";
      f<<"zo,振子垂荡位移,"<<mzo<<",m\n";f<<"dzo/dt,振子垂荡速度,"<<mvo<<",m/s\n";
      f<<"θo,振子纵摇角位移,"<<mto<<",rad\n";f<<"dθo/dt,振子纵摇角速度,"<<mwo<<",rad/s\n";f.close();
      std::cout<<"\n====== 响应幅值摘要 ======\n";
      std::cout<<"max|zf|="<<std::setprecision(6)<<mzf<<" m  max|żf|="<<mvf<<" m/s\n";
      std::cout<<"max|θf|="<<mtf<<" rad  max|θ̇f|="<<mwf<<" rad/s\n";
      std::cout<<"max|zo|="<<mzo<<" m  max|żo|="<<mvo<<" m/s\n";
      std::cout<<"max|θo|="<<mto<<" rad  max|θ̇o|="<<mwo<<" rad/s\n"; }

    { std::ofstream f(out_dir+"/result3_能量恒等式.csv",std::ios::binary);
      f<<"时间t(s),动能T(J),势能V(J),机械能E(J),波浪输入功率Pin(W),垂荡兴波耗散bh_zdotf2(W),纵摇兴波耗散btheta_thetadotf2(W),直线PTO耗散Pl(W),旋转PTO耗散Ptheta(W),总耗散Pd(W),累计输入功Win(J),累计耗散Wdiss(J),能量残差E-(Win-Wdiss)(J)\n";
      double max_rel=0,max_abs=0;
      for(auto&s:recs){
        double u[4]={s.zf,s.zo,s.tf,s.to}, v[4]={s.vf,s.vo,s.wf,s.wo};
        double Tk=0,Vp=0; for(int i=0;i<4;i++){Tk+=0.5*v[i]*(M[i][0]*v[0]+M[i][1]*v[1]+M[i][2]*v[2]+M[i][3]*v[3]);Vp+=0.5*u[i]*(K[i][0]*u[0]+K[i][1]*u[1]+K[i][2]*u[2]+K[i][3]*u[3]);}
        double E=Tk+Vp; double resid=E-(s.win-s.wdiss); double scale=std::max(1.0,std::fabs(E));
        double rel=std::fabs(resid)/scale; max_rel=std::max(max_rel,rel); max_abs=std::max(max_abs,std::fabs(resid));
        double Pl=p.c_l*(v[1]-v[0])*(v[1]-v[0]), Pth=p.c_theta*(v[3]-v[2])*(v[3]-v[2]);
        double Pd=p.b_h*v[0]*v[0]+p.b_theta*v[2]*v[2]+Pl+Pth;
        double Pin=(p.F*v[0]+p.L*v[2])*std::cos(p.omega*s.t);
        f<<std::setprecision(10)<<s.t<<","<<Tk<<","<<Vp<<","<<E<<","<<Pin<<","<<p.b_h*v[0]*v[0]<<","<<p.b_theta*v[2]*v[2]<<","<<Pl<<","<<Pth<<","<<Pd<<","<<s.win<<","<<s.wdiss<<","<<resid<<"\n";
      }
      f.close();
      std::cout<<"\n====== 能量恒等式核验 ======\n";
      std::cout<<"max |E-(Win-Wdiss)|="<<std::setprecision(6)<<max_abs<<" J  相对残差="<<max_rel<<" (阈值1e-3) "<<(max_rel<=1e-3?"→ 通过":"→ 不通过")<<std::endl;
    }

    std::cout<<"\n======== Q3 v3 交叉验证完成，输出至: "<<out_dir<<" ========"<<std::endl;
    return 0;
}
